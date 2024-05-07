#include <stdio.h>
#include "esp_log.h"
#include "esp_adc/adc_continuous.h"
#include "esp_timer.h"
#include "freertos/FreeRTOS.h"
#include "freertos/task.h"
#include "FFT.h"
#include <string.h>
#include <string>

#include "U8g2lib.h"
extern "C" {
#include "u8g2_esp32_hal.h"
}

#define ADC_FRAME_SIZE 2048

#define THRESHOLD_AMPLITUDE 2000

#define AVG_VALUES_NUM 5

#define SIGNAL_THRESHOLD_FRONT 4
#define SIGNAL_THRESHOLD_BACK 2

#define DOT_DUR_MIN 10
#define DOT_DUR_MAX 25

#define DASH_DUR_MIN 26
#define DASH_DUR_MAX 100

char data[100];
uint32_t dataCouter = 0; 

#define SAMPLE_FREQ 100000
static const char* TAG = "Main";

U8G2 display = U8G2();

const uint16_t samples = 4096; //This value MUST ALWAYS be a power of 2
const float signalFrequency = 1000;
const float samplingFrequency = SAMPLE_FREQ;

float vReal[2][samples];
float vImag[2][samples];

uint16_t fillingIndex = 0;
uint8_t fillingBuffer = 0;
uint8_t processingBuffer = 0;

std::string answer = "";

bool dataAvailable = false;


FFT fft;


adc_channel_t channel[1] = {ADC_CHANNEL_0};

adc_continuous_handle_t adc_handle = NULL;

volatile uint16_t convCount = 0;
volatile uint16_t ovfCount = 0; 

static bool IRAM_ATTR s_conv_done_cb(adc_continuous_handle_t handle, const adc_continuous_evt_data_t *edata, void *user_data) {
    convCount++;
    return false;
}

static bool IRAM_ATTR s_pool_ovf_cb(adc_continuous_handle_t handle, const adc_continuous_evt_data_t *edata, void *user_data) {
    ovfCount++;
    return false;
}


void configure_adc() {
    adc_continuous_handle_cfg_t adc_config = {
        .max_store_buf_size = 8192,
        .conv_frame_size = ADC_FRAME_SIZE,
    };

    ESP_ERROR_CHECK(adc_continuous_new_handle(&adc_config, &adc_handle));


    adc_continuous_config_t dig_cfg = {
        .pattern_num = 1,
        .sample_freq_hz = SAMPLE_FREQ,
        .conv_mode = ADC_CONV_SINGLE_UNIT_1,
        .format = ADC_DIGI_OUTPUT_FORMAT_TYPE1,
    };

    adc_digi_pattern_config_t adc_pattern[1] = {0};

    adc_pattern[0].atten = ADC_ATTEN_DB_11;
    adc_pattern[0].channel = channel[0] & 0x7;
    adc_pattern[0].unit = ADC_UNIT_1;
    adc_pattern[0].bit_width = SOC_ADC_DIGI_MAX_BITWIDTH;

    dig_cfg.adc_pattern = adc_pattern;

    ESP_ERROR_CHECK(adc_continuous_config(adc_handle, &dig_cfg));

    adc_continuous_evt_cbs_t cbs = {
        .on_conv_done = s_conv_done_cb,
        .on_pool_ovf = s_pool_ovf_cb
    };
    ESP_ERROR_CHECK(adc_continuous_register_event_callbacks(adc_handle, &cbs, NULL));
}

adc_digi_output_data_t results[ADC_FRAME_SIZE / sizeof(adc_digi_output_data_t)] = {0};

void adcStats(void* args) {
    while (true) {
        ESP_LOGI(TAG, "ADC conv count = %d, overflow count = %d", convCount, ovfCount);
        vTaskDelay(pdMS_TO_TICKS(100));
    }
}

void morseDetect(uint16_t sampleCount) {
    if (sampleCount >= DOT_DUR_MIN && sampleCount <= DOT_DUR_MAX) {
        ESP_LOGI(TAG, "DOT");
        strcat(data, ".");
    }
    else if (sampleCount >= DASH_DUR_MIN && sampleCount <= DASH_DUR_MAX) {
        ESP_LOGI(TAG, "DASH");
        strcat(data, "-");
    }
}

void decodeMorze(){
    if (strcmp(data, ".-") == 0)
        answer += "а";
    if (strcmp(data, "-...") == 0)
        answer += "б";
    if (strcmp(data, ".--") == 0)
        answer += "в";
    if (strcmp(data, "--.") == 0)
        answer += "г";
    if (strcmp(data, "-..") == 0)
        answer += "д";
    if (strcmp(data, ".") == 0)
        answer += "е";
    if (strcmp(data, "...-") == 0)
        answer += "ж";
    if (strcmp(data, "--..") == 0)
        answer += "з";
    if (strcmp(data, "..") == 0)
        answer += "и";
    if (strcmp(data, ".---") == 0)
        answer += "й";
    if (strcmp(data, "-.-") == 0)
        answer += "к";
    if (strcmp(data, ".-..") == 0)
        answer += "л";
    if (strcmp(data, "--") == 0)
        answer += "м";
    if (strcmp(data, "-.") == 0)
        answer += "н";
    if (strcmp(data, "---") == 0)
        answer += "о";
    if (strcmp(data, ".--.") == 0)
        answer += "п";
    if (strcmp(data, ".-.") == 0)
        answer += "р";
    if (strcmp(data, "...") == 0)
        answer += "с";
    if (strcmp(data, "-") == 0)
        answer += "т";
    if (strcmp(data, "..-") == 0)
        answer += "у";
    if (strcmp(data, "..-.") == 0)
        answer += "ф";
    if (strcmp(data, "....") == 0)
        answer += "х";
    if (strcmp(data, "-.-.") == 0)
        answer += "ц";
    if (strcmp(data, "---.") == 0)
        answer += "ч";
    if (strcmp(data, "----") == 0)
        answer += "ш";
    if (strcmp(data, "--.-") == 0)
        answer += "щ";
    if (strcmp(data, ".--.-.") == 0)
        answer += "ъ";
    if (strcmp(data, "-.--") == 0)
        answer += "ы";
    if (strcmp(data, "-..-") == 0)
        answer += "ь";
    if (strcmp(data, "..-..") == 0)
        answer += "э";
    if (strcmp(data, "--..") == 0)
        answer += "ю";
    if (strcmp(data, ".-.-") == 0)
        answer += "я";
}

void pauseHandler(){   
    decodeMorze();
    ESP_LOGI(TAG, "%s", answer.c_str());
    data[0] = '\0';
}


void processAmplitude(float amplitude) {
    static int8_t avgValue = 0;
    static int32_t sampleCounter = 0;
    if (amplitude >= THRESHOLD_AMPLITUDE && avgValue < AVG_VALUES_NUM) {
        avgValue++;
    } 
    else if (amplitude < THRESHOLD_AMPLITUDE && avgValue > 0) {
        avgValue--;
    }
    if (sampleCounter > 0 && avgValue <= SIGNAL_THRESHOLD_BACK) {
        //ESP_LOGI(TAG, "Signal fall detected, duration = %ld", sampleCounter);
        morseDetect(sampleCounter);
        sampleCounter = 0;
    }
    else if (sampleCounter > 0 || avgValue >= SIGNAL_THRESHOLD_FRONT) {
        if(sampleCounter < 0) {
           // ESP_LOGI(TAG, "Pause detected, duration = %ld", -sampleCounter);
            if(sampleCounter < -50){
                pauseHandler();
            }
            sampleCounter = 1;
        } else {
            sampleCounter++;
        }
    } 
    else if (sampleCounter <= 0 && avgValue <= SIGNAL_THRESHOLD_BACK) {
        if(sampleCounter < -100){
            pauseHandler();
            sampleCounter = -1;
        }
        sampleCounter--;
    } 

}

void fftTask(void* args) {
    while (true) {
        if (dataAvailable) {
            dataAvailable = false;
            fft.dcRemoval(vReal[processingBuffer], samples);
            fft.windowing(vReal[processingBuffer], samples, FFTWindow::Hamming, FFTDirection::Forward);
            fft.compute(vReal[processingBuffer], vImag[processingBuffer], samples, FFTDirection::Forward);
            fft.complexToMagnitude(vReal[processingBuffer], vImag[processingBuffer], samples);

            processAmplitude(vReal[processingBuffer][50]);

            //ESP_LOGI(TAG, "1000 Hz amplitude = %lu", (uint32_t)vReal[processingBuffer][50]);
        }
        vTaskDelay(pdMS_TO_TICKS(5));
    }
}

void displayTask(void *arg) {
    while (true) {
        display.drawUTF8(0, 10, "Распознано:");
        display.drawUTF8(0, 20, answer.c_str());
        display.sendBuffer();
        vTaskDelay(pdMS_TO_TICKS(500));
    }
    
}

extern "C" void app_main(void) {
    u8g2_esp32_hal_t u8g2_esp32_hal = U8G2_ESP32_HAL_DEFAULT;
    u8g2_esp32_hal.bus.spi.clk = GPIO_NUM_25;
    u8g2_esp32_hal.bus.spi.mosi = GPIO_NUM_33;
    u8g2_esp32_hal.bus.spi.cs = GPIO_NUM_32,
    u8g2_esp32_hal.reset = GPIO_NUM_26;
    u8g2_esp32_hal.positiveCS = true;
    u8g2_esp32_hal_init(u8g2_esp32_hal);

    
    display.getU8g2();

    u8g2_Setup_st7920_s_128x64_f( display.getU8g2(), U8G2_R0, u8g2_esp32_spi_byte_cb, u8g2_esp32_gpio_and_delay_cb);
    ESP_LOGI(TAG, "Using ST7920 LCD");


    display.initDisplay();
    display.setPowerSave(0);
    display.setContrast(45);

    display.enableUTF8Print();
    display.setFont(u8g2_font_6x12_t_cyrillic);



    configure_adc();

    // xTaskCreatePinnedToCore(adcStats, "ADC stats", 4096, NULL, 2, NULL, 1);
    xTaskCreatePinnedToCore(displayTask, "Display task", 8192, NULL, 2, NULL, 1);
    xTaskCreatePinnedToCore(fftTask, "FFT task", 4096, NULL, 2, NULL, 1);
    
    ESP_ERROR_CHECK(adc_continuous_start(adc_handle));

    while (1) {
        uint32_t ret_num = 0;
        esp_err_t ret = adc_continuous_read(adc_handle, (uint8_t*)results, ADC_FRAME_SIZE, &ret_num, 0);
        if (ret == ESP_OK) {
            for (int i = 0; i < ret_num / sizeof(adc_digi_output_data_t); i++) {
                vReal[fillingBuffer][fillingIndex] = results[i].type1.data;
                vImag[fillingBuffer][fillingIndex] = 0.0;
                fillingIndex++;
                if (fillingIndex == samples) {
                    fillingIndex = 0;
                    processingBuffer = fillingBuffer;
                    fillingBuffer = (fillingBuffer + 1) % 2;
                    dataAvailable = true;
                }
            }
            vTaskDelay(1);
        } else if (ret == ESP_ERR_TIMEOUT) {
            // ESP_LOGI(TAG, "ADC data timeout");
            vTaskDelay(1);
        }
        
    }

    ESP_LOGI(TAG, "app_main finished");

}