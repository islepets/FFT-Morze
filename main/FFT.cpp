#include "FFT.h"
#include "cmath"

FFT::FFT() {}

FFT::FFT(float* vReal, float* vImag, uint_fast16_t samples,
    float samplingFrequency, bool windowingFactors)
    : _samples(samples), _samplingFrequency(samplingFrequency), _vImag(vImag),
    _vReal(vReal) {
    if (windowingFactors) {
        _precompiledWindowingFactors = new float[samples / 2];
    }
    _power = exponent(samples);
#ifdef FFT_SPEED_OVER_PRECISION
    _oneOverSamples = 1.0 / samples;
#endif
}

FFT::~FFT(void) {
    // Destructor
    if (_precompiledWindowingFactors) {
        delete[] _precompiledWindowingFactors;
    }
}

void FFT::complexToMagnitude(void) const {
    complexToMagnitude(this->_vReal, this->_vImag, this->_samples);
}

void FFT::complexToMagnitude(float* vReal, float* vImag,
    uint_fast16_t samples) const {
    // vM is half the size of vReal and vImag
    for (uint_fast16_t i = 0; i < samples; i++) {
        vReal[i] = sqrt_internal(sq(vReal[i]) + sq(vImag[i]));
    }
}

void FFT::compute(FFTDirection dir) const {
    compute(this->_vReal, this->_vImag, this->_samples, exponent(this->_samples),
        dir);
}

void FFT::compute(float* vReal, float* vImag, uint_fast16_t samples,
    FFTDirection dir) const {
    compute(vReal, vImag, samples, exponent(samples), dir);
}

// Computes in-place complex-to-complex FFT
void FFT::compute(float* vReal, float* vImag, uint_fast16_t samples,
    uint_fast8_t power, FFTDirection dir) const {
#ifdef FFT_SPEED_OVER_PRECISION
    T oneOverSamples = this->_oneOverSamples;
    if (!this->_oneOverSamples)
        oneOverSamples = 1.0 / samples;
#endif
    // Reverse bits
    uint_fast16_t j = 0;
    for (uint_fast16_t i = 0; i < (samples - 1); i++) {
        if (i < j) {
            swap(&vReal[i], &vReal[j]);
            if (dir == FFTDirection::Reverse)
                swap(&vImag[i], &vImag[j]);
        }
        uint_fast16_t k = (samples >> 1);

        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j += k;
    }
    // Compute the FFT
    float c1 = -1.0;
    float c2 = 0.0;
    uint_fast16_t l2 = 1;
    for (uint_fast8_t l = 0; (l < power); l++) {
        uint_fast16_t l1 = l2;
        l2 <<= 1;
        float u1 = 1.0;
        float u2 = 0.0;
        for (j = 0; j < l1; j++) {
            for (uint_fast16_t i = j; i < samples; i += l2) {
                uint_fast16_t i1 = i + l1;
                float t1 = u1 * vReal[i1] - u2 * vImag[i1];
                float t2 = u1 * vImag[i1] + u2 * vReal[i1];
                vReal[i1] = vReal[i] - t1;
                vImag[i1] = vImag[i] - t2;
                vReal[i] += t1;
                vImag[i] += t2;
            }
            float z = ((u1 * c1) - (u2 * c2));
            u2 = ((u1 * c2) + (u2 * c1));
            u1 = z;
        }

#if defined(__AVR__) && defined(USE_AVR_PROGMEM)
        c2 = pgm_read_float_near(&(_c2[l]));
        c1 = pgm_read_float_near(&(_c1[l]));
#else
        float cTemp = 0.5 * c1;
        c2 = sqrt_internal(0.5 - cTemp);
        c1 = sqrt_internal(0.5 + cTemp);
#endif

        if (dir == FFTDirection::Forward) {
            c2 = -c2;
        }
    }
    // Scaling for reverse transform
    if (dir == FFTDirection::Reverse) {
        for (uint_fast16_t i = 0; i < samples; i++) {
#ifdef FFT_SPEED_OVER_PRECISION
            vReal[i] *= oneOverSamples;
            vImag[i] *= oneOverSamples;
#else
            vReal[i] /= samples;
            vImag[i] /= samples;
#endif
        }
    }
}

void FFT::dcRemoval(void) const {
    dcRemoval(this->_vReal, this->_samples);
}

void FFT::dcRemoval(float* vData, uint_fast16_t samples) const {
    // calculate the mean of vData
    float  mean = 0;
    for (uint_fast16_t i = 0; i < samples; i++) {
        mean += vData[i];
    }
    mean /= samples;
    // Subtract the mean from vData
    for (uint_fast16_t i = 0; i < samples; i++) {
        vData[i] -= mean;
    }
}

float FFT::majorPeak(void) const {
    return majorPeak(this->_vReal, this->_samples, this->_samplingFrequency);
}

void FFT::majorPeak(float* f, float* v) const {
    majorPeak(this->_vReal, this->_samples, this->_samplingFrequency, f, v);
}

float FFT::majorPeak(float* vData, uint_fast16_t samples,
    float samplingFrequency) const {
    float  frequency;
    majorPeak(vData, samples, samplingFrequency, &frequency, nullptr);
    return frequency;
}

void FFT::majorPeak(float* vData, uint_fast16_t samples,
    float samplingFrequency, float* frequency,
    float* magnitude) const {
    float maxY = 0;
    uint_fast16_t IndexOfMaxY = 0;
    findMaxY(vData, (samples >> 1) + 1, &maxY, &IndexOfMaxY);

    float delta = 0.5 * ((vData[IndexOfMaxY - 1] - vData[IndexOfMaxY + 1]) /
        (vData[IndexOfMaxY - 1] - (2.0 * vData[IndexOfMaxY]) +
            vData[IndexOfMaxY + 1]));
    float interpolatedX = ((IndexOfMaxY + delta) * samplingFrequency) / (samples - 1);
    if (IndexOfMaxY == (samples >> 1)) // To improve calculation on edge values
        interpolatedX = ((IndexOfMaxY + delta) * samplingFrequency) / (samples);
    // returned value: interpolated frequency peak apex
    *frequency = interpolatedX;
    if (magnitude != nullptr) {
#if defined(ESP8266) || defined(ESP32)
        *magnitude = fabs(vData[IndexOfMaxY - 1] - (2.0 * vData[IndexOfMaxY]) +
            vData[IndexOfMaxY + 1]);
#else
        * magnitude = std::abs(vData[IndexOfMaxY - 1] - (2.0 * vData[IndexOfMaxY]) +
            vData[IndexOfMaxY + 1]);
#endif
    }
}

float FFT::majorPeakParabola(void) const {
    float freq = 0;
    majorPeakParabola(this->_vReal, this->_samples, this->_samplingFrequency,
        &freq, nullptr);
    return freq;
}

void FFT::majorPeakParabola(float* frequency, float* magnitude) const {
    majorPeakParabola(this->_vReal, this->_samples, this->_samplingFrequency,
        frequency, magnitude);
}

float FFT::majorPeakParabola(float* vData, uint_fast16_t samples,
    float samplingFrequency) const {
    float freq = 0;
    majorPeakParabola(vData, samples, samplingFrequency, &freq, nullptr);
    return freq;
}

void FFT::majorPeakParabola(float* vData, uint_fast16_t samples,
    float samplingFrequency, float* frequency,
    float* magnitude) const {
    float maxY = 0;
    uint_fast16_t IndexOfMaxY = 0;
    findMaxY(vData, (samples >> 1) + 1, &maxY, &IndexOfMaxY);

    *frequency = 0;
    if (IndexOfMaxY > 0) {
        // Assume the three points to be on a parabola
        float a, b, c;
        parabola(IndexOfMaxY - 1, vData[IndexOfMaxY - 1], IndexOfMaxY,
            vData[IndexOfMaxY], IndexOfMaxY + 1, vData[IndexOfMaxY + 1], &a,
            &b, &c);

        // Peak is at the middle of the parabola
        float x = -b / (2 * a);

        // And magnitude is at the extrema of the parabola if you want It...
        if (magnitude != nullptr) {
            *magnitude = a * x * x + b * x + c;
        }

        // Convert to frequency
        *frequency = (x * samplingFrequency) / samples;
    }
}

uint8_t FFT::revision(void) {
    return (FFT_LIB_REV);
}

// Replace the data array pointers
void FFT::setArrays(float* vReal, float* vImag, uint_fast16_t samples) {
    _vReal = vReal;
    _vImag = vImag;
    if (samples) {
        _samples = samples;
#ifdef FFT_SPEED_OVER_PRECISION
        _oneOverSamples = 1.0 / samples;
#endif
        if (_precompiledWindowingFactors) {
            delete[] _precompiledWindowingFactors;
        }
        _precompiledWindowingFactors = new float[samples / 2];
    }
}

void FFT::windowing(FFTWindow windowType, FFTDirection dir,
    bool withCompensation) {
    // The windowing function is the same, precompiled values can be used, and
    // precompiled values exist
    if (this->_precompiledWindowingFactors && this->_isPrecompiled &&
        this->_windowFunction == windowType &&
        this->_precompiledWithCompensation == withCompensation) {
        windowing(this->_vReal, this->_samples, FFTWindow::Precompiled, dir,
            this->_precompiledWindowingFactors, withCompensation);
        // Precompiled values must be generated. Either the function changed or the
        // precompiled values don't exist
    }
    else if (this->_precompiledWindowingFactors) {
        windowing(this->_vReal, this->_samples, windowType, dir,
            this->_precompiledWindowingFactors, withCompensation);
        this->_isPrecompiled = true;
        this->_precompiledWithCompensation = withCompensation;
        this->_windowFunction = windowType;
        // Don't care about precompiled windowing values
    }
    else {
        windowing(this->_vReal, this->_samples, windowType, dir, nullptr,
            withCompensation);
    }
}

void FFT::windowing(float* vData, uint_fast16_t samples,
    FFTWindow windowType, FFTDirection dir,
    float* windowingFactors, bool withCompensation) {
    // Weighing factors are computed once before multiple use of FFT
    // The weighing function is symmetric; half the weighs are recorded
    if (windowingFactors != nullptr && windowType == FFTWindow::Precompiled) {
        for (uint_fast16_t i = 0; i < (samples >> 1); i++) {
            if (dir == FFTDirection::Forward) {
                vData[i] *= windowingFactors[i];
                vData[samples - (i + 1)] *= windowingFactors[i];
            }
            else {
#ifdef FFT_SPEED_OVER_PRECISION
                T inverse = 1.0 / windowingFactors[i];
                vData[i] *= inverse;
                vData[samples - (i + 1)] *= inverse;
#else
                vData[i] /= windowingFactors[i];
                vData[samples - (i + 1)] /= windowingFactors[i];
#endif
            }
        }
    }
    else {
        float samplesMinusOne = (float(samples) - 1.0);
        float compensationFactor;
        if (withCompensation) {
            compensationFactor =
                _WindowCompensationFactors[static_cast<uint_fast8_t>(windowType)];
        }
        for (uint_fast16_t i = 0; i < (samples >> 1); i++) {
            float indexMinusOne = float(i);
            float ratio = (indexMinusOne / samplesMinusOne);
            float weighingFactor = 1.0;
            // Compute and record weighting factor
            switch (windowType) {
            case FFTWindow::Hamming: // hamming
                weighingFactor = 0.54 - (0.46 * cos(twoPi * ratio));
                break;
            case FFTWindow::Hann: // hann
                weighingFactor = 0.54 * (1.0 - cos(twoPi * ratio));
                break;
            case FFTWindow::Triangle: // triangle (Bartlett)
#if defined(ESP8266) || defined(ESP32)
                weighingFactor =
                    1.0 - ((2.0 * fabs(indexMinusOne - (samplesMinusOne / 2.0))) /
                        samplesMinusOne);
#else
                weighingFactor =
                    1.0 - ((2.0 * std::abs(indexMinusOne - (samplesMinusOne / 2.0))) /
                        samplesMinusOne);
#endif
                break;
            case FFTWindow::Nuttall: // nuttall
                weighingFactor = 0.355768 - (0.487396 * (cos(twoPi * ratio))) +
                    (0.144232 * (cos(fourPi * ratio))) -
                    (0.012604 * (cos(sixPi * ratio)));
                break;
            case FFTWindow::Blackman: // blackman
                weighingFactor = 0.42323 - (0.49755 * (cos(twoPi * ratio))) +
                    (0.07922 * (cos(fourPi * ratio)));
                break;
            case FFTWindow::Blackman_Nuttall: // blackman nuttall
                weighingFactor = 0.3635819 - (0.4891775 * (cos(twoPi * ratio))) +
                    (0.1365995 * (cos(fourPi * ratio))) -
                    (0.0106411 * (cos(sixPi * ratio)));
                break;
            case FFTWindow::Blackman_Harris: // blackman harris
                weighingFactor = 0.35875 - (0.48829 * (cos(twoPi * ratio))) +
                    (0.14128 * (cos(fourPi * ratio))) -
                    (0.01168 * (cos(sixPi * ratio)));
                break;
            case FFTWindow::Flat_top: // flat top
                weighingFactor = 0.2810639 - (0.5208972 * cos(twoPi * ratio)) +
                    (0.1980399 * cos(fourPi * ratio));
                break;
            case FFTWindow::Welch: // welch
                weighingFactor = 1.0 - sq((indexMinusOne - samplesMinusOne / 2.0) /
                    (samplesMinusOne / 2.0));
                break;
            default:
                // This is Rectangle windowing which doesn't do anything
                // and Precompiled which shouldn't be selected
                break;
            }
            if (withCompensation) {
                weighingFactor *= compensationFactor;
            }
            if (windowingFactors) {
                windowingFactors[i] = weighingFactor;
            }
            if (dir == FFTDirection::Forward) {
                vData[i] *= weighingFactor;
                vData[samples - (i + 1)] *= weighingFactor;
            }
            else {
#ifdef FFT_SPEED_OVER_PRECISION
                T inverse = 1.0 / weighingFactor;
                vData[i] *= inverse;
                vData[samples - (i + 1)] *= inverse;
#else
                vData[i] /= weighingFactor;
                vData[samples - (i + 1)] /= weighingFactor;
#endif
            }
        }
    }
}

// Private functions

uint_fast8_t FFT::exponent(uint_fast16_t value) const {
    // Calculates the base 2 logarithm of a value
    uint_fast8_t result = 0;
    while (value >>= 1)
        result++;
    return result;
}

void FFT::findMaxY(float* vData, uint_fast16_t length, float* maxY,
    uint_fast16_t* index) const {
    *maxY = 0;
    *index = 0;
    // If sampling_frequency = 2 * max_frequency in signal,
    // value would be stored at position samples/2
    for (uint_fast16_t i = 1; i < length; i++) {
        if ((vData[i - 1] < vData[i]) && (vData[i] > vData[i + 1])) {
            if (vData[i] > vData[*index]) {
                *index = i;
            }
        }
    }
    *maxY = vData[*index];
}

void FFT::parabola(float x1, float y1, float x2, float y2, float x3, float y3, float* a, float* b,
    float* c) const {
    // const T reversed_denom = 1 / ((x1 - x2) * (x1 - x3) * (x2 - x3));
    // This is a special case in which the three X coordinates are three positive,
    // consecutive integers. Therefore the reverse denominator will always be -0.5
    const float reversed_denom = -0.5;

    *a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) * reversed_denom;
    *b = (x3 * x3 * (y1 - y2) + x2 * x2 * (y3 - y1) + x1 * x1 * (y2 - y3)) *
        reversed_denom;
    *c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 +
        x1 * x2 * (x1 - x2) * y3) *
        reversed_denom;
}

void FFT::swap(float* a, float* b) const {
    float temp = *a;
    *a = *b;
    *b = temp;
}

#ifdef FFT_SQRT_APPROXIMATION
// Fast inverse square root aka "Quake 3 fast inverse square root", multiplied
// by x. Uses one iteration of Halley's method for precision. See:
// https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Iterative_methods_for_reciprocal_square_roots
// And: https://github.com/HorstBaerbel/approx
template <typename T> float FFT::sqrt_internal(float x) const {
    union // get bits for floating point value
    {
        float x;
        int32_t i;
    } u;
    u.x = x;
    u.i = 0x5f375a86 - (u.i >> 1); // gives initial guess y0.
    float xu = x * u.x;
    float xu2 = xu * u.x;
    // Halley's method, repeating increases accuracy
    u.x = (0.125 * 3.0) * xu * (5.0 - xu2 * ((10.0 / 3.0) - xu2));
    return u.x;
}

template <typename T> double FFT::sqrt_internal(double x) const {
    // According to HosrtBaerbel, on the ESP32 the approximation is not faster, so
    // we use the standard function
#ifdef ESP32
    return sqrt(x);
#else
    union // get bits for floating point value
    {
        double x;
        int64_t i;
    } u;
    u.x = x;
    u.i = 0x5fe6ec85e7de30da - (u.i >> 1); // gives initial guess y0.
    double xu = x * u.x;
    double xu2 = xu * u.x;
    // Halley's method, repeating increases accuracy
    u.x = (0.125 * 3.0) * xu * (5.0 - xu2 * ((10.0 / 3.0) - xu2));
    return u.x;
#endif
}
#endif

const float FFT::_WindowCompensationFactors[10] = {
    1.0000000000 * 2.0, // rectangle (Box car)
    1.8549343278 * 2.0, // hamming
    1.8554726898 * 2.0, // hann
    2.0039186079 * 2.0, // triangle (Bartlett)
    2.8163172034 * 2.0, // nuttall
    2.3673474360 * 2.0, // blackman
    2.7557840395 * 2.0, // blackman nuttall
    2.7929062517 * 2.0, // blackman harris
    3.5659039231 * 2.0, // flat top
    1.5029392863 * 2.0  // welch
};
