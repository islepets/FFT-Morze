#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "defs.h"
//#include "types.h"
#include <math.h>
#include <stdint.h>

// This definition uses a low-precision square root approximation instead of the
// regular sqrt() call
// This might only work for specific use cases, but is significantly faster.

#ifndef FFT_SQRT_APPROXIMATION
#define sqrt_internal sqrt
#endif

enum class FFTDirection { Forward, Reverse };

enum class FFTWindow {
    Rectangle,        // rectangle (Box car)
    Hamming,          // hamming
    Hann,             // hann
    Triangle,         // triangle (Bartlett)
    Nuttall,          // nuttall
    Blackman,         // blackman
    Blackman_Nuttall, // blackman nuttall
    Blackman_Harris,  // blackman harris
    Flat_top,         // flat top
    Welch,            // welch
    Precompiled       // Placeholder for using custom or precompiled window values
};
#define FFT_LIB_REV 0x20
/* Custom constants */
/* These defines keep compatibility with pre 2.0 code */
#define FFT_FORWARD FFTDirection::Forward
#define FFT_REVERSE FFTDirection::Reverse

/* Windowing type */
#define FFT_WIN_TYP_RECTANGLE FFTWindow::Rectangle /* rectangle (Box car) */
#define FFT_WIN_TYP_HAMMING FFTWindow::Hamming     /* hamming */
#define FFT_WIN_TYP_HANN FFTWindow::Hann           /* hann */
#define FFT_WIN_TYP_TRIANGLE FFTWindow::Triangle   /* triangle (Bartlett) */
#define FFT_WIN_TYP_NUTTALL FFTWindow::Nuttall     /* nuttall */
#define FFT_WIN_TYP_BLACKMAN FFTWindow::Blackman   /* blackman */
#define FFT_WIN_TYP_BLACKMAN_NUTTALL                                           \
  FFTWindow::Blackman_Nuttall /* blackman nuttall */
#define FFT_WIN_TYP_BLACKMAN_HARRIS                                            \
  FFTWindow::Blackman_Harris                    /* blackman harris*/
#define FFT_WIN_TYP_FLT_TOP FFTWindow::Flat_top /* flat top */
#define FFT_WIN_TYP_WELCH FFTWindow::Welch      /* welch */
/* End of compatibility defines */

/* Mathematial constants */
#define twoPi 6.28318531
#define fourPi 12.56637061
#define sixPi 18.84955593

class FFT {
public:
    FFT();
    FFT(float* vReal, float* vImag, uint_fast16_t samples, float samplingFrequency,
        bool windowingFactors = false);

    ~FFT();

    void complexToMagnitude(void) const;
    void complexToMagnitude(float* vReal, float* vImag, uint_fast16_t samples) const;

    void compute(FFTDirection dir) const;
    void compute(float* vReal, float* vImag, uint_fast16_t samples,
        FFTDirection dir) const;
    void compute(float* vReal, float* vImag, uint_fast16_t samples, uint_fast8_t power,
        FFTDirection dir) const;

    void dcRemoval(void) const;
    void dcRemoval(float* vData, uint_fast16_t samples) const;

    float  majorPeak(void) const;
    void majorPeak(float* f, float* v) const;
    float majorPeak(float* vData, uint_fast16_t samples, float samplingFrequency) const;
    void majorPeak(float* vData, uint_fast16_t samples, float samplingFrequency,
        float* frequency, float* magnitude) const;

    float majorPeakParabola(void) const;
    void majorPeakParabola(float* frequency, float* magnitude) const;
    float majorPeakParabola(float* vData, uint_fast16_t samples,
        float samplingFrequency) const;
    void majorPeakParabola(float* vData, uint_fast16_t samples, float samplingFrequency,
        float* frequency, float* magnitude) const;

    uint8_t revision(void);

    void setArrays(float* vReal, float* vImag, uint_fast16_t samples = 0);

    void windowing(FFTWindow windowType, FFTDirection dir,
        bool withCompensation = false);
    void windowing(float* vData, uint_fast16_t samples, FFTWindow windowType,
        FFTDirection dir, float* windowingFactors = nullptr,
        bool withCompensation = false);

private:
    static const float _WindowCompensationFactors[10];
    bool _isPrecompiled = false;
    bool _precompiledWithCompensation = false;
    uint_fast8_t _power = 0;
    float* _precompiledWindowingFactors;
    uint_fast16_t _samples;
    float _samplingFrequency;
    float* _vImag;
    float* _vReal;
    FFTWindow _windowFunction;

    uint_fast8_t exponent(uint_fast16_t value) const;
    void findMaxY(float* vData, uint_fast16_t length, float* maxY,
        uint_fast16_t* index) const;
    void parabola(float  x1, float  y1, float  x2, float  y2, float  x3, float  y3, float* a, float* b, float* c) const;
    void swap(float* a, float* b) const;
};
