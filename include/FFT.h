/*
=======================================================================
Simple Linear Algebra Header (SLACH)
The library provides some useful linear algebra algorithms implementations
for ANSI C:
Matrix and Vector
Element-wise math functions
Matrix multiplication, add, transpose, inverse, vector dot, norm, slice
Random functions: uniform distr., Gaussian distri., Exp distri., random numbers
                   generation seed settings, integer interval random numbers generation
Matrix decomposition: LU decomposition, QR decomposition, SVD decomposition and eigenvalue
                      decomposition
                      solve linear equations use LUD or QRD
Fast Fourier Transform
Some utilities: floor, ceil, round, divide, perr, printv, printvArr, printm, printmArr, MAX, MIN,
                swap, safe malloc, safe free


Author: cltian
Email: tianchunlin123@gmail.com
Version: 0.1
========================================================================


Copyright cltian

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*/
#ifndef FFT_H_
#define FFT_H_

#ifdef __cplusplus
    extern "C" {
#endif
#include "base.h"

/*
interface are abs and phase after fft
*/
typedef struct _complex_{
    float re;
    float im;
}complex;

void fftAbs(float* src, size_t len1, float* dest, size_t len2, int N1, int N2);
void fftPhase(float* src, size_t len1, float* dest, size_t len2, int N1, int N2);
void fftshift(float* src, int len, int N);
complex* DFT_naive(complex* x, int N);
complex* FFT_CooleyTukey(complex* input, int N, int N1, int N2);
#ifdef __cplusplus
}
#endif

#endif
