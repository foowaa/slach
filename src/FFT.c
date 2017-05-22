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
#include "../include/FFT.h"

/*
FFT
https://github.com/jtfell/c-fft

/**< private functions for complex operations */
complex _conv_from_polar(double r, double radians) {
    complex result;
    result.re = r * cos(radians);
    result.im = r * sin(radians);
    return result;
}
complex _cadd(complex left, complex right) {
    complex result;
    result.re = left.re + right.re;
    result.im = left.im + right.im;
    return result;
}
complex _cmultiply(complex left, complex right) {
    complex result;
    result.re = left.re*right.re - left.im*right.im;
    result.im = left.re*right.im + left.im*right.re;
    return result;
}
/** \brief abs of complex
 *
 * \param complex
 * \return float
 *
 */

float cAbs(complex a){
    return  sqrt(a.re * a.re + a.im * a.im);
}
/** \brief phase of complex
 *
 * \param complex
 * \return rad
 *
 */

float cPhase(complex a){
    return atan(a.im/a.re);
}
/** \brief naive Discrete Fourier Transform
 *
 * \param complex*
 * \param point
 * \return complex*
 *
 */

complex* DFT_naive(complex* x, int N) {
    complex* X = slach_malloc(complex, N);
    int k, n;
    for(k = 0; k < N; k++) {
        X[k].re = 0.0;
        X[k].im = 0.0;
        for(n = 0; n < N; n++) {
            X[k] = _cadd(X[k], _cmultiply(x[n], _conv_from_polar(1, -2*PI*n*k/N)));
        }
    }
    return X;
}
/** \brief Implements the Cooley-Tukey FFT algorithm.
 *   Cooley-Tukey FFT algorithm re-express DFT of an arbitrary composite size N = N1*N2
 *   in terms of N1 smaller DFTs of sizes N2, recursively.
 * \param complex*, points-N
 * \param N=N1*N2, ref: https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
 * \return complex*
 *
 */

complex* FFT_CooleyTukey(complex* input, int N, int N1, int N2) {
    int k1, k2;
    /* Allocate columnwise matrix */
    complex** columns = slach_malloc(complex*, N1);
    complex** rows;
    complex* output;
    for(k1 = 0; k1 < N1; k1++) {
        columns[k1] = slach_malloc(complex, N2);
    }
    /* Allocate rowwise matrix */
    rows = slach_malloc(complex*, N2);
    for(k2 = 0; k2 < N2; k2++) {
        rows[k2] = slach_malloc(complex, N1);
    }
    /* Reshape input into N1 columns */
    for (k1 = 0; k1 < N1; k1++) {
        for(k2 = 0; k2 < N2; k2++) {
            columns[k1][k2] = input[N1*k2 + k1];
        }
    }
    /* Compute N1 DFTs of length N2 using naive method */
    for (k1 = 0; k1 < N1; k1++) {
        columns[k1] = DFT_naive(columns[k1], N2);
    }
    /* Multiply by the twiddle factors  ( e^(-2*pi*j/N * k1*k2)) and transpose */
    for(k1 = 0; k1 < N1; k1++) {
        for (k2 = 0; k2 < N2; k2++) {
            rows[k2][k1] = _cmultiply(_conv_from_polar(1, -2.0*PI*k1*k2/N), columns[k1][k2]);
        }
    }
    /* Compute N2 DFTs of length N1 using naive method */
    for (k2 = 0; k2 < N2; k2++) {
        rows[k2] = DFT_naive(rows[k2], N1);
    }
    /* Flatten into single output */
    output = slach_malloc(complex, N);
    for(k1 = 0; k1 < N1; k1++) {
        for (k2 = 0; k2 < N2; k2++) {
            output[N2*k1 + k2] = rows[k2][k1];
        }
    }
    /* Free all alocated memory except output and input arrays */
    for(k1 = 0; k1 < N1; k1++) {
        slach_free(columns[k1]);
    }
    for(k2 = 0; k2 < N2; k2++) {
        slach_free(rows[k2]);
    }
    slach_free(columns);
    slach_free(rows);
    return output;
}

/** \brief interface of FFT to calculate abs
 *
 * \param 1-dim array, len
 * \param 1-dim array to save result, len
 * \return N = N1*N2
 *
 */

void fftAbs(float* src, size_t len1, float* dest, size_t len2, int N1, int N2){
    int i;
    complex* input = slach_malloc(complex, len1);
    complex* res;
    for (i=0; i<len1; i++){
        input[i].re = i;
        input[i].im = 0;
    }
    res = FFT_CooleyTukey(input, len1, N1, N2);
    for (i=0; i<len2; i++)
        dest[i] = cAbs(res[i]);
    slach_free(res);
}
/** \brief interface of FFT to calculate phase
 *
 * \param 1-dim array, len
 * \param 1-dim array to save result, len
 * \return N = N1*N2
 *
 */
void fftPhase(float* src, size_t len1, float* dest, size_t len2, int N1, int N2){
    int i;
    complex* input = slach_malloc(complex, len1);
    complex* res;
    for (i=0; i<len1; i++){
        input[i].re = i;
        input[i].im = 0;
    }
    res = FFT_CooleyTukey(input, len1, N1, N2);

    for (i=0; i<len2; i++)
        dest[i] = cPhase(res[i]);
    slach_free(res);
}
/** \brief center the DC
 *
 * \param 1-dim array, len
 * \param points
 * \return
 *
 */

void fftshift(float* sig, int len, int N){
    int i;
    int N1 = N>>1;
    float temp[N];
    if (len < N)
        perr("in fftshift, len < N\n");
    for (i=0; i<N1; i++){
        temp[i] = sig[i+N1];
        sig[i+N1] = sig[i];
        sig[i] = temp[i];
    }

}
