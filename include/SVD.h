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
#ifndef SVD_H_
#define SVD_H_

#ifdef __cplusplus
    extern "C" {
#endif
#include "base.h"
/*
SVD decomposition. Besides, SVD is the generalization of eigenvalue decomposition
*/
void getS(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t len);
void getV(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width);
void getUs(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width);

#ifdef __cplusplus
}
#endif

#endif
