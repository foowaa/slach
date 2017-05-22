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

#ifndef OPERATION_H_
#define OPERATION_H_

#ifdef __cplusplus
    extern "C" {
#endif
#include "base.h"

/*
slice a matrix to a small matrix OR vector
*/
void slicev (INOUT float* arr, size_t row, size_t col, int isRow, int loc, int start, int end,
                                                                  OUT float* dest, size_t len);
void slicem (INOUT float* arr, size_t row, size_t col, int startRow, int endRow, int startCol, int endCol,
                                                                  OUT float* dest, size_t height, size_t width);
/*
element-wise math functions
*/
void absv(INOUT float* arr, size_t len, OUT float* dest, size_t lend);
void absm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width);
void sinv(INOUT float* arr, size_t len, OUT float* dest, size_t lend);
void sinm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width);
void cosv(INOUT float* arr, size_t len, OUT float* dest, size_t lend);
void cosm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width);
void asinv(INOUT float* arr, size_t len, OUT float* dest, size_t lend);
void asinm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width);
void acosv(INOUT float* arr, size_t len, OUT float* dest, size_t lend);
void acosm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width);
void atanv(INOUT float* arr, size_t len, OUT float* dest, size_t lend);
void atanm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width);
void expv(INOUT float* arr, size_t len, OUT float* dest, size_t lend);
void expm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width);
void logv(INOUT float* arr, size_t len, OUT float* dest, size_t lend);
void logm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width);
void powv(INOUT float* arr, size_t len, double order, OUT float* dest, size_t lend);
void powm(INOUT float* arr, size_t row, size_t col, double order, OUT float* dest, size_t height, size_t width);
void sqrtv(INOUT float* arr, size_t len, OUT float* dest, size_t lend);
void sqrtm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width);

/*
matrix multiplication, matrix addition, inner product of vectors, matrix transpose
*/
void mmMul(INOUT float* arr1, size_t row1, size_t col1, INOUT float* arr2, size_t row2, size_t col2,
                                                        OUT float* dest, size_t height, size_t width);
void mvMul(INOUT float* arr1, size_t row1, size_t col1, INOUT float* arr2, size_t row2, size_t col2,
                                                        OUT float* dest, size_t len);
void mmAdd(INOUT float* arr1, size_t row1, size_t col1, INOUT float* arr2, size_t row2, size_t col2,
                                                        OUT float* dest, size_t height, size_t width);
void vvAdd(INOUT float* arr1, size_t len1, INOUT float* arr2, size_t len2,
                              OUT float* dest, size_t len);
float dot(INOUT float* arr1, size_t len1, INOUT float* arr2, size_t len2);

/*
vector l-p norm and matrix Frobenius norm
*/
float vNorm (char* type, float* arr, size_t len);
float mNorm (char* type, float* arr, size_t row, size_t col);





#ifdef __cplusplus
}
#endif

#endif
