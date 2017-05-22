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

#ifndef BASE_H_
#define BASE_H_

#ifdef __cplusplus
    extern "C" {
#endif

//All C89 libraries
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

/*
some useful memory control functions: slach_malloc(type, size), slach_free(ptr)
 */
#define slach_malloc_(type, size) _slach_malloc_(size, sizeof(type))
#define slach_malloc(type, size) (type*)slach_malloc_(type, size)
#define slach_free(ptr) _slach_free(ptr)
void* _slach_malloc(size_t size);
void _slach_free(void* ptr);

/*
random variables generation
 */
//set random seed
void slach_rand_seed(unsigned int seed);
//(a,b) integer
int slach_rand_int_range_1(int min, int max);
//[a,b] integer
int slach_rand_int_range_2(int min, int max);
//(a,b] integer
int slach_rand_int_range_3(int min, int max);
//[a,b) integer
int slach_rand_int_range_4(int min, int max);
//uniform distribution
float uRand(float low, float high);
//Gaussian distribution
float gaussRand(float mu, float sigma);
//exponential distribution
float expRand(float lambda);
//swap
void swap(float* x, float* y);


//Macros
#define IN
#define OUT
#define INOUT
#define MIN(x,y) (((x)>(y))?(y):(x))
#define MAX(x,y) (((x)<(y))?(y):(x))
#define SIZE_ARRAY_1(arr) sizeof(arr)/sizeof(arr[0])
#define SIZE_ARRAY_2_ROW(arr) sizeof(arr)/sizeof(arr[0])
#define SIZE_ARRAY_2_COL(arr) sizeof(arr[0])/sizeof(arr[0][0])
#define FLOAT_MAX FLT_MAX
#define FLOAT_MIN FLT_MIN
#define FLOAT_EPSILON FLT_EPSILON
#define PI 3.14159265359
#define E  2.71828182846
#define FLOAT_EQUY(x, y) fabs((x)-(y))<=FLOAT_EPSILON


/*
Base types
*/
//The Matrix base
typedef struct _Matrix_
{
    size_t mHeight;
    size_t mWidth;
    float** mData;
}Matrix;

Matrix* createMatrix(IN size_t mHeight, IN size_t mWidth);
void destroyMatrix(INOUT Matrix* mPtr);
void copyMatrix(IN Matrix* src, OUT Matrix* dest);
void arrayToMatrix(IN float* src, OUT Matrix* dest, size_t height, size_t weight);
void matrixToArray(IN Matrix* src, OUT float* dest, size_t height, size_t weight);
void matrixToArrayWithoutFree(IN Matrix* src, OUT float* dest, size_t height, size_t width);

//The Vector base
typedef struct _Vector_
{
    size_t vLength;
    float *vData;
}Vector;

Vector* createVector(IN size_t vLength);
Vector* copyToVector(IN float* src, int len);
void destroyVector(INOUT Vector *vptr);
void copyVector(IN Vector* src, OUT Vector* dest);
void arrayToVector(IN float* src, OUT Vector* dest, size_t len);
void vectorToArray(IN Vector* src, OUT float* dest, size_t len);
void vectorToArrayWithoutFree(IN Vector* src, OUT float* dest, size_t len);

/*
some utilities functions
*/
typedef struct _REAMINDER_RES_{
    int res;
    int remainder;
}remainderAndRes;
void perr(char* str); //print error and abort program
void printmArr(float* arr, size_t row, size_t col); //print matrix as 2-dim array
void printvArr(float* arr, size_t len); //print vector as 1-dim array
void printm(Matrix* m); //print matrix
void printv(Vector* v); //print vector
int ceilInt( int m, int n );
int floorInt( int m, int n);
int roundInt (int m, int n);
remainderAndRes divide(int m, int n); //exact division

#ifdef __cplusplus
}
#endif

#endif
