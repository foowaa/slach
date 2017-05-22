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
#include "base.h"

/*
some utilities functions
*/

/** \brief print error
 *
 * \param str: message
 * \return no-return
 *
 */

void perr(char* str){
    fprintf(stderr, str);
    exit(1);
}
/** \brief ceil for m/n
 *
 * \param m
 * \param  n
 * \return integer
 *
 */
int ceilInt(int m, int n){
    int q;
    if(n != 0){
        q = m / n;
        if(m%n != 0)
            q += 1;
        return q;
    }
    else{
        perr("The dividend shouldn't be zero.\n");
        return 0;
    }
}
/** \brief floor for m/n
 *
 * \param m
 * \param n
 * \return integer
 *
 */

int floorInt (int m, int n){
    if (n != 0){
        return (int)m/n;
    }

    else{
        perr("The dividend shouldn't be zero.\n");
        return 0;
    }
}
/** \brief round for m/n
 *
 * \param m
 * \param n
 * \return integer
 *
 */

int roundInt (int m, int n){
    float q;
    if (n != 0){
        q = m/n;
        if (q-0.5 >= FLOAT_EPSILON)
            return q++;
        else{
            return q--;
        }
    }
    else{
        perr("The dividend shouldn't be zero.\n");
        return 0;
    }
}
/** \brief exact division for m/n
 *
 * \param m
 * \param n
 * \return a struct reminderAndRes=>res, reminder
 *
 */

remainderAndRes divide(int m, int n){
    int remainder;  //remainder
    int res;         //quotient
    if ( n!= 0){
        remainder = m%n;
        res = floorInt(m, n);
        remainderAndRes set;
        set.res = res;
        set.remainder = remainder;
        return set;
    }
    else{
        perr("The dividend shouldn't be zero.\n");
    }
}

/** \brief print matrix as 2-dim array
 *
 * \param array
 * \param row, col
 * \return no-return
 *
 */
void printmArr(float* arr, size_t row, size_t col){
    int limit = 20;  //the maximum number displays per line
    remainderAndRes set = divide(col, limit);
    int i,j,k;
    puts("Matrix:\n");
    for (i=0; i<row; i++){
        for (j=0; j<set.res+1; j++){
            if (j<set.res){
                for (k=0; k<limit; k++){
                    printf("%f ",arr[i*col+(k+j*limit)]);
                    if (k == limit-1){
                        puts("...\n");
                    }
                }
            }
            else{
                for (k=0; k<set.remainder; k++){
                    printf("%f ",arr[i*col+(k+j*limit)]);
                    if (k == set.remainder-1){
                        puts("\n");
                    }
                }
            }
        }
    }
    puts("\n");
}
/** \brief print vector as 1-dim array
 *
 * \param array
 * \param len
 * \return no-return
 *
 */

void printvArr (float* arr, size_t len){
    int limit = 20;
    remainderAndRes set = divide(len, limit);
    int i,j;
    puts("Vector:\n");
    for (i=0; i<set.res+1; i++){
        if (i<set.res){
            for (j=0; j<limit; j++){
                printf("%f ",arr[j+i*limit]);
                if (j == limit-1){
                    puts("...\n");
                }
            }
        }
        else{
            for (j=0; j<set.remainder; j++){
                printf("%f ", arr[j+i*limit]);
                if (j == set.remainder-1){
                    puts("\n");
                }
            }
        }
    }
    puts("\n");
}
/** \brief print matrix
 *
 * \param Matrix* matrix
 * \return no-return
 *
 */

void printm(Matrix* m){
    int limit = 20;
    size_t row = m->mHeight; size_t col = m->mWidth;
    remainderAndRes set = divide(col, limit);
    int i,j,k;
    puts("Matrix:\n");
    for (i=0; i<row; i++){
        for (j=0; j<set.res+1; j++){
            if (j<set.res){
                for (k=0; k<limit; k++){
                    printf("%f ", m->mData[i][k+j*limit]);
                    if (k == limit-1){
                        puts("...\n");
                    }
                }
            }
            else{
                for (k=0; k<set.remainder; k++){
                    printf("%f ", m->mData[i][k+j*limit]);
                    if (k == set.remainder-1){
                        puts("\n");
                    }
                }
            }
        }
    }
    puts("\n");
}
/** \brief print vector
 *
 * \param Vector* v
 * \return no-return
 *
 */

void printv(Vector* v){
    int limit = 20;
    int len = v->vLength;
    remainderAndRes set = divide(len, limit);
    int i,j;
    puts("Vector:\n");
    for (i=0; i<set.res+1; i++){
        if (i<set.res){
            for (j=0; j<limit; j++){
                printf("%f ",v->vData[j+i*limit]);
                if (j == limit-1){
                    puts("...\n");
                }
            }
        }
        else{
            for (j=0; j<set.remainder; j++){
                printf("%f ",v->vData[j+i*limit]);
                if (j == set.remainder-1){
                    puts("\n");
                }
            }
        }
    }
    puts("\n");
}

/** \brief safe malloc, private function
 *
 * \param n: number of malloc
 * \param size: sizeof(type)
 * \return void* ptr
 *
 */

void* _slach_malloc_(size_t n, size_t size){
    void* ptr = NULL;
    ptr = calloc(n, size);
    if(ptr == NULL){
        perr("Fail to malloc!\n");
    }
    return ptr;
}

/** \brief safe free, private function
 *
 * \param void* ptr: waiting for free
 * \return no-return
 *
 */

void _slach_free(void* ptr){
    if (ptr == NULL){
        perr("Fail to free!\n");
    }
    free(ptr);
    ptr = NULL;
}
/** \brief Integer interval r.v. generation. It is recommended that when use r.v. initialize seed
 *
 * \param unsigned int: seed
 * \return no-return
 *
 */

void slach_rand_seed(unsigned int seed){  //set seed, recommendations: when using r.v., it's better to reset seed
    if (seed != 0)
        srand(seed);
    else
        srand(time(0));
}
/** \brief (a,b) r.v. integer
 *
 * \param min
 * \param max
 * \return integer
 *
 */

int slach_rand_int_range_1(int min, int max){
    return rand()%(max-min)+min;
}
/** \brief [a,b] r.v. integer
 *
 * \param min
 * \param max
 * \return integer
 *
 */

int slach_rand_int_range_2(int min, int max){
    return rand()%(max-min+1)+min;
}
/** \brief (a,b] r.v. integer
 *
 * \param min
 * \param max
 * \return integer
 *
 */
int slach_rand_int_range_3(int min, int max){
    return rand()%(max-min)+min+1;
}
/** \brief [a,b) r.v. integer
 *
 * \param min
 * \param max
 * \return integer
 *
 */
int slach_rand_int_range_4(int min, int max){
    return rand()%(max-min+1)+min-1;
}

/** \brief normalize rand(), private function
 *
 * \param empty
 * \return float
 *
 */

float _surand()
{
  return( (float) rand()/RAND_MAX );
}
/** \brief uniform distribution r.v.
 *
 * \param low
 * \param high
 * \return [low, high] float
 *
 */

float uRand(float low, float high)
{
  return(low+(high-low)*_surand());
}
/** \brief exponential distribution
 *
 * \param lambda
 * \return float
 *
 */

float expRand(float lambda)
{
  float u,x;
  u=_surand();
  x=(-1/lambda)*log(u);
  return(x);
}
/** \brief Box-Muller Transform to generate Gaussian distribution r.v.
 *
 * \param mu
 * \param sigma
 * \return float
 *
 */

float gaussRand(float mu, float sigma)
{
  float theta,rsq,x;
  theta=uRand(0,2*PI);
  rsq=expRand(0.5);
  x=sqrt(rsq)*cos(theta);
  return x*sigma+mu;
}
/** \brief swap 2 float var.
 *
 * \param pointer to float
 * \param pointer to float
 * \return no-return
 *
 */

void swap(float* x, float* y){
    float temp;
    temp = *x;
    *x = *y;
    *y = temp;
}


/**<Matrix  */
/** \brief create matrix
 *
 * \param height
 * \param width
 * \return Matrix*
 *
 */

Matrix* createMatrix(IN size_t mHeight, IN size_t mWidth){
	Matrix* mPtr;
	size_t i;
	if (mHeight == 0 || mWidth == 0){
		perr("height != width\n");
	}
	else{
		mPtr = slach_malloc (Matrix, 1);
		mPtr->mData = slach_malloc(float*, mHeight);
		for (i = 0; i<mHeight; i++){
			mPtr->mData[i] = slach_malloc(float, mWidth);
		}
		mPtr->mHeight = mHeight;
		mPtr->mWidth = mWidth;
		return mPtr;
	}
}

/** \brief free matrix
 *
 * \param Matrix* ptr
 * \return no-return
 *
 */

void destroyMatrix(INOUT Matrix* mPtr){
	size_t i;
	if (mPtr == NULL){
		perr("ptr is NULL is free!\n");
	}
	else{
		for (i = 0; i<mPtr->mHeight; i++){
			slach_free(mPtr->mData[i]);
		}
		slach_free(mPtr->mData);
		slach_free(mPtr);
	}
}
/** \brief deep copy of src and dest
 *
 * \param Matrix* src
 * \param Matrix* dest
 * \return no-return
 *
 */

void copyMatrix(IN Matrix* src, OUT Matrix* dest){
	size_t i,j;
	if (src == NULL || dest == NULL){
		perr("src or dest is NULL in copy!\n");
	}
	else if (src->mHeight != dest->mHeight || src->mWidth != dest->mWidth){
		perr("The size of src and dest is mismatched! \n");
	}
	else{
		for (i=0; i<dest->mHeight; i++){
			for (j=0; j<dest->mWidth; j++){
                dest->mData[i][j] = src->mData[i][j];
			}
		}
	}
}
/** \brief 2-dim array to matrix
 *
 * \param 2-dim array src, height, width
 * \param Matrix* dest
 * \return no-return
 *
 */

void arrayToMatrix(IN float* src, OUT Matrix* dest, size_t height, size_t width){
	size_t i,j;
	if (src == NULL || dest == NULL){
		perr("src or dest is NULL in copy!\n");
	}
	else if (height != dest->mHeight || width != dest->mWidth){
		perr("The size of src and dest is mismatched! \n");
	}
	else{
		for (i=0; i<height; i++){
			for (j=0; j<width; j++){
				dest->mData[i][j] = src[width*i+j];
			}
		}
	}
}

/** \brief matrix to 2-dim array, then free matrix
 *
 * \param Matrix* src
 * \param 2-dim array, height, width
 * \return
 *
 */

void matrixToArray(IN Matrix* src, OUT float* dest, size_t height, size_t width){
	size_t i,j;
	if (src == NULL || dest == NULL){
		perr("src or dest is NULL in copy!\n");
	}
	else if (height != src->mHeight || width != src->mWidth){
		perr("The size of src and dest is mismatched! \n");
	}
	else{
		for (i=0; i<height; i++){
			for (j=0; j<width; j++){
				dest[i*width+j] = src->mData[i][j];
			}
		}
		destroyMatrix(src);
	}
}

/** \brief matrix to 2-dim array without free
 *
 * \param Matrix* src
 * \param 2-dim array, height, width
 * \return
 *
 */
void matrixToArrayWithoutFree(IN Matrix* src, OUT float* dest, size_t height, size_t width){
	size_t i,j;
	if (src == NULL || dest == NULL){
		perr("src or dest is NULL in copy!\n");
	}
	else if (height != src->mHeight || width != src->mWidth){
		perr("The size of src and dest is mismatched! \n");
	}
	else{
		for (i=0; i<height; i++){
			for (j=0; j<width; j++){
				dest[i*width+j] = src->mData[i][j];
			}
		}
	}
}
/** \brief create and assign matrix with all num, private function
 *
 * \param row, col
 * \param num
 * \return Matrix*
 *
 */

Matrix* _assignm(size_t row, size_t col, float num){
    size_t i,j;
    Matrix* m = createMatrix(row, col);
    for (i=0; i<row; i++){
        for (j=0; j<col; j++){
            m->mData[i][j] = num;
        }
    }
    return m;
}

/** \brief create a eye matrix
 *
 * \param n
 * \return Matrix*
 *
 */

Matrix* _eyem(size_t n){
    size_t i;
    Matrix* eye = createMatrix(n, n);
    for (i=0; i<n; i++){
        eye->mData[i][i] = 1;
    }
    return eye;
}


/**< Vector */
/** \brief create vector
 *
 * \param len
 * \return Vector*
 *
 */

Vector* createVector(IN size_t vLength){
    Vector* vPtr;
    if (vLength == 0){
        perr("The size of src and dest is mismatched! \n");
    }
    else{
        vPtr = slach_malloc(Vector, 1);
        vPtr->vData = slach_malloc(float, vLength);
        vPtr->vLength = vLength;
        return vPtr;
    }
}

/** \brief free vector
 *
 * \param Vector* ptr
 * \return
 *
 */

void destroyVector(INOUT Vector* vPtr){
    if (vPtr == NULL){
        perr("ptr is NULL is free!\n");
    }
    else{
        slach_free(vPtr->vData);
        slach_free(vPtr);
    }
}
/** \brief deep copy src to dest
 *
 * \param Vector* src
 * \param Vector* dest
 * \return
 *
 */

void copyVector(IN Vector* src, OUT Vector* dest){
    size_t i;
    if (src == NULL || dest == NULL){
        perr("src or dest is NULL in copy!\n");
    }
    else if (src->vLength != dest->vLength){
        perr("The size of src and dest is mismatched! \n");
    }
    else{
        for (i=0; i<dest->vLength; i++)
            dest->vData[i] = src->vData[i];
    }
}
/** \brief 1-dim array to vector
 *
 * \param 1-dim array src, len
 * \param Vector* dest
 * \return
 *
 */

void arrayToVector(IN float *src, OUT Vector* dest,size_t len){
    size_t i;
    if (src == NULL || dest == NULL){
        perr("src or dest is NULL in copy!\n");
    }
    else if (len != dest->vLength){
        perr("The size of src and dest is mismatched! \n");
    }
    else{
        for (i = 0; i<len; i++){
            dest->vData[i] = src[i];
        }
    }
}

/** \brief vector to matrix, then free vector
 *
 * \param Vector* src
 * \param 1-dim array, len
 * \return
 *
 */

void vectorToArray(IN Vector* src, OUT float* dest, size_t len){
    size_t i;
    if (src == NULL || dest == NULL){
        perr("src or dest is NULL in copy!\n");
    }
    else if (src->vLength != len){
        perr("The size of src and dest is mismatched! \n");
    }
    else{
        for (i = 0; i<len; i++){
            dest[i] = src->vData[i];
        }
        destroyVector(src);
    }
}

/** \brief vector to matrix without free vector
 *
 * \param Vector* src
 * \param 1-dim array, len
 * \return
 *
 */

void vectorToArrayWithoutFree(IN Vector* src, OUT float* dest, size_t len){
    size_t i;
    if (src == NULL || dest == NULL){
        perr("src or dest is NULL in copy!\n");
    }
    else if (src->vLength != len){
        perr("The size of src and dest is mismatched! \n");
    }
    else{
        for (i = 0; i<len; i++){
            dest[i] = src->vData[i];
        }
    }
}


