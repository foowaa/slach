/*
此处主要是定义了一些安全的内存管理机制，如需内存调试，请使用slach_malloc_log, slach_free_log，这将会将所有使用过的记录放在文件中

*/

#ifndef __SLACH_H__
#define __SLACH_H__

#ifdef __cplusplus
    extern "C" {
#endif


#include <stdio.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <setjmp.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>


/*
some useful memory control functions: slach_malloc(type, size), slach_free(ptr)
 */
#define slach_malloc_(type, size) _slach_malloc(size, sizeof(type))
#define slach_malloc(type, size) (type*)slach_malloc_(type, size)
#define slach_free(ptr) _slach_free(ptr)

void* _slach_malloc(size_t size); 
void _slach_free(void* ptr);


//http://blog.csdn.net/hackbuteer1/article/details/6789164
//http://www.cnblogs.com/QG-whz/p/5140930.html
void* _slach_malloc(size_t n, size_t size){
    void* ptr = NULL;
    ptr = calloc(n, size);
    if(ptr == NULL){
        fprintf(stderr, "Fail to malloc!\n");
        exit(1);
    }
    return ptr;
}

void _slach_free(void* ptr){
    if (ptr == NULL){
        fprintf(stderr, "Fail to free!\n");
        exit(1);
    }
    free(ptr);
    ptr = NULL;
}


/*** Random number generation ***/

void slach_rand_seed(unsigned int seed);
int slach_rand_int_range(int min, int max);

void slach_rand_seed(unsigned int seed){
    if (seed != NULL)
        srand(seed);
    else
        srand(time(0));
}

int slach_rand_int_range(int max, int min){
    return rand()%min+max;
}

/*
Matrix and Vector
 */
#define T float
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



typedef struct _Matrix_
{
    unsigned int mHeight;
    unsigned int mWidth;
    T** mData;
}Matrix;

Matrix* createMatrix(IN unsigned int mHeight, IN unsigned int mWidth);
void destroyMatrix(INOUT Matrix* mPtr);
void copyMatrix(IN Matrix* src, OUT Matrix* dest);
void arrayToMatrix(IN T** src, OUT Matrix* dest, unsigned int height, unsigned int weight);
void matrixToArray(IN Matrix* src, OUT T** dest, unsigned int height, unsigned int weight);

typedef struct _Vector_
{
    unsigned int vLength;
    T *vData;
}Vector;

Vector* createVector(IN unsigned int vLength);
Vector* copyToVector(IN T* src, int len);
void destroyVector(INOUT Vector *vptr);
void copyVector(IN Vector* src, OUT Vector* dest); 
void arrayToVector(IN T* src, OUT Vector* dest);
void vectorToArray(IN Vector* src, OUT T* dest);




#ifdef __cplusplus
}
#endif