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

#include "operation.h"

/**< Matrix operations */
/** \brief matrix*matrix, private function
 *
 * \param 2-dim array, row, col
 * \param 2-dim array, row, col
 * \return Matrix*
 *
 */

Matrix* _mmMul(INOUT float* arr1, size_t row1, size_t col1, INOUT float* arr2, size_t row2, size_t col2){
    Matrix* m1;
    Matrix* m2;
    Matrix* result;
    size_t i,j,k;
    float sum;
    if (col1 != row2){
        perr("In mmMul(), col1 != row2!\n");
    }

    if (col1<=1 || row1<=1 || col2<=1 || row2<=1){
        perr("In mmMul(), col or row has problems!\n");
    }

    else{
        m1 = createMatrix(row1, col1);
        m2 = createMatrix(row2, col2);
        result = createMatrix(row1, col2);
        arrayToMatrix(arr1, m1, row1, col1);
        arrayToMatrix(arr2, m2, row2, col2);

        for (i = 0; i<row1; i++){
            for (j = 0; j<col2; j++){
                sum = 0;
                for (k = 0; k<row2; k++){
                    sum += m1->mData[i][k]*m2->mData[k][j];
                }
                result->mData[i][j] = sum;
            }
        }
        destroyMatrix(m1);destroyMatrix(m2);
        return result;
    }
}
/** \brief interface of matrix*matrix
 *
 * \param 2-dim array, row, col
 * \param 2-dim array, row, col
 * \param 2-dim array to save result, row, col
 * \return
 *
 */

void mmMul(INOUT float* arr1, size_t row1, size_t col1, INOUT float* arr2, size_t row2, size_t col2,
                                                        OUT float* dest, size_t height, size_t width){

    matrixToArray(_mmMul(arr1, row1, col1, arr2, row2, col2), dest, height, width);
}
/** \brief matrix*vector OR vector*matrix, private function
 *
 * \param 2-dim array, row, col OR 1-dim array, row, 1
 * \param 1-dim array, 1, col OR 2-dim array, row, col
 * \return Vector*
 *
 */

Vector* _mvMul(INOUT float* arr1, size_t row1, size_t col1, INOUT float* arr2, size_t row2, size_t col2){
    Matrix* m;
    Vector* v;
    Vector* result;
    size_t i,k;
    float sum;
    if (col1 != row2){
        perr("In mvMul(), col1 != row2!\n");
    }

    if ((col1 > 1 && row1>1) && (col2 == 1 && row2 > 1)){
        m = createMatrix(row1, col1);
        v = createVector(row2);
        result = createVector(row1);
        arrayToMatrix(arr1, m, row1, col1);
        arrayToVector(arr2, v, row2);

        for (i = 0; i<row1; i++){
            sum = 0;
            for (k = 0; k<row2; k++){
                sum += m->mData[i][k]*v->vData[k];
            }
            result->vData[i] = sum;
        }
        destroyMatrix(m); destroyVector(v);
        return result;
    }

    else if ((col1 > 1 && row1 == 1) && (col2 >1 || row2 > 1)){
        v = createVector(col1);
        m = createMatrix(row2, col2);
        result = createVector(col2);
        arrayToVector(arr1, v, col2);
        arrayToMatrix(arr2, m, row2, col2);

        for (i = 0; i<row2; i++){
            sum = 0;
            for (k = 0; k<col2; k++){
                sum += v->vData[i]*m->mData[i][k];
            }
            result->vData[i] = sum;
        }
        destroyVector(v); destroyMatrix(m);
        return result;
    }
    else{
        perr("In mvMul(), auguments are illegal!\n");
    }
}
/** \brief interface of matrix*vector OR vector*matrix
 *
 * \param 2-dim array, row, col OR 1-dim array, row, 1
 * \param 1-dim array, 1, col OR 2-dim array, row, col
 * \param 1-dim array to save result, len
 * \return
 *
 */

void mvMul(INOUT float* arr1, size_t row1, size_t col1, INOUT float* arr2, size_t row2, size_t col2,
                                                        OUT float* dest, size_t len){

    vectorToArray(_mvMul(arr1, row1, col1, arr2, row2, col2), dest, len);
}

/** \brief matrix+matrix, private function
 *
 * \param 2-dim array, row, col
 * \param 2-dim array, row, col
 * \return Matrix*
 *
 */

Matrix* _mmAdd(INOUT float* arr1, size_t row1, size_t col1, INOUT float* arr2, size_t row2, size_t col2){
    Matrix* m1;
    Matrix* m2;
    Matrix* result;
    size_t i,j;
    if (row1 != row2 || col1 != col2){
        perr("In mmAdd(), col or row is mismatched!\n");
    }
    m1 = createMatrix(row1, col1);
    m2 = createMatrix(row2, col2);
    result = createMatrix(row1, col1);
    arrayToMatrix(arr1, m1, row1, col1);
    arrayToMatrix(arr2, m2, row2, col2);

    for (i=0; i<row1; i++){
        for (j=0; j<col1; j++){
            result->mData[i][j] = m1->mData[i][j]+m2->mData[i][j];
        }
    }
    destroyMatrix(m1);destroyMatrix(m2);
    return result;
}
/** \brief interface of matrix+matrix
 *
 * \param 2-dim array, row, col
 * \param 2-dim array, row, col
 * \param 2-dim array to save result, row, col
 * \return
 *
 */
void mmAdd (INOUT float* arr1, size_t row1, size_t col1, INOUT float* arr2, size_t row2, size_t col2,
                                                        OUT float* dest, size_t height, size_t width){
    matrixToArray(_mmAdd(arr1, row1, col1, arr2, row2, col2), dest, height, width);
}
/** \brief vector+vector, private function
 *
 * \param 1-dim array, len
 * \param 1-dim array, len
 * \return Vector*
 *
 */

Vector* _vvAdd(INOUT float* arr1, size_t len1, INOUT float* arr2, size_t len2){
    Vector* v1;
    Vector* v2;
    Vector* result;
    size_t i;
    if (len1 != len2){
        perr("len1 != len2\n");
    }
    v1 = createVector(len1);
    v2 = createVector(len2);
    result = createVector(len1);
    arrayToVector(arr1, v1, len1);
    arrayToVector(arr2, v2, len2);
    for (i=0; i<len1; i++){
        result->vData[i] = v1->vData[i]+v2->vData[i];
    }
    destroyVector(v1); destroyVector(v2);
    return result;
}
/** \brief interface of vector+vector
 *
 * \param 1-dim array, len
 * \param 1-dim array, len
 * \param 1-dim array to save result, len
 * \return
 *
 */

void vvAdd (INOUT float* arr1, size_t len1, INOUT float* arr2, size_t len2,
                              OUT float* dest, size_t len){
    vectorToArray(_vvAdd(arr1, len1, arr2, len2), dest, len);
}
/** \brief interface of dot(vector, vector)
 *
 * \param 1-dim array, len
 * \param 1-dim array, len
 * \param 1-dim array to save result, len
 * \return float
 *
 */
float dot(INOUT float* arr1, size_t len1, INOUT float* arr2, size_t len2){
    Vector* v1;
    Vector* v2;
    float sum;
    size_t i;
    if (len1 != len2){
        perr("len1 != len2\n");
    }
    v1 = createVector(len1);
    v2 = createVector(len2);
    arrayToVector(arr1, v1, len1);
    arrayToVector(arr2, v2, len2);
    sum = 0;

    for (i=0; i<len1; i++){
        sum += (v1->vData[i])*(v2->vData[i]);
    }
    destroyVector(v1); destroyVector(v2);
    return sum;
}
/** \brief transpose(matrix), private function
 *
 * \param 2-dim array, row, col
 * \return Matrix*
 *
 */

Matrix* _mT (INOUT float* arr, size_t row, size_t col){
    Matrix* m = createMatrix(row, col);
    Matrix* result = createMatrix(col, row);
    size_t i,j;
    arrayToMatrix(arr, m, row, col);

    for (i=0; i<row; i++){
        for (j=0; j<col; j++){
            result->mData[j][i] = m->mData[i][j];
        }
    }
    destroyMatrix(m);
    return result;
}
/** \brief interface of transpose(matrix)
 *
 * \param 2-dim array, row, col
 * \param 2-dim array, row, col
 * \return
 *
 */

void mT (INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    matrixToArray(_mT(arr, row, col), dest, height, width);
}

/** \brief slice of matrix to vector, private function
 *
 * \param 2-dim array, row, col
 * \param 0/1 isRow, slice from row or col
 * \param start, end
 * \return Vector*
 *
 */

Vector* _slicev (INOUT float* arr, size_t row, size_t col, int isRow, int loc, int start, int end){
    Matrix* m;
    Vector* result;
    size_t i,j;
    if (end <= start || end < 0 || start < 0 || loc < 0 ){
        perr("In slicev, end or start has problems!\n");
    }
    m = createMatrix(row, col);
    arrayToMatrix(arr, m, row, col);
    result = createVector(end-start+1);
    if (isRow == 1){
        for (i=start, j=0; i<=end; i++, j++){
            result->vData[j] = m->mData[loc][i];
        }
        destroyMatrix(m);
        return result;
    }
    else{
        for (i=start, j=0; i<=end; i++, j++){
            result->vData[j] = m->mData[i][loc];
        }
        destroyMatrix(m);
        return result;
    }
}
/** \brief interface of slice of matrix to vector
 *
 * \param 2-dim array, row, col
 * \param 0/1 isRow, slice from row or col
 * \param start, end
 * \param 1-dim array, len
 * \return
 *
 */

void slicev(INOUT float* arr, size_t row, size_t col, int isRow, int loc, int start, int end,
                                                                 OUT float* dest, size_t len){
    vectorToArray(_slicev(arr, row, col, isRow, loc, start, end), dest, len);
}
/** \brief slice a matrix to a small matrix, private function
 *
 * \param 2-dim array, row, col
 * \param startRow, endRow, startCol, endCol
 * \return Matrix*
 *
 */

Matrix* _slicem (INOUT float* arr, size_t row, size_t col, int startRow, int endRow, int startCol, int endCol){
    Matrix* m;
    Matrix* result;
    size_t i,j,k,z;
    if (endRow <= startRow || endCol <= startCol || startRow < 0 || startCol < 0 || endRow < 0 || endCol < 0){
        perr("In slicem, end or start has problems!\n");
    }
    m = createMatrix(row, col);
    result = createMatrix(endRow-startRow+1, endCol-startCol+1);
    arrayToMatrix(arr, m, row, col);

    for (i=startRow, z=0; i<=endRow; i++, z++){
        for (j=startCol, k=0; j<=endCol; j++, k++){
            result->mData[z][k] = m->mData[i][j];
        }
    }
    destroyMatrix(m);
    return result;
}
/** \brief interface of slice a matrix to a small matrix
 *
 * \param 2-dim array, row, col
 * \param startRow, endRow, startCol, endCol
 * \param 2-dim array, row, col
 *
 */
void slicem(INOUT float* arr, size_t row, size_t col, int startRow, int endRow, int startCol, int endCol,
                                                                  OUT float* dest, size_t height, size_t width){
   matrixToArray(_slicem(arr, row, col, startRow, endRow, startCol, endCol), dest, height, width);
}

/** \brief vector l-p norm
 *
 * \param type: "inf" OR a number
 * \param 1-dim array, len
 * \return float
 *
 */

float vNorm (char* type, float* arr, size_t len){
    Vector* v;
    size_t i;
    float temp;
    int order;
    if (!strcmp(type, "0")) perr("type is not 0\n");
    v = createVector(len);
    arrayToVector(arr, v, len);
    //inf norm
    if (!strcmp(type, "inf")){
        temp = arr[0];
        for (i=1; i<len; i++){
            if (fabs((double)v->vData[i])>fabs(temp)){
                temp = (float)fabs((double)v->vData[i]);
            }
        }
        destroyVector(v);
        return temp;
    }
    else{
        temp = 0;
        order = atoi(type);
        for (i=0; i<len; i++){
            temp += pow(fabs((double)v->vData[i]), (double)order);
        }
        destroyVector(v);
        return (float)pow(temp, (double)1/order);
    }
}

/** \brief matrix norm. NOTE: at present, only Frobenius norm implements
 *
 * \param type: "F"
 * \param 2-dim array, row, col
 * \return float
 *
 */

float mNorm (char* type, float* arr, size_t row, size_t col){
    Matrix* m = createMatrix(row, col);
    size_t i,j;
    float sum;
    arrayToMatrix(arr, m, row, col);
    if (strcmp(type, "F")){
        //PASS
        perr("this type is undefined in mNorm\n");
    }
    else {
        sum = 0;
        for (i=0; i<row; i++){
            for (j=0; j<col; j++){
                sum += m->mData[i][j]*m->mData[i][j];
            }
        }
        destroyMatrix(m);
        return (float)pow(sum, 0.5);
    }
}

/** \brief element-wise math functions of matrix or vector
 *
 * \param 1-dim array, len; 2-dim array, row, col
 * \param 1-dim array to save result, len; 2-dim array to save result, row, col
 * \return
 *
 */

Vector* _absv (INOUT float* arr, size_t len){
    Vector* v = createVector(len);
    Vector* result = createVector(len);
    size_t i;
    arrayToVector(arr, v, len);

    for (i=0; i<len; i++){
        result->vData[i] = (float)fabs((double)v->vData[i]);
    }
    destroyVector(v);
    return result;
}

void absv(INOUT float* arr, size_t len, OUT float* dest, size_t lend){
    vectorToArray(_absv(arr, len), dest, lend);
}

Matrix* _absm (INOUT float* arr, size_t row, size_t col){
    Matrix* m = createMatrix(row, col);
    size_t i,j;
    Matrix* result = createMatrix(row, col);
    arrayToMatrix(arr, m, row, col);
    for (i=0; i<row; i++){
        for (j=0; j<col; j++){
            result->mData[i][j] = (float)fabs((double)m->mData[i][j]);
        }
    }
    destroyMatrix(m);
    return result;
}

void absm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    matrixToArray(_absm(arr, row, col), dest, height, width);
}

//element-wise sin
Vector* _sinv (INOUT float* arr, size_t len){
    Vector* v = createVector(len);
    Vector* result = createVector(len);
    size_t i;
    arrayToVector(arr, v, len);

    for (i=0; i<len; i++){
        result->vData[i] = (float)sin((double)v->vData[i]);
    }
    destroyVector(v);
    return result;
}

void sinv(INOUT float* arr, size_t len, OUT float* dest, size_t lend){
    vectorToArray(_sinv(arr, len), dest, lend);
}

Matrix* _sinm (INOUT float* arr, size_t row, size_t col){
    Matrix* m = createMatrix(row, col);
    Matrix* result = createMatrix(row, col);
    arrayToMatrix(arr, m, row, col);
    size_t i,j;
    for (i=0; i<row; i++){
        for (j=0; j<col; j++){
            result->mData[i][j] = (float)sin((double)m->mData[i][j]);
        }
    }
    destroyMatrix(m);
    return result;
}

void sinm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    matrixToArray(_sinm(arr, row, col), dest, height, width);
}

//element-wise cos
Vector* _cosv (INOUT float* arr, size_t len){
    Vector* v = createVector(len);
    Vector* result = createVector(len);
    size_t i;
    arrayToVector(arr, v, len);

    for (i=0; i<len; i++){
        result->vData[i] = (float)cos((double)v->vData[i]);
    }
    destroyVector(v);
    return result;
}

void cosv(INOUT float* arr, size_t len, OUT float* dest, size_t lend){
    vectorToArray(_cosv(arr, len), dest, lend);
}

Matrix* _cosm (INOUT float* arr, size_t row, size_t col){
    Matrix* m = createMatrix(row, col);
    Matrix* result = createMatrix(row, col);
    size_t i,j;
    arrayToMatrix(arr, m, row, col);
    for (i=0; i<row; i++){
        for (j=0; j<col; j++){
            result->mData[i][j] = (float)cos((double)m->mData[i][j]);
        }
    }
    destroyMatrix(m);
    return result;
}

void cosm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    matrixToArray(_cosm(arr, row, col), dest, height, width);
}

//element-wise tan
Vector* _tanv (INOUT float* arr, size_t len){
    Vector* v = createVector(len);
    Vector* result = createVector(len);
    size_t i;
    arrayToVector(arr, v, len);

    for (i=0; i<len; i++){
        result->vData[i] = (float)tan((double)v->vData[i]);
    }
    destroyVector(v);
    return result;
}

void tanv(INOUT float* arr, size_t len, OUT float* dest, size_t lend){
    vectorToArray(_tanv(arr, len), dest, lend);
}

Matrix* _tanm (INOUT float* arr, size_t row, size_t col){
    Matrix* m = createMatrix(row, col);
    Matrix* result = createMatrix(row, col);
    size_t i,j;
    arrayToMatrix(arr, m, row, col);
    for (i=0; i<row; i++){
        for (j=0; j<col; j++){
            result->mData[i][j] = (float)tan((double)m->mData[i][j]);
        }
    }
    destroyMatrix(m);
    return result;
}

void tanm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    matrixToArray(_tanm(arr, row, col), dest, height, width);
}

//element-wise asin
Vector* _asinv (INOUT float* arr, size_t len){
    Vector* v = createVector(len);
    Vector* result = createVector(len);
    size_t i;
    arrayToVector(arr, v, len);

    for (i=0; i<len; i++){
        result->vData[i] = (float)asin((double)v->vData[i]);
    }
    destroyVector(v);
    return result;
}

void asinv(INOUT float* arr, size_t len, OUT float* dest, size_t lend){
    vectorToArray(_asinv(arr, len), dest, lend);
}

Matrix* _asinm (INOUT float* arr, size_t row, size_t col){
    Matrix* m = createMatrix(row, col);
    Matrix* result = createMatrix(row, col);
    size_t i,j;
    arrayToMatrix(arr, m, row, col);
    for (i=0; i<row; i++){
        for (j=0; j<col; j++){
            result->mData[i][j] = (float)asin((double)m->mData[i][j]);
        }
    }
    destroyMatrix(m);
    return result;
}

void asinm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    matrixToArray(_asinm(arr, row, col), dest, height, width);
}

//element-wise acos
Vector* _acosv (INOUT float* arr, size_t len){
    Vector* v = createVector(len);
    Vector* result = createVector(len);
    size_t i;
    arrayToVector(arr, v, len);

    for (i=0; i<len; i++){
        result->vData[i] = (float)acos((double)v->vData[i]);
    }
    destroyVector(v);
    return result;
}

void acosv(INOUT float* arr, size_t len, OUT float* dest, size_t lend){
    vectorToArray(_acosv(arr, len), dest, lend);
}

Matrix* _acosm (INOUT float* arr, size_t row, size_t col){
    Matrix* m = createMatrix(row, col);
    Matrix* result = createMatrix(row, col);
    size_t i,j;
    arrayToMatrix(arr, m, row, col);
    for (i=0; i<row; i++){
        for (j=0; j<col; j++){
            result->mData[i][j] = (float)acos((double)m->mData[i][j]);
        }
    }
    destroyMatrix(m);
    return result;
}

void acosm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    matrixToArray(_acosm(arr, row, col), dest, height, width);
}

//element-wise atan
Vector* _atanv (INOUT float* arr, size_t len){
    Vector* v = createVector(len);
    Vector* result = createVector(len);
    size_t i;
    arrayToVector(arr, v, len);

    for (i=0; i<len; i++){
        result->vData[i] = (float)atan((double)v->vData[i]);
    }
    destroyVector(v);
    return result;
}

void atanv(INOUT float* arr, size_t len, OUT float* dest, size_t lend){
    vectorToArray(_atanv(arr, len), dest, lend);
}

Matrix* _atanm (INOUT float* arr, size_t row, size_t col){
    Matrix* m = createMatrix(row, col);
    Matrix* result = createMatrix(row, col);
    size_t i,j;
    arrayToMatrix(arr, m, row, col);
    for (i=0; i<row; i++){
        for (j=0; j<col; j++){
            result->mData[i][j] = (float)atan((double)m->mData[i][j]);
        }
    }
    destroyMatrix(m);
    return result;
}

void atanm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    matrixToArray(_atanm(arr, row, col), dest, height, width);
}

//element-wise exp
Vector* _expv (INOUT float* arr, size_t len){
    Vector* v = createVector(len);
    Vector* result = createVector(len);
    size_t i;
    arrayToVector(arr, v, len);

    for (i=0; i<len; i++){
        result->vData[i] = (float)exp((double)v->vData[i]);
    }
    destroyVector(v);
    return result;
}

void expv(INOUT float* arr, size_t len, OUT float* dest, size_t lend){
    vectorToArray(_expv(arr, len), dest, lend);
}

Matrix* _expm (INOUT float* arr, size_t row, size_t col){
    Matrix* m = createMatrix(row, col);
    Matrix* result = createMatrix(row, col);
    size_t i,j;
    arrayToMatrix(arr, m, row, col);
    for (i=0; i<row; i++){
        for (j=0; j<col; j++){
            result->mData[i][j] = (float)exp((double)m->mData[i][j]);
        }
    }
    destroyMatrix(m);
    return result;
}

void expm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    matrixToArray(_expm(arr, row, col), dest, height, width);
}

//element-wise log
Vector* _logv (INOUT float* arr, size_t len){
    Vector* v = createVector(len);
    Vector* result = createVector(len);
    size_t i;
    arrayToVector(arr, v, len);

    for (i=0; i<len; i++){
        result->vData[i] = (float)log((double)v->vData[i]);
    }
    destroyVector(v);
    return result;
}

void logv(INOUT float* arr, size_t len, OUT float* dest, size_t lend){
    vectorToArray(_logv(arr, len), dest, lend);
}

Matrix* _logm (INOUT float* arr, size_t row, size_t col){
    Matrix* m = createMatrix(row, col);
    Matrix* result = createMatrix(row, col);
    size_t i,j;
    arrayToMatrix(arr, m, row, col);
    for (i=0; i<row; i++){
        for (j=0; j<col; j++){
            result->mData[i][j] = (float)log((double)m->mData[i][j]);
        }
    }
    destroyMatrix(m);
    return result;
}

void logm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    matrixToArray(_logm(arr, row, col), dest, height, width);
}

//element-wise pow
Vector* _powv (INOUT float* arr, size_t len, double order){
    Vector* v = createVector(len);
    Vector* result = createVector(len);
    size_t i;
    arrayToVector(arr, v, len);

    for (i=0; i<len; i++){
        result->vData[i] = (float)pow((double)v->vData[i], order);
    }
    destroyVector(v);
    return result;
}

void powv(INOUT float* arr, size_t len, double order, OUT float* dest, size_t lend){
    vectorToArray(_powv(arr, len, order), dest, lend);
}

Matrix* _powm (INOUT float* arr, size_t row, size_t col, double order){
    Matrix* m = createMatrix(row, col);
    Matrix* result = createMatrix(row, col);
    size_t i,j;
    arrayToMatrix(arr, m, row, col);
    for (i=0; i<row; i++){
        for (j=0; j<col; j++){
            result->mData[i][j] = (float)pow((double)m->mData[i][j], order);
        }
    }
    destroyMatrix(m);
    return result;
}

void powm(INOUT float* arr, size_t row, size_t col, double order, OUT float* dest, size_t height, size_t width){
    matrixToArray(_powm(arr, row, col, order), dest, height, width);
}

//element-wise sqrt
Vector* _sqrtv (INOUT float* arr, size_t len){
    Vector* v = createVector(len);
    Vector* result = createVector(len);
    size_t i;
    arrayToVector(arr, v, len);

    for (i=0; i<len; i++){
        result->vData[i] = (float)sqrt((double)v->vData[i]);
    }
    destroyVector(v);
    return result;
}

void sqrtv(INOUT float* arr, size_t len, OUT float* dest, size_t lend){
    vectorToArray(_sqrtv(arr, len), dest, lend);
}

Matrix* _sqrtm (INOUT float* arr, size_t row, size_t col){
    Matrix* m = createMatrix(row, col);
    Matrix* result = createMatrix(row, col);
    size_t i,j;
    arrayToMatrix(arr, m, row, col);
    for (i=0; i<row; i++){
        for (j=0; j<col; j++){
            result->mData[i][j] = (float)sqrt((double)m->mData[i][j]);
        }
    }
    destroyMatrix(m);
    return result;
}

void sqrtm(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    matrixToArray(_sqrtm(arr, row, col), dest, height, width);
}
