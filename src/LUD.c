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
#include "LUD.h"

/** \brief LUD result data structure
 *
 * \param LU
 * \param pivot
 * \return
 *
 */

typedef struct _LUDecRes_{
    Matrix* LU;
    Vector* piv;
}LUDecRes;

/** \brief LUD implementation, private function
 *
 * \param 2-dim array, row, col
 * \return LUDecRes
 *
 */

LUDecRes _LUdec(INOUT float* arr, size_t row, size_t col){
    LUDecRes result;
    Matrix* LU = createMatrix(row, col);
    size_t i,j,k;
    Vector* piv = createVector(row);
    Vector* LUrowi = createVector(row);
    Vector* LUcolj = createVector(col);
    size_t kmax,p;
    float s;
    if (row != col){
        perr("row != col in LUD!\n");
    }
    arrayToMatrix(arr, LU, row, col);
    for (i=0; i<row; i++){
        piv->vData[i] = i;
    }
    for (j=0; j<col; j++){
        for (i=0; i<row; i++){
            LUcolj->vData[i] = LU->mData[i][j];
        }
        for (i=0; i<row; i++){
            for (k=0; k<col; k++){
                LUrowi->vData[k] = LU->mData[i][k];
            }
            kmax = MAX(i,j);
            s = 0;
            for (k=0; k<kmax; k++){
                s += LUcolj->vData[k]*LUrowi->vData[k];
            }
            LUrowi->vData[j] -= s;
            LUcolj->vData[i] -= s;
        }
        p =j;
        for (i=j+1; i<row; i++){
            if (fabs((double)LUcolj->vData[i] > fabs((double)LUcolj->vData[p])))
                p = i;
        }
        if (p != j){
            for (k=0; k<col; k++){
                swap(&LU->mData[p][k], &LU->mData[j][k]);
            }
            swap(&piv->vData[p], &piv->vData[j]);
        }
        if ((j < row) && LU->mData[j][j] <= FLOAT_EPSILON && LU->mData[j][j] >= -FLOAT_EPSILON){
            for (i=j+1; i<row; i++){
                LU->mData[i][j] /= LU->mData[j][j];
            }
        }
    }
    destroyVector(LUrowi); destroyVector(LUcolj);
    result.LU = LU;
    result.piv = piv;
    return result;
}

/** \brief permute copy for LUsolvev, private function
 *
 * \param Vector* b
 * \param LUD result
 * \return Vector*
 *
 */

Vector* _permuteCopy1(Vector* v, LUDecRes res){
    size_t pivLen = res.piv->vLength;
    Vector* x = createVector(pivLen);
    size_t i;
    if (pivLen != v->vLength){
        perr("problem in _permuteCopy1\n");
    }
    for (i=0; i<pivLen; i++){
        x->vData[i] = v->vData[(int)res.piv->vData[i]];
    }
    return x;
}
/** \brief permute copy for LUsolvem, private function
 *
 * \param Matrix* B
 * \param LUD result
 * \return Matrix*
 *
 */

Matrix* _permuteCopy2(Matrix* m, LUDecRes res, size_t j0, size_t j1){
    int pivLen = res.piv->vLength;
    Matrix* X = createMatrix(pivLen, j1-j0+1);
    size_t i,j;
    for (i=0; i<pivLen; i++){
        for (j=j0; j<=j1; j++){
            X->mData[i][j-j0] = m->mData[(int)res.piv->vData[i]][j];
        }
    }
    return X;
}
/** \brief determine whether matrix is non-singular, private function
 *
 * \param LUD result
 * \return 0/1
 *
 */

int _isLUNonsingular(LUDecRes temp){
    size_t j;
    for (j=0; j<temp.LU->mHeight; j++){
        if ((float)fabs((double)temp.LU->mData[j][j]) <= FLOAT_EPSILON)
            return 0;
    }
    return 1;
}

/** \brief interface to get L
 *
 * \param 2-dim array, row, col
 * \param 2-dim array, row, col
 * \return
 *
 */

void getL(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    size_t p = MIN(row, col);
    Matrix* temp = createMatrix(row, p);
    LUDecRes LU = _LUdec(arr, row, col);
    size_t i,j;
    for (i=0; i<row; i++){
        for (j=0; j<i && j<col; j++){
            temp->mData[i][j] = LU.LU->mData[i][j];
        }
    }
    for (i=0; i<p; i++)
        temp->mData[i][i] = 1;

    destroyMatrix(LU.LU); destroyVector(LU.piv);
    matrixToArray(temp, dest, height, width);
}

/** \brief interface to get U
 *
 * \param 2-dim array, row, col
 * \param 2-dim array, row, col
 * \return
 *
 */
void getU(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    size_t p = MIN(row, col);
    Matrix* temp = createMatrix(p, col);
    LUDecRes LU = _LUdec(arr, row, col);
    size_t i,j;

    for (i=0; i<row; i++){
        for (j=0; j<col; j++){
            temp->mData[i][j] = LU.LU->mData[i][j];
        }
    }
    destroyMatrix(LU.LU); destroyVector(LU.piv);
    matrixToArray(temp, dest, height, width);
}

/** \brief interface to solve equations. AX = b
 *
 * \param 2-dim array, row, col
 * \param 1-dim array, len
 * \param 2-dim array to save result, row, col
 * \return
 *
 */
void LUsolvev(INOUT float* arr1, size_t row, size_t col, INOUT float* arr2, size_t len1,
              OUT float* dest, size_t len2){
    Vector* b = createVector(len1);
    arrayToVector(arr2, b, len1);
    LUDecRes temp = _LUdec(arr1, row, col);
    Vector* x = _permuteCopy1(b, temp);
    int i,k;

    if (len1 != row){
        perr("In LUsolvev, len1 != row\n");
    }
    if (!_isLUNonsingular(temp)){
        perr("In LUsolvev, arr1 is singular.\n");
    }
    for (k=0; k<col; k++){
        for (i=k+1; i<row; i++){
            x->vData[i] -= x->vData[k]*temp.LU->mData[i][k];
        }
    }
    for (k=col-1; k>=0; k--){
        x->vData[k] /= temp.LU->mData[k][k];
        for (i=0; i<k; i++){
            x->vData[i] -= x->vData[k]*temp.LU->mData[i][k];
        }
    }
    destroyMatrix(temp.LU); destroyVector(temp.piv); destroyVector(b);
    vectorToArray(x, dest, len2);
}
/** \brief interface to solve equations. AX = B.
 *
 * \param 2-dim array, row, col
 * \param 2-dim array, row, col
 * \param 2-dim array to save result, row, col
 * \return
 *
 */
void LUsolvem(INOUT float* arr1, size_t row1, size_t col1, INOUT float* arr2, size_t row2, size_t col2,
              OUT float* dest, size_t height, size_t width){
    // dimensions: A is mxn, X is nxk, B is mxk
    Matrix* B = createMatrix(row2, col2);
    arrayToMatrix(arr2, B, row2, col2);
    LUDecRes temp = _LUdec(arr1, row1, col1);
    size_t nx;
    Matrix* X;
    int i,j,k;
    if (col2 != row1) perr("In LUsolvem, len1 != row\n");
    if (!_isLUNonsingular(temp)){
        perr("In LUsolvem, arr1 is singular.\n");
    }
    nx = B->mWidth;
    X = _permuteCopy2(B, temp, 0, nx-1);

    for (k=0; k<row1; k++){
        for (i=k+1; i<col1; i++)
            for (j=0; j<nx; j++){
                X->mData[i][j] -= X->mData[k][j]*temp.LU->mData[i][k];
        }
    }

    for (k=col1-1; k>=0; k--){
        for (j=0; j<nx; j++){
            X->mData[k][j] /= temp.LU->mData[k][k];
        }
        for (i=0; i<k; i++)
            for (j=0; j<nx; j++)
                X->mData[i][j] -= X->mData[k][j]*temp.LU->mData[i][k];
    }

    destroyMatrix(temp.LU); destroyVector(temp.piv);destroyMatrix(B);
    matrixToArray(X, dest, height, width);
}


/** \brief inverse of matrix
 *
 * \param 2-dim array, row, col
 * \param 2-dim array to save result, row, col
 * \return
 *
 */

void inv(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    Matrix* B = _eyem(col);
    float b[col][col];
    if (row != col)  perr("inv needs squared matrix!\n");
    matrixToArray(B, b, col, col);
    LUsolvem(arr, row, col, b, col, col,  dest, height, width);
}
