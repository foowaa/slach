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


#include "../include/QRD.h"

/** \brief QRD result data structure
 *
 * \param QR
 * \param R diag
 * \return
 *
 */
typedef struct _QRDecRes_{
    Matrix* QR;
    Vector* RDiag;
}QRDecRes;
/** \brief QRD implementation, private function
 *
 * \param 2-dim array, row, col
 * \return LUDecRes
 *
 */
QRDecRes _QRdec(INOUT float* arr, size_t row, size_t col){
    size_t m = row;
    size_t n = col;
    size_t p = MIN(m, n);
    QRDecRes result;
    Matrix* QR;
    Vector* RDiag;
    int i,j,k;
    float nrm;
    float s;
    if (row != col){
        perr("row != col in QRD!\n");
    }
    QR = createMatrix(row, col);
    arrayToMatrix(arr, QR, row, col);
    RDiag = createVector(p);

    for (k=0; k<p; k++){
        nrm = 0;
        for (i=k; i<m; i++)
            nrm = (float)sqrt(nrm*nrm+QR->mData[i][k]*QR->mData[i][k]);
        if (fabs(nrm)>=FLOAT_EPSILON){
            if (QR->mData[k][k] < 0)
                nrm = -nrm;
            for (i=k; i<m; i++){
                QR->mData[i][k] /= nrm;
            }
            QR->mData[k][k] += 1;

            for (j=k+1; j<n; j++){
                s = 0;
                for (i=k; i<m; i++)
                    s += QR->mData[i][k]*QR->mData[i][j];
                s = -s/QR->mData[k][k];
                for (i=k; i<m; i++){
                    QR->mData[i][j] += s*QR->mData[i][k];
                }
            }
        }
        RDiag->vData[k] = -nrm;
    }
    result.QR = QR; result.RDiag = RDiag;
    return result;
}
/** \brief determine whether matrix is full rank, private function
 *
 * \param QRD result
 * \return 0/1
 *
 */
int _isFullRank(QRDecRes temp){
    size_t j;
    for (j=0; j<temp.RDiag->vLength; j++){
        if (temp.RDiag->vData[j] == 0)
            return 0;
    }
    return 1;
}
/** \brief interface to get Q
 *
 * \param 2-dim array, row, col
 * \param 2-dim array, row, col
 * \return
 *
 */
void getQ(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    size_t m = row;
    size_t p = MIN(row, col);
    int i,j,k;
    float s;
    Matrix* Q = createMatrix(m, p);
    QRDecRes temp = _QRdec(arr, row, col);
    for (k=p-1; k>=0; k--){
        for (i=0; i<m; i++)
            Q->mData[i][k] = 0;
        Q->mData[k][k] = 1;
        for (j=k; j<p; j++){
            if (fabs(temp.QR->mData[k][k])>FLOAT_EPSILON){
                s = 0;
                for (i=k; i<m; i++){
                    s += temp.QR->mData[i][k]*Q->mData[i][j];
                }
                s = -s/temp.QR->mData[k][k];
                for (i=k; i<m; i++)
                    Q->mData[i][j] += s*temp.QR->mData[i][k];
            }
        }
    }
    matrixToArray(Q, dest, height, width);
    destroyMatrix(temp.QR);
    destroyVector(temp.RDiag);
}
/** \brief interface to get R
 *
 * \param 2-dim array, row, col
 * \param 2-dim array, row, col
 * \return
 *
 */
void getR(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    size_t n = col;
    size_t p = MIN(row, col);
    int i,j;
    Matrix* R = createMatrix(p, n);
    QRDecRes temp = _QRdec(arr, row, col);
    for (i=0; i<p; i++){
        for (j=0; j<col; j++){
            if (i<j){
                R->mData[i][j] = temp.QR->mData[i][j];
            }
            if (i==j){
                R->mData[i][j] = temp.RDiag->vData[i];
            }
        }
    }
    matrixToArray(R, dest, height, width);
    destroyMatrix(temp.QR);destroyVector(temp.RDiag);
}


/** \brief interface to solve equations. AX = b
 *
 * \param 2-dim array, row, col
 * \param 1-dim array, len
 * \param 2-dim array to save result, row, col
 * \return
 *
 */
void QRsolvev(INOUT float* arr1, size_t row, size_t col, INOUT float* arr2, size_t len1,
              OUT float* dest, size_t len2){
    QRDecRes temp = _QRdec(arr1, row, col);
    size_t m = temp.QR->mHeight;
    size_t n = temp.QR->mWidth;
    int i,k,s;
    Vector* x_;
    if (!_isFullRank(temp))
        perr("in QRD, arr1 is full rank!\n");
    Vector* x = createVector(len1);
    arrayToVector(arr2, x, len1);

    for (k=0; k<n; k++){
        s = 0;
        for (i=k; i<m; i++)
            s += temp.QR->mData[i][k]*x->vData[i];
        s = -s/temp.QR->mData[k][k];
        for (i=k; i<m; i++)
            x->vData[i] += s*temp.QR->mData[i][k];
        }

    for (k=n-1; k>=0; k--){
        x->vData[k] /= temp.RDiag->vData[k];
        for (i=0; i<k; i++)
            x->vData[i] -= x->vData[k]*temp.QR->mData[i][k];
    }

    x_ = createVector(n);
    for (i=0; i<n; i++)
        x_->vData[i] = x->vData[i];
    vectorToArray(x_, dest, len2);
    destroyMatrix(temp.QR);destroyVector(temp.RDiag);destroyVector(x_);

}

/** \brief interface to solve equations. AX = B
 *
 * \param 2-dim array, row, col
 * \param 2-dim array, row, col
 * \param 2-dim array to save result, row, col
 * \return
 *
 */
void QRsolvem(INOUT float* arr1, size_t row1, size_t col1, INOUT float* arr2, size_t row2, size_t col2,
              OUT float* dest, size_t height, size_t width){
    QRDecRes temp = _QRdec(arr1, row1, col1);
    size_t m = temp.QR->mHeight;
    size_t n = temp.QR->mWidth;
    size_t nx = col2;
    int i,j,k;
    Matrix* X = createMatrix(row2, col2);
    arrayToMatrix(arr2, X, row2, col2);
    Matrix* X_;
    float s;
    if (!_isFullRank(temp))
        perr("in QRD, arr1 is full rank!\n");
    for (k=0; k<n; k++){
        for (j=0; j<nx; j++){
            s = 0;
            for (i=k; i<m; i++){
                s += temp.QR->mData[i][k]*X->mData[i][j];
            }
            s = -s/temp.QR->mData[k][k];
            for (i=k; i<m; i++)
                X->mData[i][j] += s*temp.QR->mData[i][k];
        }
    }

    for (k=n-1; k>=0; k--){
        for (j=0; j<nx; j++)
            X->mData[k][j] /= temp.RDiag->vData[k];
        for (i=0; i<k; i++)
            for (j=0; j<nx; j++)
                X->mData[i][j] -= X->mData[k][j]*temp.QR->mData[i][k];
    }

    X_ = createMatrix(n, nx);
    for (i=0; i<n; i++)
        for (j=0; j<nx; j++)
            X_->mData[i][j] = X->mData[i][j];

    matrixToArray(X_, arr2, height, width);
    destroyMatrix(temp.QR);destroyVector(temp.RDiag);destroyMatrix(X);

}

