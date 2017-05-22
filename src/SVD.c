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
#include "SVD.h"

/** \brief SVD result data structure
 *
 * \param S
 * \param V
 * \param U
 * \return
 *
 */
typedef struct _SVD_{
    Matrix* U;
    Matrix* V;
    Vector* S;
}SVD;
/** \brief private function to generate V
 *
 * \param Matrix* A, Vector* v_
 * \return
 *
 */

void _svd_1d(Matrix* A, Vector* v_){
	int n = A->mHeight;
	int m = A->mWidth;
	int i,j,k;
	Matrix* B ;
	float sum, norm;
	Vector* currentV = createVector(m);
    int iter = 0;
	float epsilon = pow(10,-10);
	Vector* lastV;
	float norm2;
	slach_rand_seed(0);
	for (i=0; i<m; i++){
        currentV->vData[i] = gaussRand(0,1);
	}
	//currentV->vData[0] = -0.86; currentV->vData[1] = -0.04; currentV->vData[2] = 0.5;
	sum = 0;
	for (i=0; i<m; i++){
        sum += currentV->vData[i]*currentV->vData[i];
	}
	sum = sqrt(sum);
	for (i=0; i<m; i++){
        currentV->vData[i] = currentV->vData[i]/sum;
	}
	if (n>m){
		B = createMatrix(m,m);
		for (i=0; i<m; i++){
			for (j=0; j<m; j++){
                sum = 0;
                for (k=0; k<n; k++){
                    sum += A->mData[k][i]*A->mData[k][j];
                }
				B->mData[i][j] = sum;
			}
		}
	}
	else{
		B = createMatrix(n,n);
		for (i=0; i<n; i++){
			for (j=0; j<n; j++){
                    sum = 0;
                    for (k=0; k<m; k++){
                        sum += A->mData[i][k]*A->mData[j][k];
                    }
				B->mData[i][j] = sum;
			}
		}
	}

    lastV = createVector(currentV->vLength);
	while (1){
		norm2 = 0;
		copyVector(currentV, lastV);
		for (i=0; i<B->mHeight; i++){
			sum = 0;
			for (j=0; j<B->mWidth; j++){
				sum += lastV->vData[j] * B->mData[i][j];
			}
			norm2 += sum*sum;
			currentV->vData[i] = sum;
		}
		norm = sqrt(norm2);
		for (i=0; i<currentV->vLength; i++){
			currentV->vData[i] /= norm;
		}
		sum = 0;
		for (i=0; i<currentV->vLength; i++){
			sum += currentV->vData[i]*lastV->vData[i];
		}
		iter++;

		if (fabs(sum) > 1-epsilon){
            copyVector(currentV, v_);
			destroyVector(currentV); destroyVector(lastV);
			break;
		}

	}
}
/** \brief SVD implementation, private function
 *
 * \param 2-dim array, row, col
 * \return SVD
 *
 */

SVD _SVDdec(float* src, size_t row, size_t col){
	size_t n = row;
	size_t m = col;
	size_t k = MIN(n, m);
	SVD svd;
	svd.S = createVector(k);
	svd.U = createMatrix(row, k);
	svd.V = createMatrix(k, col);
	typedef struct _svdSoFar_{
		Vector* u;
		Vector* v;
		float singularValue;
	}svdSoFar_[k];
	Matrix* matrixFor1D = createMatrix(row, col);
	svdSoFar_ svdSoFar;
	int i,j;
	int p,q;
	Vector* v;
	Vector* u;
	Vector* v_;
	Vector* u_;
	float singularValue;
    float u_unnormalized_val;
    float sigma2;
    float sigma;
	for (i=0; i<k; i++){
		arrayToMatrix (src, matrixFor1D, row, col);
		for (j=0; j<i; j++){
			v = svdSoFar[j].v;
			u = svdSoFar[j].u;
			singularValue = svdSoFar[j].singularValue;
			for (p=0; p<row; p++){
				for (q=0; q<col; q++){
					matrixFor1D->mData[p][q] -= singularValue*u->vData[p]*v->vData[q];
				}
			}
		}

		v_ = createVector(col);
		_svd_1d(matrixFor1D, v_);
		u_ = createVector(row);
		sigma2 = 0;
		for (p=0; p<row; p++){
			u_unnormalized_val = 0;
			for (q=0; q<col; q++){
				u_unnormalized_val += src[p*col+q]*v_->vData[q];
			}
			sigma2 += u_unnormalized_val*u_unnormalized_val;
			u_->vData[p] = u_unnormalized_val;
		}
		sigma = sqrt(sigma2);
		for (j=0; j<row; j++)
			u_->vData[j] /= sigma;
		svdSoFar[i].singularValue = sigma;
		svdSoFar[i].v = v_;
		svdSoFar[i].u = u_;
	}

	destroyMatrix(matrixFor1D);
	//unzip
	for (i=0; i<k; i++){
        //S
        svd.S->vData[i] = svdSoFar[i].singularValue;
        //U
        for (j=0; j<row; j++){
            svd.U->mData[j][i] = svdSoFar[i].u->vData[j];
        }
        destroyVector(svdSoFar[i].u);
        //V
        for (j=0; j<col; j++){
            svd.V->mData[i][j] = svdSoFar[i].v->vData[j];
        }
        destroyVector(svdSoFar[i].v);
	}
	return svd;
}

/** \brief interface to get S
 *
 * \param 2-dim array, row, col
 * \param 1-dim array, len
 * \return
 *
 */
void getS(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t len){
    SVD svd = _SVDdec(arr, row, col);
    Vector* s = svd.S;
    vectorToArray(s, dest, len);
    destroyMatrix(svd.U);destroyMatrix(svd.V);

}
/** \brief interface to get V
 *
 * \param 2-dim array, row, col
 * \param 2-dim array, row, col
 * \return
 *
 */
void getV(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    SVD svd = _SVDdec(arr, row, col);
    Matrix* v = svd.V;
    matrixToArray(v, dest, height, width);
    destroyMatrix(svd.U); destroyVector(svd.S);

}
/** \brief interface to get U
 *
 * \param 2-dim array, row, col
 * \param 2-dim array, row, col
 * \return
 *
 */
void getUs(INOUT float* arr, size_t row, size_t col, OUT float* dest, size_t height, size_t width){
    SVD svd = _SVDdec(arr, row, col);
    Matrix* u = svd.U;
    matrixToArray(u, dest, height, width);
    destroyMatrix(svd.V);destroyVector(svd.S);
}
