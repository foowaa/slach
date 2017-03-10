#include "slach.h"

Matrix* createMatrix(IN unsigned int mHeight, IN unsigned int mWidth){
	Matrix* mPtr;
	int i;

	if (mHeight == 0 || mWidth == 0){
		return NULL;
	}
	else{
		mPtr = slach_malloc (Matrix, Matrix);
		mPtr->mData = slach_malloc(T*, mHeight);
		for (i = 0; i<mHeight; i++){
			mPtr->mData[i] = slach_malloc(T, mWidth);
		}

		mPtr->mHeight = mHeight;
		mPtr->mWidth = mWidth;
	}
}


void destroyMatrix(INOUT Matrix* mPtr){
	int i;
	if (mPtr == NULL){
		return;
	}
	else{
		for (i = 0; i<mPtr->mHeight; i++){
			slach_free(mPtr->mData[i]);
		}

		slach_free(mPtr->mData);
		slach_free(mPtr);
	}
}

//deep copy
void copyMatrix(IN Matrix* src, OUT Matrix* dest){
	int i;
	if (src == NULL || dest == NULL){
		return;
	}
	else if (src->mHeight != dest->mheight || src->mWidth != dest->mWidth){
		return;
	}
	else{
		for (i=0; i<mHeight; i++){
			memcpy(src->mData[i], dest->mData[i], mWidth); 
		}
	}
}

//自行确保 array和matrix的维度相同
void arrayToMatrix(IN T** src, OUT Matrix* dest, unsigned int height, unsigned int weight){
	int i,j;
	if (src == NULL || dest == NULL){
		return;
	}
	else if (height != dest->mHeight || width != dest->mHeight){
		return;
	}
	else{
		for (i=0; i<height; i++){
			for (j=0; j<width; j++){
				*(dest->mData[i]+j) = src[i][j];
			}
		}
	}
}

void matrixToArray(IN Matrix* src, OUT T** dest, unsigned int height, unsigned int weight){
	int i,j;
	if (src == NULL || dest == NULL){
		return;
	}
	else if (height != dest->mHeight || width != dest->mHeight){
		return;
	}
	else{
		for (i=0; i<height; i++){
			for (j=0; j<width; j++){
				*(dest->mData[i]+j) = src[i][j];
			}
		}
	}
}


