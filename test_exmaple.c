#include "./include/base.h"
#include "./include/operation.h"
#include "./include/LUD.h"
#include "./include/QRD.h"
#include "./include/SVD.h"
#include "./include/FFT.h"

/*
This is an example, and test only whether it can run or not. The validity can be verified by Matlab-like software.
 */

int main(){
    float a1[3][3];
    Matrix* m1;Vector* v1;
    float a2[3];
    int i;float f;float sum=0;
    float a3[3]; float a4[2][2];
    float a5[3][3];
    float vn,mn;
    float a6[3][3] = {{1,2,3},{4,5,6},{7,8,9}};
    float a7[3] = {1,2,3};
    float Dot;
    float a8[8][3] = {{2,5,3},{1,2,1},{4,1,1},{3,5,2},{5,3,1},{4,5,5},{2,4,2},{2,2,5}};
    float a9[3];
    float a10[3][3];
    float a11[8][3];
    float a12[16] = {1,2,3,4,56,3,4,5,6,23,2,4,1,2,8,9};
    float a13[16];
	/*
	Test base
	 */
	m1 = createMatrix(3,3);
	printm(m1);
	assert(FLOAT_EQUY(m1->mData[1][0], 0));
	matrixToArrayWithoutFree(m1,a1,3,3);
	printmArr(a1,3,3);
	assert(FLOAT_EQUY(a1[1][0], 0));
    a1[0][0] = 1.1;a1[0][1] = 2.3; a1[1][1] = 5.2; a1[2][2] = 8;
    arrayToMatrix(a1,m1,3,3);
    printm(m1);
    assert(FLOAT_EQUY(m1->mData[0][1], 2.3));
    destroyMatrix(m1);
    v1 = createVector(3);
    printv(v1);
    assert((FLOAT_EQUY(v1->vData[1],0)));
    vectorToArrayWithoutFree(v1,a2,3);
    printvArr(a2,3);
    assert(FLOAT_EQUY(a2[1],0));
    a2[0] = 1.1; a2[2] = 2.3;
    arrayToVector(a2,v1,3);
    printv(v1);
    assert(FLOAT_EQUY(v1->vData[2], 2.3));
    destroyVector(v1);
    //some r.v. generation
    slach_rand_seed(0);
    for (i=0; i<100; i++){
        f = gaussRand(1,1);
        sum += f;
    }
    printf("sum of gaussian:%f\n\n",sum/100);

    /*
    Test operation
    */

    slicev(a1,3,3,1,2,0,2,a3,3);
    printvArr(a3,3);
    assert(FLOAT_EQUY(a3[2], 8));
    slicem(a1,3,3,1,2,1,2,a4,2,2);
    printmArr(a4,2,2);
    assert(FLOAT_EQUY(a4[1][1], 8));

    //Test element-wise math functions

    absv(a2,3,a3,3);
    printvArr(a3,3);
    assert(FLOAT_EQUY(a3[2], 2.3));
    sinv(a2,3,a3,3);
    printvArr(a3,3);
    cosv(a2,3,a3,3);
    printvArr(a3,3);
    tanv(a2,3,a3,3);
    printvArr(a3,3);
    expv(a2,3,a3,3);
    printvArr(a3,3);
    logv(a2,3,a3,3);
    printvArr(a3,3);
    powv(a2,3,3,a3,3);
    printvArr(a3,3);
    sqrtv(a2,3,a3,3);
    printvArr(a3,3);
    absm(a1,3,3,a5,3,3);
    printmArr(a5,3,3);
    assert(FLOAT_EQUY(a5[0][1], 2.3));
    sinm(a1,3,3,a5,3,3);
    printmArr(a5,3,3);
    cosm(a1,3,3,a5,3,3);
    printmArr(a5,3,3);
    tanm(a1,3,3,a5,3,3);
    printmArr(a5,3,3);
    expm(a1,3,3,a5,3,3);
    printmArr(a5,3,3);
    logm(a1,3,3,a5,3,3);
    printmArr(a5,3,3);
    powm(a1,3,3,3,a5,3,3);
    printmArr(a5,3,3);
    sqrtm(a1,3,3,a5,3,3);
    printmArr(a5,3,3);

    //norm
    vn = vNorm("2",a2,3);
    printf("vnorm: %f\n",vn);
    mn = mNorm("F",a1,3,3);
    printf("mnorm: %f\n",mn);

    //Test matrix multiply, add, transpose, vector dot,
    mmMul(a1,3,3,a6,3,3,a5,3,3);
    printmArr(a5,3,3);
    mvMul(a1,3,3,a2,3,1,a3,3);
    printvArr(a3,3);
    mmAdd(a1,3,3,a6,3,3,a5,3,3);
    printmArr(a5,3,3);
    vvAdd(a7,3,a2,3,a3,3);
    printvArr(a3,3);
    mT(a1,3,3,a5,3,3);
    printmArr(a5,3,3);
    printvArr(a7,3);printvArr(a2,3);
    Dot = dot(a7,3,a2,3);
	printf("dot: %f\n",Dot);


    /*
    Test LUD,  inverse
    */
    getL(a1,3,3,a5,3,3);
    printmArr(a5,3,3);
    getU(a1,3,3,a5,3,3);
    printmArr(a5,3,3);
    LUsolvev(a1,3,3,a2,3,a3,3);
    printvArr(a3,3);
    inv(a1,3,3,a5,3,3);
    printmArr(a5,3,3);

    /*
    Test QRD
    */
    getQ(a1,3,3,a5,3,3);
    printmArr(a5,3,3);
    getR(a1,3,3,a6,3,3);
    printmArr(a6,3,3);
    QRsolvev(a1,3,3,a2,3,a3,3);
    printvArr(a3,3);
    QRsolvem(a6,3,3,a1,3,3,a5,3,3);
    printmArr(a5,3,3);

    /*
    Test SVD
   */

    getS(a8,8,3,a9,3);
    printvArr(a9,3);
    getV(a8,8,3,a10,3,3);
    printmArr(a10,3,3);
    getUs(a8,8,3,a11,8,3);
    printmArr(a11,8,3);


    /*
    Test FFT
    */

    fftAbs(a12, 16,a13, 16,4,4);
    printvArr(a13,16);
	return 0;
}

