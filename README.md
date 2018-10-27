# Simple Linear Algebra C Header (SLACH)

[![Build Status](https://travis-ci.org/foowaa/slach.svg?branch=master)](https://travis-ci.org/foowaa/slach)

slach is a simple linear algebra libraries desigined for *low-coupling* used C89, and it is convient for portability. slach has main parts: base, operation, QR decomposition, LU decomposition, SVD decomposition and Fast Fourier Transform(FFT). All other parts except for base all depend on base, rather than compelx dependency. If users want to just use the library, the users don't need to know inner data structures, instead you can use arraies to achieve results. If users want to implement the algorithms, it's easy to transfer programs. In the program, I seperate the inner and outer clearly, so a lot of wrapper functions employed in the program. 

Usage
------
The interface of slach basically like: 

* `float* src`
* `size_t row`, `size_t col` OR `size_t len`
* `float* dest`
* `size_t height`, `size_t width` OR `size_t len` 
* other paramaters

The meaning is: use this function on `src` and save in `dest`.

base
------
base declares and defines basic data structures: vector and matrix, it can be easily used in other applications. base also defines some utilities: safe malloc and free, print function and random numbers generation. 

1. `Matrix` and `Vector` provide some basic functions: create, destroy, deep copy, array to matrix, matrix to array, vector to array, array to vector.
2. `slach_malloc` and `slach_free` are safe memory control functions.
3. `slach_rand_seed` sets rand seed, `slach_rand_int_range_*` generates integer r.v. in different range, `uRand` generates uniform distribution, `gaussrand` generates Gaussian distribution, `expRand` generates exponential distribution.
4. `perr` print error and exit program, `print*` print vectors and matrices.

operation
------
operation declares and defines operation functions: element-wise math function, matrix multiplication, add, transpose, vector inner product, vector l-p norm and matrix norm, slice like matlab.

1. `*m` and `*v` are element-wise math functions.
2. `slicev` and `slicem` do slice like matlab.
3. `mmMul`, `mvMul`, `mmAdd`, `vvAdd`, `dot`, `vnorm` and `mnorm` do matrix multiplication, add, transpose, vector inner product, vector l-p norm and matrix norm.

LUD
------
LUD implements LU decomposition, and offers LUD interface, solving linear equations using LUD and matrix inverse using LUD.

1. `getL` and `getU` get `L` and `U` matrix.
2. `LUsolvem` and `LUsolvev` solve linear equations: AX=b, AX=B
3. `inv` inverse square matrix using LU decomposition.

QRD
------
QRD implements QR decomposition, and offers QRD interface, solving linear equations using QRD.

1. `getQ` and `getR` get `Q` and `R` matrix.
2. `QRsolvem` and `QRsolvev` solve linear equations: AX=b, AX=B

SVD
------
SVD implements SVD decomposition, and offers SVD interface. And also we can do eigenvalue decomposition using SVD.

`getS`, `getV` and `getUs` get `S` vector, `V` and `U` matrix.

FFT
-----
FFT implements naive *Discrete Fourier Transform* and *Cooley-Turkey FFT*. Besides, we provide abs and phase using FFT--often we use in reality is abs and pahse after FFT. And we also provide `DFT_naive` and `FFT_CooleyTukey`.

1. `fftAbs` and `fftPhase` do FFT and then calculate abs and pahse.
2. `DFT_naive` and `FFT_CooleyTukey` return complex value of FFT.
3. `fftshift` do like matlab `fftshift`.
