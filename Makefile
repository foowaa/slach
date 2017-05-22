all:
	$(CC) -lm ./src/base.c ./src/operation.c ./src/LUD.c ./src/QRD.c ./src/SVD.c ./src/FFT.c test_example.c -o test_example

test:
	 ./test_example || exit 1