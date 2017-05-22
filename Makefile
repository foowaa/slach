all:
	$(CC) ./src/base.c ./src/operation.c ./src/LUD.c ./src/QRD.c ./src/SVD.c ./src/FFT.c test_example.c -o test_example -lm

test:
	 ./test_example || exit 1