#makefile
all:
	$(CC) -Wall -Wextra ./include/base.h ./include/operation.h ./include/LUD.h ./include/QRD.h ./include/SVD.h ./include/FFT.h ./src/base.c ./src/operation.c ./src/LUD.c ./src/QRD.c ./src/SVD.c ./src/FFT.c test_example.c -o test_example

test:
	 ./test_example || exit 1