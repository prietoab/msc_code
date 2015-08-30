/*
* multNoShare.h
*
* Robert Hochberg
* January 24, 2012
*
* Based nearly entirely on the code from the CUDA C Programming Guide
*/

// Matrices are stored in row-major order:
// M(row, col) = *(M.elements + row * M.width + col)
/*
typedef struct {
	int width;
	int height;
	float* elements;
} Matrix;
*/

//__global__ void MatMulKernel(Matrix A, Matrix B, Matrix C);

__global__ void MatMulKernel(	float *matrizA, unsigned int qtdeLinhasMatrizA, unsigned int qtdeColunasMatrizA,
										float *matrizB, unsigned int qtdeLinhasMatrizB, unsigned int qtdeColunasMatrizB,
										float *matrizC);