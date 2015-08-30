//*******************************************************************************************************
__device__ unsigned int dev_Multiplica_Matrizes(	float *matrizA,
																	int numLinhasMatrizA,
																	int numColunasMatrizA,
																	float *matrizB,
																	int numLinhasMatrizB,
																	int numColunasMatrizB,
																	float *matrizResultado) {
//*******************************************************************************************************

	int i, j, k;

	if(numColunasMatrizA != numLinhasMatrizB) {
		//printf("Soh pode haver multiplicacao se numColunasMatrizA = numLinhasMatrizB");
		return 1;
	}

	for (i = 0; i < numLinhasMatrizA; i++) {
		for (j = 0; j < numColunasMatrizB; j++) {
			matrizResultado[i*numColunasMatrizB + j] = 0.0;
			for (k = 0; k < numColunasMatrizA; k++) {
				*(matrizResultado + i*numColunasMatrizB + j) += (*(matrizA + i*numColunasMatrizA + k)) * (*(matrizB + k*numColunasMatrizB + j));
			}
		}
	}
	return 0;
}

//*****************************************************************************************
__device__ void dev_multiplica_matriz_por_escalar(	float *Matriz,
																	float Escalar,
																	unsigned int qtde_Elementos_da_Matriz,
																	float *Matriz_Resultado) {
//*****************************************************************************************

	unsigned int iElemento;

	for (iElemento = 0; iElemento < qtde_Elementos_da_Matriz; iElemento++) {
				*(Matriz_Resultado + iElemento) = *(Matriz + iElemento) * Escalar;
	}
}


//*****************************************************************************************
__device__ void dev_Subtrai_Matrizes(float *matrizA, float *matrizB, unsigned int numElementos, float *matrizResultado) {
//*****************************************************************************************

	unsigned int iElemento;

	for (iElemento = 0; iElemento < numElementos; iElemento++) {
				*(matrizResultado + iElemento) = *(matrizA + iElemento) - *(matrizB + iElemento);
	}
}

// As duas procedures abaixo são utilizadas para multiplicação de matrizes.

// Matrix multiplication - Host code
// Matrix dimensions are assumed to be multiples of BLOCK_SIZE


//*****************************************************************************************
__global__ void dev_Subtrai_Matrizes(float *dev_matrizA,
												 float *dev_matrizB,
												 unsigned int numLinhas,
												 unsigned int numColunas,
												 float *dev_Matriz_Resultado) {
//*****************************************************************************************

	unsigned int iLinha = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int iColuna = threadIdx.y + blockIdx.y*blockDim.y;

	if ( (iLinha < numLinhas) && (iColuna < numColunas)) {	
		unsigned int iElementoGlobal = iColuna + numColunas*iLinha;
		*(dev_Matriz_Resultado + iElementoGlobal) = *(dev_matrizA + iElementoGlobal) - *(dev_matrizB + iElementoGlobal);
	}
}


//*****************************************************************************************
__global__ void dev_multiplica_matriz_por_escalar_kernel(float *dev_Matriz,
																  unsigned int numLinhas,
																  unsigned int numColunas,
																  char flag_localEscalar,
																  float *dev_Escalar,
																  float host_Escalar,
																  float *dev_Matriz_Resultado) {
//*****************************************************************************************

	unsigned int iLinha = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int iColuna = threadIdx.y  + blockIdx.y*blockDim.y;
	
	if ( (iLinha < numLinhas) && (iColuna < numColunas)) {	
		float escalar;
		if (flag_localEscalar == 0) {
			// o escalar está na CPU/host
			escalar = host_Escalar;
		} 
		else {
			// o escalar está na GPU
			escalar = *dev_Escalar;
		}
		unsigned int iElementoGlobal = iColuna + numColunas*iLinha;
		*(dev_Matriz_Resultado + iElementoGlobal) = escalar * (*(dev_Matriz + iElementoGlobal));
	}
}

//*****************************************************************************************
void multiplica_matriz_por_escalar(	float *dev_Matriz,
												unsigned int numLinhas,
												unsigned int numColunas,
												char flag_localEscalar,
												float *dev_Escalar,
												float host_Escalar,
												float *dev_Matriz_Resultado) {
//*****************************************************************************************

	// TESTE - inicio.

	/*
	float *host_Matriz;
	unsigned int qtdeElementos = numLinhas*numColunas;
	host_Matriz = (float *)malloc(qtdeElementos*sizeof(float));
	cudaErrorCheck( cudaMemcpy(host_Matriz, dev_Matriz, qtdeElementos*sizeof(float), cudaMemcpyDeviceToHost ) );
	printf("\n"); printf("Escalar = %f", );
	printf("\n"); printf("Matriz antes = [");
	printf("\n"); imprimeMatriz(host_Matriz, numLinhas, numColunas);	
	printf("\n"); printf("]");
	*/

	// TESTE - fim.

	// flag_localEscalar;
	//		0: o escalar está na CPU
	//		1: o escalar está na GPU


	unsigned int numBlocosLinha = (numLinhas + nThreadsPorBloco)/nThreadsPorBloco;
	unsigned int numBlocosColuna = (numColunas + nThreadsPorBloco)/nThreadsPorBloco;

	dim3 numThreads(nThreadsPorBloco,nThreadsPorBloco);
	dim3 numBlocos(numBlocosLinha,numBlocosColuna);

	dev_multiplica_matriz_por_escalar_kernel<<<numBlocos,numThreads>>>(
			dev_Matriz,
			numLinhas,
			numColunas,
			flag_localEscalar,
			dev_Escalar,
			host_Escalar,
			dev_Matriz_Resultado);
	
	// TESTE - inicio.

	/*
	cudaErrorCheck( cudaMemcpy(host_Matriz, dev_Matriz, qtdeElementos*sizeof(float), cudaMemcpyDeviceToHost ) );
	printf("\n"); printf("Matriz depois = [");
	printf("\n"); imprimeMatriz(host_Matriz, numLinhas, numColunas);	
	printf("\n"); printf("]");
	host_Matriz = NULL; free(host_Matriz);
	*/

	// TESTE - FIM.

}

//+==============================================================================================
__global__ void MatMulKernel(	float *matrizA, unsigned int qtdeLinhasMatrizA, unsigned int qtdeColunasMatrizA,
										float *matrizB, unsigned int qtdeLinhasMatrizB, unsigned int qtdeColunasMatrizB,
										float *matrizC) {
//+==============================================================================================
	
	// Matrix multiplication kernel called by MatMul()

	// Each thread computes one element of C
	// by accumulating results into Cvalue

	float Cvalue = 0.0;
	int rowC = blockIdx.y * blockDim.y + threadIdx.y;
	int colC = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int	qtdeColunasMatrizC = qtdeColunasMatrizB;
	
	if ( ( rowC > (qtdeLinhasMatrizA - 1) ) || ( colC > (qtdeColunasMatrizB - 1) ) ) return;
	//if(row > A.height || col > B.width) return;
	
	for (int e = 0; e < qtdeColunasMatrizA; ++e)
		Cvalue += ( *(matrizA + rowC*qtdeColunasMatrizA + e) ) * ( *(matrizB + e*qtdeColunasMatrizB + colC) );
		//Cvalue += (A.elements[row * A.width + e]) * (B.elements[e * B.width + col]);
	
	*(matrizC + rowC * qtdeColunasMatrizC + colC) = Cvalue;
	//C.elements[row * C.width + col] = Cvalue;
//-------------------------------------------------------------------------------------
}
//-------------------------------------------------------------------------------------

//=====================================================================================
void MatMul(float *dev_matrizA, unsigned int numLinhasMatrizA, unsigned int numColunasMatrizA,
				float *dev_matrizB, unsigned int numLinhasMatrizB, unsigned int numColunasMatrizB,
				float *dev_matrizResultado, struct parametros *host_ParametrosGA) {
//=====================================================================================

	//dim3 dimBlock(nThreadsPorBloco, nThreadsPorBloco);
	dim3 dimBlock(256, 256);

	/*		A variável dimGrid dá o número de blocos, no grid, para cada dimensão. Note que o número de
		blocos é uma função do número de colunas da matriz B e das linhas da matriz A. Então, indiretamente,
		o número de threads é uma função do número de genes quando a matriz A ou B é o elemento 'gene' de 
		um indivíduo. Portanto, não há necessidade de alteração dessa função.
	*/

	// blocos por grid
	dim3 dimGrid((numColunasMatrizB + dimBlock.x - 1) / dimBlock.x,	(numLinhasMatrizA + dimBlock.y - 1) / dimBlock.y);
	
	MatMulKernel<<<dimGrid, dimBlock>>>(dev_matrizA, numLinhasMatrizA, numColunasMatrizA,
													dev_matrizB, numLinhasMatrizB, numColunasMatrizB,
													dev_matrizResultado);
	
	cudaError_t err = cudaThreadSynchronize();
	//printf("\n\n // Chamada do kernel: %s\n", cudaGetErrorString(err));
}

//===========================================================================================
__global__ void dev_GeraMatrizIdentidade(	float *dev_MatrizIdentidade,
														struct parametros *dev_parametrosGA) {
//===========================================================================================
	
	unsigned int iLinha = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int iColuna = threadIdx.y + blockIdx.y*blockDim.y;

	unsigned int numLinhas = dev_parametrosGA->numGenes;
	unsigned int numColunas = dev_parametrosGA->numGenes;
	float valorElemento;

	//for (iElemento = 0; iElemento < qtdeElementos; iElemento++) {
	if ((iLinha < numLinhas) && (iColuna < numColunas)) {
		if (iLinha == iColuna) {
			valorElemento = 1.0F;
		}
		else {
			valorElemento = 0.0F;
		}
		*(dev_MatrizIdentidade + iColuna + iLinha*numColunas) = valorElemento;
	}
}