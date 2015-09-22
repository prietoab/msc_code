
//*********************************************************************************************************
// Cabeçalhos
//*********************************************************************************************************

void Multiplica_Matrizes(float *matrizA, unsigned short int numLinhasMatrizA, unsigned short int numColunasMatrizA, float *matrizB, unsigned short int numLinhasMatrizB, unsigned short int numColunasMatrizB, float *matrizResultado);
void testa_Multiplica_Matrizes(void);

void multiplica_matriz_por_escalar(float *Matriz, float Escalar, unsigned long int qtde_Elementos_da_Matriz, float *Matriz_Resultado);
void testa_multiplica_matriz_por_escalar(void);

void Subtrai_Matrizes(float *matrizA, float *matrizB, unsigned long int numElementos, float *matrizResultado);
void testa_Subtrai_Matrizes(void);

void imprimeMatriz(float *matriz, unsigned short int numLinhas, unsigned short int numColunas);
void testa_imprimeMatriz(void);


//*********************************************************************************************************
// Funções
//*********************************************************************************************************

//=========================================================================================
void Multiplica_Matrizes(float *matrizA,
								 unsigned short int numLinhasMatrizA,
								 unsigned short int numColunasMatrizA,
								 float *matrizB,
								 unsigned short int numLinhasMatrizB,
								 unsigned short int numColunasMatrizB,
								 float *matrizResultado) {
//=========================================================================================

	unsigned short int i, j, k;

	if(numColunasMatrizA != numLinhasMatrizB) {
		printf("Soh pode haver multiplicacao se numColunasMatrizA = numLinhasMatrizB");
		return;
	}

	for (i = 0; i < numLinhasMatrizA; i++) {
		for (j = 0; j < numColunasMatrizB; j++) {
			matrizResultado[i*numColunasMatrizB + j] = 0.0;
			for (k = 0; k < numColunasMatrizA; k++) {
				*(matrizResultado + i*numColunasMatrizB + j) += (*(matrizA + i*numColunasMatrizA + k)) * (*(matrizB + k*numColunasMatrizB + j));
			}
		}
	}
}

//=================================================================================
void imprimeMatriz(float *matriz,
						 unsigned short int numLinhas,
						 unsigned short int numColunas) {
//=================================================================================

	unsigned short int iLinha, iColuna;

	for (iLinha = 0; iLinha < numLinhas; iLinha++) {
		printf("\n");
		printf("%f", *(matriz + 0 + numColunas*iLinha));
		for (iColuna = 1; iColuna < numColunas; iColuna++) {
			printf("\t"); printf("%f", *(matriz + iColuna + numColunas*iLinha));
		}
	}
}

//=================================================================================
void testa_imprimeMatriz(void) {
//=================================================================================

	unsigned short int numLinhas, numColunas;
	unsigned long int numElementos, iElemento;

	float *matriz;

	numLinhas = 132;
	numColunas = 3;

	numElementos = numLinhas * numColunas;

	matriz = (float *)malloc(numElementos*sizeof(float));

	for (iElemento = 0; iElemento < numElementos; iElemento++) {
		*(matriz + iElemento) = (float)iElemento;
	}
	
	printf("matriz =");
	imprimeMatriz(matriz, numLinhas, numColunas);

}

//=================================================================================
void testa_Multiplica_Matrizes(void) {
//=================================================================================
	unsigned short int numLinhasA = 100;
	unsigned short int numColunasA = 200;

	unsigned short int numLinhasB = numColunasA;
	unsigned int numColunasB = 300;

	unsigned short int numLinhasC = numLinhasA;
	unsigned short int numColunasC = numColunasB;

	float *matrizA, *matrizB, *matrizC;

	unsigned long int	iElemento, qtdeElementos;

	srand((unsigned int)time (NULL));

	unsigned short int parametro = 1; // RAND_MAX;

	matrizA = (float *)malloc(numLinhasA*numColunasA*sizeof(float));
	matrizB = (float *)malloc(numLinhasB*numColunasB*sizeof(float));
	matrizC = (float *)malloc(numLinhasC*numColunasC*sizeof(float));

	qtdeElementos = numLinhasA*numColunasA;
	for (iElemento = 0; iElemento < qtdeElementos; iElemento++) {
		*(matrizA + iElemento) = (float)rand()/parametro; // /;
	}
	
	printf("\n"); printf("A = [");
	printf("\n"); imprimeMatriz(matrizA, numLinhasA, numColunasA);
	printf("];");

	qtdeElementos = numLinhasB*numColunasB;
	for (iElemento =0; iElemento < qtdeElementos; iElemento++) {
		*(matrizB + iElemento) = (float)rand()/parametro; // /RAND_MAX;
	}

	printf("\n"); printf("B = [");
	printf("\n"); imprimeMatriz(matrizB, numLinhasB, numColunasB);
	printf("];");

	Multiplica_Matrizes(	matrizA, numLinhasA, numColunasA,
								matrizB, numLinhasB, numColunasB,
								matrizC);

	printf("\n"); printf("C = [");
	printf("\n"); imprimeMatriz(matrizC, numLinhasC, numColunasC);
	printf("];");

	printf("\n"); printf("C_Scilab = A*B;");
	printf("\n"); printf("diferenca = C - C_Scilab;"); 
	printf("\n"); printf("diferenca"); 
	printf("\n"); printf("max(diferenca)");
	printf("\n"); printf("min(diferenca)");

	printf("\n"); printf("max(diferenca) / min(C_Scilab)");
	printf("\n"); printf("(max(diferenca) / min(C_Scilab))*100");

	free(matrizA);
	free(matrizB);
	free(matrizC);
}
//*****************************************************************************************
void multiplica_matriz_por_escalar(float *Matriz, float Escalar, unsigned long int qtde_Elementos_da_Matriz, float *Matriz_Resultado) {
//*****************************************************************************************

	unsigned long int iElemento;

	for (iElemento = 0; iElemento < qtde_Elementos_da_Matriz; iElemento++) {
				*(Matriz_Resultado + iElemento) = *(Matriz + iElemento) * Escalar;
	}
}

//=================================================================================
void testa_multiplica_matriz_por_escalar(void) {
//=================================================================================
	
	unsigned short int numLinhasA = 10;
	unsigned short int numColunasA = 20;
	
	float *matrizA, escalar = 2.4F;

	unsigned long int	iElemento, qtdeElementos;

	srand((unsigned int)time (NULL));

	unsigned short int parametro = 1; // RAND_MAX;

	matrizA = (float *)malloc(numLinhasA*numColunasA*sizeof(float));

	qtdeElementos = numLinhasA*numColunasA;
	for (iElemento = 0; iElemento < qtdeElementos; iElemento++) {
		*(matrizA + iElemento) = (float)rand()/parametro; // /;
	}
	
	
	printf("\n"); printf("A_antes = [");
	printf("\n"); imprimeMatriz(matrizA, numLinhasA, numColunasA);
	printf("];");

	printf("\n"); printf("escalar = %f", escalar);

	multiplica_matriz_por_escalar(matrizA,
											escalar,
											qtdeElementos,
											matrizA);

	printf("\n"); printf("A_depois = [");
	printf("\n"); imprimeMatriz(matrizA, numLinhasA, numColunasA);
	printf("];");

	printf("\n"); printf("A_depois_Scilab = A_antes * escalar;");
	printf("\n"); printf("diferenca = A_depois - A_depois_Scilab;"); 
	printf("\n"); printf("diferenca"); 
	printf("\n"); printf("max(diferenca)");
	printf("\n"); printf("min(diferenca)");

	printf("\n"); printf("max(diferenca) / min(A_depois_Scilab)");
	printf("\n"); printf("(max(diferenca) / min(A_depois_Scilab))*100");

	free(matrizA);
}

//=================================================================================
void Subtrai_Matrizes(float *matrizA, float *matrizB, unsigned long int numElementos, float *matrizResultado) {
//=================================================================================

	unsigned long int iElemento;

	for (iElemento = 0; iElemento < numElementos; iElemento++) {
				*(matrizResultado + iElemento) = *(matrizA + iElemento) - *(matrizB + iElemento);
	}
}


//=================================================================================
void Subtrai_Matrizes_Geral(float *C, float a, float *A, float b, float *B, unsigned long int numElementos) {
//=================================================================================

	/*
		Obtém a seguinte subtração de matrizes:

			C = a*A - b*B
	*/

	unsigned long int iElemento;

	for (iElemento = 0; iElemento < numElementos; iElemento++) {
		*(C + iElemento) = a*(*(A + iElemento)) - b*(*(B + iElemento));
	}
}

//=================================================================================
void testa_Subtrai_Matrizes(void) {
//=================================================================================
	
	unsigned short int numLinhasA = 5;
	unsigned short int numColunasA = 7;

	unsigned short int numLinhasB = numLinhasA;
	unsigned short int numColunasB = numColunasA;

	unsigned short int numLinhasC = numLinhasA;
	unsigned short int numColunasC = numColunasA;

	float *matrizA, *matrizB, *matrizC;

	unsigned long int	iElemento, qtdeElementos = numLinhasA * numColunasA;

	srand((unsigned int)time (NULL));

	unsigned short int parametro = 1; // RAND_MAX;

	matrizA = (float *)malloc(numLinhasA*numColunasA*sizeof(float));
	matrizB = (float *)malloc(numLinhasB*numColunasB*sizeof(float));
	matrizC = (float *)malloc(numLinhasC*numColunasC*sizeof(float));

	for (iElemento = 0; iElemento < qtdeElementos; iElemento++) {
		*(matrizA + iElemento) = (float)rand()/parametro; // /;
	}
	
	printf("\n"); printf("A = [");
	printf("\n"); imprimeMatriz(matrizA, numLinhasA, numColunasA);
	printf("];");

	qtdeElementos = numLinhasB*numColunasB;
	for (iElemento =0; iElemento < qtdeElementos; iElemento++) {
		*(matrizB + iElemento) = (float)rand()/parametro; // /RAND_MAX;
	}

	printf("\n"); printf("B = [");
	printf("\n"); imprimeMatriz(matrizB, numLinhasB, numColunasB);
	printf("];");

	qtdeElementos = numLinhasA * numColunasA;

	Subtrai_Matrizes(matrizA, matrizB, qtdeElementos, matrizC);

	printf("\n"); printf("C = [");
	printf("\n"); imprimeMatriz(matrizC, numLinhasC, numColunasC);
	printf("];");

	printf("\n"); printf("C_Scilab = A - B;");
	printf("\n"); printf("diferenca = C - C_Scilab;"); 
	printf("\n"); printf("diferenca"); 
	printf("\n"); printf("max(diferenca)");
	printf("\n"); printf("min(diferenca)");

	printf("\n"); printf("max(diferenca) / min(C_Scilab)");
	printf("\n"); printf("(max(diferenca) / min(C_Scilab))*100");

	free(matrizA);
	free(matrizB);
	free(matrizC);
}

//===========================================================================================
void GeraMatrizIdentidade(	float *host_MatrizIdentidade,
									unsigned short int ordem_da_matriz) {
//===========================================================================================
	
	unsigned short int iLinha, iColuna;
	unsigned short int numLinhas = ordem_da_matriz;
	unsigned short int numColunas = ordem_da_matriz;
	float valorElemento;

	for (iLinha = 0; iLinha < numLinhas; iLinha++) {

		for (iColuna = 0; iColuna < numColunas; iColuna++) {

			if (iLinha == iColuna) {
				valorElemento = 1.0F;
			}
			else {
				valorElemento = 0.0F;
			}

			*(host_MatrizIdentidade + iColuna + iLinha*numColunas) = valorElemento;

		}
	}
}


//===========================================================================================
void testaMatrizIdentidade(void) {
//===========================================================================================

	unsigned short int ordemMatriz;
	float *matrizI;

	for (ordemMatriz = 100; ordemMatriz < 1000; ordemMatriz += 100) {
		matrizI = (float *)malloc(ordemMatriz*ordemMatriz*sizeof(float));
		GeraMatrizIdentidade(matrizI, ordemMatriz);
		printf("\n"); printf("I (ordem %d) = ", ordemMatriz);
		imprimeMatriz(matrizI, ordemMatriz, ordemMatriz);
		free(matrizI);
	}
}