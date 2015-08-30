#include <stdio.h> // biblioteca input/output padrão do C 
#include <stdlib.h>
#include <time.h>
#include <malloc.h>
#include <math.h>
#include <curand_kernel.h>

#define nThreadsPorBloco 16

#include "C:\Users\Comum\Documents\Unicamp\Mestrado\Programação\_bibliotecas\Estruturas.h"
#include "C:\Users\Comum\Documents\Unicamp\Mestrado\Programação\_bibliotecas\Algebra_Linear_serial.h"
#include "C:\Users\Comum\Documents\Unicamp\Mestrado\Programação\_bibliotecas\Algebra_Linear_paralelo.h"
#include "C:\Users\Comum\Documents\Unicamp\Mestrado\Programação\_bibliotecas\Auxiliares_serial.h"
#include "C:\Users\Comum\Documents\Unicamp\Mestrado\Programação\_bibliotecas\Auxiliares_paralelo.h"
#include "C:\Users\Comum\Documents\Unicamp\Mestrado\Programação\_bibliotecas\GA_Paralelo.h"


//---------------------------------------------------------------------------------------
int main (void) {
//---------------------------------------------------------------------------------------
	
	clock_t			time_i = 999, time_f = 999
					,	time_pgm_i = 999, time_pgm_f = 999;
	
	// Obtem o momento inicial
	time_pgm_i = clock();
		
	struct parametros				parametrosGA, *devParametrosGA;
	struct parametros_Metodo	parametrosMetodo, *devParametrosMetodo;
	struct generation				geracao[2], *devGeracao0, *devGeracao1;
	unsigned int					iGeracao, qtdeElementos;

	cudaErrorCheck(cudaMalloc((void **)&devGeracao0, sizeof(struct generation)));
	cudaErrorCheck(cudaMalloc((void **)&devGeracao1, sizeof(struct generation)));

	inicializa_Parametros(&parametrosGA, &parametrosMetodo);
	inicializa_Semente();

	// imprime o cabeçalho
	//imprimeTempo(0, 0, 0 , 0, 1, 2, time_i, time_f);
	
	// -----------------------------------------------------------------------
	// Copia para a GPU os parâmetros do GA e do Método.
	// -----------------------------------------------------------------------
	cudaErrorCheck(cudaMalloc((void **)&devParametrosGA, sizeof(struct parametros)));
	cudaErrorCheck(cudaMemcpy(devParametrosGA, &parametrosGA, sizeof(struct parametros), cudaMemcpyHostToDevice));

	cudaErrorCheck(cudaMalloc((void **)&devParametrosMetodo, sizeof(struct parametros_Metodo)));
	cudaErrorCheck(cudaMemcpy(devParametrosMetodo, &parametrosMetodo, sizeof(struct parametros_Metodo), cudaMemcpyHostToDevice));

	// -----------------------------------------------------------------------
	// HAMILTONIANO
	// -----------------------------------------------------------------------
	time_i = clock();
	float *hostHamiltoniano, *devHamiltoniano;
	hostHamiltoniano = (float *)malloc(parametrosGA.numGenes * parametrosGA.numGenes*sizeof(float));
	cudaErrorCheck(cudaMalloc((void **)&devHamiltoniano, parametrosGA.numGenes * parametrosGA.numGenes * sizeof(float)));
	unsigned int numBlocosLinha, numBlocosColuna;
	numBlocosLinha = (parametrosGA.numGenes + nThreadsPorBloco)/nThreadsPorBloco;
	numBlocosColuna = (parametrosGA.numGenes + nThreadsPorBloco)/nThreadsPorBloco;
	dim3 ThreadsH(nThreadsPorBloco, nThreadsPorBloco);
	dim3 BlocosH(numBlocosLinha, numBlocosColuna);
	dev_Gera_Matriz_de_Coope_naGPU<<<BlocosH,ThreadsH>>>(devHamiltoniano, devParametrosGA);
	time_f = clock();
	//imprimeTempo(1, 0, 1 , 0, 1, 2, time_i, time_f);
	//testaHamiltoniano(hostHamiltoniano, devHamiltoniano, &parametrosGA);

	// no link abaixo há informações sobre a aparente não execução da função acima
	// http://stackoverflow.com/questions/21982826/cuda-kernel-function-not-called
	// http://stackoverflow.com/questions/11587227/basic-cuda-getting-kernels-to-run-on-the-device-using-c
	
	// medindo os erros:
	// http://devblogs.nvidia.com/parallelforall/how-query-device-properties-and-handle-errors-cuda-cc/	
	
	
	// -----------------------------------------------------------------------
	// NÚMEROS PSEUDO-ALEATÓRIOS
	// -----------------------------------------------------------------------
	// semente para os números aleatórios
	unsigned int hostSemente, *devSemente;
	hostSemente = (unsigned int)time(NULL);
	cudaErrorCheck(cudaMalloc((void **)&devSemente, sizeof(unsigned int)));
	cudaErrorCheck(cudaMemcpy(devSemente, &hostSemente, sizeof(unsigned int), cudaMemcpyHostToDevice));
	
	// Inicia estado para números pseudoaleatrórios.
	curandState *dev_estadosCURAND;
	unsigned int numEstadosTotal = parametrosGA.numIndividuos*parametrosGA.numGenes;
	cudaErrorCheck(cudaMalloc((void **)&dev_estadosCURAND, numEstadosTotal*sizeof(curandState)));

	unsigned int numBlocosCuRand = (numEstadosTotal + nThreadsPorBloco)/nThreadsPorBloco;

	dim3 numThreadsCURAND(nThreadsPorBloco);
	dim3 numBlocosCURAND(numBlocosCuRand);

	devInicializaCURAND<<<numBlocosCURAND,numThreadsCURAND>>>(	devGeracao0,
																					devSemente,
																					dev_estadosCURAND,
																					devParametrosGA);
	cudaErrorCheck(cudaFree(devSemente));
	
	// -----------------------------------------------------------------------
	// POPULACAO INICIAL
	// -----------------------------------------------------------------------
	time_i = clock();
	GeraPopulacaoInicial_paralela(devGeracao0, &parametrosGA, devParametrosGA, dev_estadosCURAND);
	time_f = clock();
	//imprimeTempo(1, 0, 1, 0, 2, 2, time_i, time_f);
	//testePopulacaoInicial(&geracao[0], devGeracao0, &parametrosGA);
	
	// -----------------------------------------------------------------------
	// MATRIZ IDENTIDADE
	// -----------------------------------------------------------------------
	time_i = clock();
	float *devMatrizIdentidade;
	qtdeElementos = parametrosGA.numGenes * parametrosGA.numGenes;
	cudaErrorCheck(cudaMalloc((void **)&devMatrizIdentidade, qtdeElementos*sizeof(float)));
	numBlocosLinha = (parametrosGA.numGenes + nThreadsPorBloco)/nThreadsPorBloco;
	numBlocosColuna = (parametrosGA.numGenes + nThreadsPorBloco)/nThreadsPorBloco;
	dim3 ThreadsI(nThreadsPorBloco,nThreadsPorBloco);
	dim3 BlocosI(numBlocosLinha,numBlocosColuna);
	dev_GeraMatrizIdentidade<<<BlocosI,ThreadsI>>>(devMatrizIdentidade, devParametrosGA);
	time_f = clock();
	imprimeTempo(1, 0, 1 , 0, 3, 2, time_i, time_f);
	//testeMatrizIdentidade_GPU(devMatrizIdentidade,&parametrosGA);

	// -----------------------------------------------------------------------
	// ITERAÇÃO PRINCIPAL
	// -----------------------------------------------------------------------

	// - teste selecao - pode ser removido depois
	struct testeGeracao_s *dev_testeGeracao_s;
	cudaErrorCheck(cudaMalloc((void **)&dev_testeGeracao_s, sizeof(testeGeracao_s)));
	//
	imprimeComportamentoFitness_GPU(0, 0, 0, 0, 0, devGeracao0, devParametrosMetodo);

	for (iGeracao = 0; iGeracao < parametrosGA.numGeracoes; iGeracao++) {
		
		// -----------------------------------------------------------------------
		// FITNESS -- INICIO
		// -----------------------------------------------------------------------

		time_i = clock();	
		Fitness(	'2',
					devHamiltoniano,
					devGeracao0,
					&parametrosGA,
					devParametrosGA,
					devParametrosMetodo,
					devMatrizIdentidade);

		//teste_Fitness_GPU(devGeracao0, &parametrosGA);
		time_f = clock();
		//imprimeTempo(1, 0, 1, iGeracao, 8, 2, time_i, time_f);
		imprimeComportamentoFitness_GPU(1, 0, 1, 2, iGeracao, devGeracao0, devParametrosMetodo);


		// -----------------------------------------------------------------------
		// SELECAO -- INICIO
		// -----------------------------------------------------------------------

		time_i = clock();
		//testeSelecao(0, devGeracao0, &parametrosGA, dev_testeGeracao_s);
		
		Selecao(	devGeracao0,
					devGeracao1,
					&parametrosGA,
					devParametrosGA,
					dev_estadosCURAND,
					dev_testeGeracao_s);
		
		//testeSelecao(1, devGeracao1, &parametrosGA, dev_testeGeracao_s);
		time_f = clock();		
		//imprimeTempo(1, 0, 1, iGeracao, 5, 2, time_i, time_f);
		
		// -----------------------------------------------------------------------
		// CROSSOVER - INICIO
		// -----------------------------------------------------------------------

		time_i = clock();
		//testeCrossOver(2, devGeracao1, &parametrosGA, dev_testeGeracao_s);
		//testeCrossOver(0, devGeracao1, &parametrosGA, dev_testeGeracao_s);
		CrossOver1Ponto(
			devGeracao1,
			devGeracao0,
			&parametrosGA,
			devParametrosGA,
			dev_estadosCURAND,
			dev_testeGeracao_s);
		//testeCrossOver(1, devGeracao0, &parametrosGA, dev_testeGeracao_s);
		time_f = clock();
		//imprimeTempo(1, 0, 1, iGeracao, 6, 2, time_i, time_f);

		// -----------------------------------------------------------------------
		// MUTACAO - INICIO
		// -----------------------------------------------------------------------

		time_i = clock();
		//testeMutacao(0, devGeracao0, &parametrosGA);
		Mutacao(
			devGeracao0,
			&parametrosGA,
			devParametrosGA,
			dev_estadosCURAND);
		//testeMutacao(1, devGeracao0, &parametrosGA);
		time_f = clock();
		//imprimeTempo(1, 0, 1, iGeracao, 7, 2, time_i, time_f);		
				
	} // ===> Fim da iteração principal
						
	
	time_i = clock();
	Fitness(	'2',
				devHamiltoniano,
				devGeracao0,
				&parametrosGA,
				devParametrosGA,
				devParametrosMetodo,
				devMatrizIdentidade);

	time_f = clock();
	//imprimeTempo(1, 0, 1, iGeracao, 8, 2, time_i, time_f);
	imprimeComportamentoFitness_GPU(1, 0, 1, 2, iGeracao, devGeracao0, devParametrosMetodo);


	// Obtem o momento final em segundos e imprime a duração total do programa
	time_pgm_f = clock();
	imprimeTempo(1, 0, 1, iGeracao, 0, 2, time_pgm_i, time_pgm_f);

	// Retorno dos valores da GPU. Aqui
	printf("\n\n");
	printf("--------------------------------------------------------\n");
	printf("Geracao %d final\n", iGeracao);
	printf("--------------------------------------------------------\n");

	cudaErrorCheck(cudaMemcpy(&geracao[0], devGeracao0, sizeof(struct generation), cudaMemcpyDeviceToHost) );
	imprimeGeracao(&geracao[0], &parametrosGA);

	hostHamiltoniano = NULL; free(hostHamiltoniano);
	cudaErrorCheck(cudaFree(devHamiltoniano));
	cudaErrorCheck(cudaFree(devGeracao0));
	cudaErrorCheck(cudaFree(devGeracao1));
	cudaErrorCheck(cudaFree(devParametrosGA));
	cudaErrorCheck(cudaFree(devParametrosMetodo));
	cudaErrorCheck(cudaFree(devMatrizIdentidade));
	//cudaErrorCheck(cudaFree(dev_testeGeracao_s));

	return 0;
}