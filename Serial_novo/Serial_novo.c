#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "C:\Users\Comum\Documents\Unicamp\Mestrado\Programação\_bibliotecas\Estruturas.h"
#include "C:\Users\Comum\Documents\Unicamp\Mestrado\Programação\_bibliotecas\Auxiliares_serial.h"
#include "C:\Users\Comum\Documents\Unicamp\Mestrado\Programação\_bibliotecas\Algebra_Linear_serial.h"
#include "C:\Users\Comum\Documents\Unicamp\Mestrado\Programação\_bibliotecas\GA_Serial.h"
#include "C:\Users\Comum\Documents\Unicamp\Mestrado\Programação\_bibliotecas\Estatistica.h"

//==========================================================================
int main(int argc, char *argv[]) {
//==========================================================================

	// Pega os parâmetros da linha de comando
	parametrosPrograma parmsPrograma;
	// coloquei na mão só pra testar a alocação de memória
	parmsPrograma.parmQtdeGenes				= 50;
	parmsPrograma.parmQtdeMaxGeracoes		= 500000;
	parmsPrograma.parmQtdeIndividuos			= 100; // número de indivíduos por geração
	parmsPrograma.parmTamanhoTorneio			= 2;
	parmsPrograma.parmProbCrossOver			= 0.8F;
	parmsPrograma.parmQtdePontosCorte		= 2;
	parmsPrograma.parmProbMutacao				= 0.3F;
	parmsPrograma.parmIntensidadeMutacao	= 0.2F;
	parmsPrograma.parmLambda					= 0.005F;
	parmsPrograma.parmRhoMinimo				= 0.2F;
	parmsPrograma.parmTolerancia				= 0.00001F;

	
	parmsPrograma.parmQtdeGenes				= atoi(argv[1]);
	parmsPrograma.parmQtdeMaxGeracoes		= atoi(argv[2]);
	parmsPrograma.parmQtdeIndividuos			= atoi(argv[3]); // número de indivíduos por geração
	parmsPrograma.parmTamanhoTorneio			= atoi(argv[4]);
	parmsPrograma.parmProbCrossOver			= (float)atof(argv[5])/100;
	parmsPrograma.parmQtdePontosCorte		= atoi(argv[6]);
	parmsPrograma.parmProbMutacao				= (float)atof(argv[7])/100;
	parmsPrograma.parmIntensidadeMutacao	= (float)atof(argv[8])/10;
	parmsPrograma.parmLambda					= (float)atof(argv[9]);
	parmsPrograma.parmRhoMinimo				= (float)atof(argv[10]);
	parmsPrograma.parmTolerancia				= (float)atof(argv[11]);
	

	// Parâmetros variáveis
	printf("parmsPrograma.parmQtdeGenes = %d	\n", parmsPrograma.parmQtdeGenes);
	printf("parmsPrograma.parmQtdeIndividuos = %d	\n", parmsPrograma.parmQtdeIndividuos);
	printf("parmsPrograma.parmProbCrossOver = %f	\n", parmsPrograma.parmProbCrossOver);
	printf("parmsPrograma.parmProbMutacao = %f	\n", parmsPrograma.parmProbMutacao);
	printf("parmsPrograma.parmIntensidadeMutacao = %f	\n", parmsPrograma.parmIntensidadeMutacao);

	// Parâmetros fixos
	printf("\n");
	printf("parmsPrograma.parmQtdeMaxGeracoes = %d	\n", parmsPrograma.parmQtdeMaxGeracoes);
	printf("parmsPrograma.parmTamanhoTorneio = %d	\n", parmsPrograma.parmTamanhoTorneio);
	printf("parmsPrograma.parmQtdePontosCorte = %d	\n", parmsPrograma.parmQtdePontosCorte);
	printf("parmsPrograma.parmLambda = %f	\n", parmsPrograma.parmLambda);
	printf("parmsPrograma.parmRhoMinimo = %f	\n", parmsPrograma.parmRhoMinimo);

	// Definição das variáveis utilizadas para marcar o tempo.
	// São do tipo clock_t, pois marcarei diretamente o clock
	// absoluto do processador a partir do início do programa.
	clock_t	time_pgm_i = 0, time_pgm_f = 0,
				time_i = 0, time_f = 0;
			
	// Obtem o momento inicial do programa
	time_pgm_i = clock();

	// Definição das estruturas e variáveis globais.
	struct parametros				parametrosGA;
	struct parametros_Metodo	parametrosMetodo;

	inicializa_Parametros(&parmsPrograma, &parametrosGA, &parametrosMetodo);

	// Alocando memória para as gerações e suas repectivas
	// estruturas internas dinâmicas.
	struct generation				*geracao0, *geracao1;
	
	geracao0 = (struct generation*)malloc(sizeof(struct generation));
	criaIndividuosNaGeracao(&parametrosGA, geracao0);
	
	geracao1 = (struct generation*)malloc(sizeof(struct generation));
	criaIndividuosNaGeracao(&parametrosGA, geracao1);

	// Definindo as outras variáveis
	//struct testeGeracao_s		host_testeGeracao_s;
	unsigned long int				iGeracao;
	float								*Hamiltoniano;

	inicializa_Semente();

	// Imprime o cabecalho da marcação de tempo
	//imprimeTempo(0, 0, 0, 0, 0, 0, 0, 0);

	// GERA HAMILTONIANO
	time_i = clock();
	Hamiltoniano = (float *)malloc(parametrosGA.numGenes * parametrosGA.numGenes*sizeof(float));
	Gera_Matriz_de_Coope_LinhaColuna(Hamiltoniano, parametrosGA.numGenes);
	time_f = clock();
	//imprimeTempo(1, 0, 0, 0, 1, 2, time_i, time_f);

	// POPULAÇÃO INICIAL
	time_i = clock();
	GeraPopulacaoInicial_serial(geracao0, &parametrosGA); //, &host_testeGeracao_s);
	//GeraPopInicial_vetsOrtonormais(geracao0, &parametrosGA);
	//printf("\nPOPULACAO INICIAL:\n\n");
	//imprimeGeracao(geracao0,&parametrosGA);
	time_f = clock();
	//imprimeTempo(1, 0, 0, 0, 2, 2, time_i, time_f);

	// MATRIZ IDENTIDADE
	time_i = clock();
	float *MatrizIdentidade;
	MatrizIdentidade = (float *)malloc(parametrosGA.numGenes * parametrosGA.numGenes * sizeof(float));
	GeraMatrizIdentidade(MatrizIdentidade, parametrosGA.numGenes);
	time_f = clock();
	//imprimeTempo(1, 0, 0, 0, 3, 2, time_i, time_f);

	// Cabeçalho dos dados de comportamento do fitness
	imprimeComportamentoFitness(	0, 0, 0, 0, 0, geracao0,
											&parametrosMetodo, &parametrosGA, &parmsPrograma);

	iGeracao = 0;
	int flagAtingiuTolerancia = 0; // falso.
	float erroAbsolutoNoRho = 1.0F;
	float tolerancia = parmsPrograma.parmTolerancia;

	while (iGeracao < parmsPrograma.parmQtdeMaxGeracoes) {
		
		time_i = clock();
		Fitness_Serial(	geracao0,
								&parametrosGA,
								&parametrosMetodo,
								Hamiltoniano,
								parametrosGA.numGenes,
								parametrosMetodo.lambda,
								parametrosMetodo.rho_minimo,
								MatrizIdentidade);
		time_f = clock();
		//imprimeTempo(1, 0, 0, iGeracao, 4, 2, time_i, time_f);

		erroAbsolutoNoRho = geracao0->difRho;
		
		imprimeComportamentoFitness(	1, 0, 0, 0, iGeracao, geracao0,
												&parametrosMetodo, &parametrosGA, &parmsPrograma);

		// if (erroAbsolutoNoRho <= toleranciaErroRho) {
		if ( (geracao0->gradRhoMedio <= tolerancia) || (geracao0->difRho <= tolerancia) ) {
			flagAtingiuTolerancia = 1; // verdadeiro
			break; // sai do while
		}

		time_i = clock();
		//testeSelecao_serial(0, geracao0, &parametrosGA, &host_testeGeracao_s);
		Selecao_Por_Torneio_serial(geracao0, geracao1, &parametrosGA); //, &host_testeGeracao_s);
		//testeSelecao_serial(1, geracao1, &parametrosGA, &host_testeGeracao_s);
		time_f = clock();
		//imprimeTempo(1, 0, 0, iGeracao, 5, 2, time_i, time_f);

		//testeCrossOver_serial(0, geracao1, &parametrosGA, &host_testeGeracao_s);
		time_i = clock();
		//CrossOver1Ponto_serial(geracao1, geracao0, &parametrosGA); //, &host_testeGeracao_s);
		CrossOver2Pontos_serial(geracao1, geracao0, &parametrosGA); //, &host_testeGeracao_s);
		time_f = clock();
		//imprimeTempo(1, 0, 0, iGeracao, 6, 2, time_i, time_f);
		//testeCrossOver_serial(1, geracao0, &parametrosGA, &host_testeGeracao_s);

		//testeMutacao(0,geracao0,&parametrosGA, &host_testeGeracao_s);
		time_i = clock();
		Mutacao_Serial(geracao0, &parametrosGA); //, &host_testeGeracao_s);
		time_f = clock();
		//imprimeTempo(1, 0, 0, iGeracao, 7, 2, time_i, time_f);
		//testeMutacao(1,geracao0,&parametrosGA, &host_testeGeracao_s);

		iGeracao++;
	}
	
	if (flagAtingiuTolerancia == 0) {
		time_i = clock();
		Fitness_Serial(	geracao0,
								&parametrosGA,
								&parametrosMetodo,
								Hamiltoniano,
								parametrosGA.numGenes,
								parametrosMetodo.lambda,
								parametrosMetodo.rho_minimo,
								MatrizIdentidade);
		time_f = clock();
		//imprimeTempo(1, 0, 0, iGeracao, 4, 2, time_i, time_f);
		imprimeComportamentoFitness(	1, 0, 0, 0, iGeracao, geracao0,
												&parametrosMetodo, &parametrosGA, &parmsPrograma);
	}

	time_pgm_f = clock();
	imprimeTempo(	1, 0, 0, iGeracao, 0, 2, time_pgm_i, time_pgm_f,
						&parametrosGA, &parametrosMetodo, &parmsPrograma);

	printf("\n"); printf("Geracao final:");
	imprimeGeracao(geracao0, &parametrosGA);

	gravaEstatistica(geracao0, &parametrosGA, &parametrosMetodo);
	
	// ==================================================================
	
	free(MatrizIdentidade);

	return 0;
}

	/*
	float x = 5.8F;
	
	printf("\n"); printf("ceil(x) = %f", ceil(x));
	printf("\n"); printf("floor(x) = %f", floor(x));
	printf("\n"); printf("x_int = %d", (unsigned short int)x);
	*/

	//testa_Multiplica_Matrizes();
	//testa_multiplica_matriz_por_escalar();
	//testa_Subtrai_Matrizes();
	//testa_Calcula_Quociente_Rayleigh();
	//testa_Gradiente_de_Rho();
	//testa_imprimeMatriz();
	//testa_gera_Matriz_de_Coope_abs();
	//testa_Calcula_Quociente_Rayleigh();
	//testa_Gradiente_de_Rho();
	//testa_Fitness_Serial();
	//testa_Randomico();
	//testa_inicializa_Parametros();
	//teste_GeraPopulacaoInicial_serial();
	//testa_Gera_Matriz_de_Coope_LinhaColuna();