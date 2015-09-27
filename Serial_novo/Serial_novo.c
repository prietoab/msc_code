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
	
	parmsPrograma.parmMaquina									= atoi(argv[1]);
	parmsPrograma.parmSerial_ou_Paralelo					= atoi(argv[2]);
	parmsPrograma.parmQtdeGenes								= atoi(argv[3]);
	parmsPrograma.parmQtdeMaxGeracoes						= atoi(argv[4]);
	parmsPrograma.parmQtdeIndividuos							= atoi(argv[5]); // número de indivíduos por geração
	parmsPrograma.parmTipoFitness								= atoi(argv[6]);
	parmsPrograma.parmTamanhoTorneio							= atoi(argv[7]);
	parmsPrograma.parmProbCrossOver							= (float)atof(argv[8])/100;
	parmsPrograma.parmQtdePontosCorte						= atoi(argv[9]);
	parmsPrograma.parmProbMutacao								= (float)atof(argv[10])/100;
	parmsPrograma.parmIntensidadeMutacao					= (float)atof(argv[11])/10;
	parmsPrograma.parmLambda									= (float)atof(argv[12]);
	parmsPrograma.parmRhoMinimo								= (float)atof(argv[13]);
	parmsPrograma.parmTolerancia								= (float)atof(argv[14]);
	parmsPrograma.parmFlagImprimeComportamentoFitness	= atoi(argv[15]);
	parmsPrograma.parmFlagImprimeTempo						= atoi(argv[16]);
	parmsPrograma.parmFlagNovaSemente						= atoi(argv[17]);
	parmsPrograma.parmSemente									= atoi(argv[18]);
	
	/*
	parmsPrograma.parmMaquina									= 0;
	parmsPrograma.parmSerial_ou_Paralelo					= 0;
	parmsPrograma.parmQtdeGenes								= 100;
	parmsPrograma.parmQtdeMaxGeracoes						= 500000;
	parmsPrograma.parmQtdeIndividuos							= 200; // número de indivíduos por geração
	parmsPrograma.parmTipoFitness								= 1;
	parmsPrograma.parmTamanhoTorneio							= 2;
	parmsPrograma.parmProbCrossOver							= 90;
	parmsPrograma.parmQtdePontosCorte						= 2;
	parmsPrograma.parmProbMutacao								= (float)30/100;
	parmsPrograma.parmIntensidadeMutacao					= (float)40/10;
	parmsPrograma.parmLambda									= (float)-1;
	parmsPrograma.parmRhoMinimo								= (float)0.2;
	parmsPrograma.parmTolerancia								= (float)0.000001;
	parmsPrograma.parmFlagImprimeComportamentoFitness	= 0;
	parmsPrograma.parmFlagImprimeTempo						= 0;
	parmsPrograma.parmFlagNovaSemente						= 1;
	parmsPrograma.parmSemente									= 1443205662;

	
	printf("parmsPrograma.parmMaquina = %d \n", parmsPrograma.parmMaquina);
	printf("parmsPrograma.parmSerial_ou_Paralelo = %d \n", parmsPrograma.parmSerial_ou_Paralelo);
	printf("parmsPrograma.parmQtdeGenes = %d \n", parmsPrograma.parmQtdeGenes);
	printf("parmsPrograma.parmQtdeMaxGeracoes = %d \n", parmsPrograma.parmQtdeMaxGeracoes);
	printf("parmsPrograma.parmQtdeIndividuos = %d \n", parmsPrograma.parmQtdeIndividuos);
	printf("parmsPrograma.parmTipoFitness = %d \n", parmsPrograma.parmTipoFitness);
	printf("parmsPrograma.parmTamanhoTorneio = %d \n", parmsPrograma.parmTamanhoTorneio);
	printf("parmsPrograma.parmProbCrossOver = %f \n", parmsPrograma.parmProbCrossOver);
	printf("parmsPrograma.parmQtdePontosCorte = %d \n", parmsPrograma.parmQtdePontosCorte);
	printf("parmsPrograma.parmProbMutacao = %f \n", parmsPrograma.parmProbMutacao);
	printf("parmsPrograma.parmIntensidadeMutacao = %f \n", parmsPrograma.parmIntensidadeMutacao);
	printf("parmsPrograma.parmLambda = %f \n", parmsPrograma.parmLambda);
	printf("parmsPrograma.parmRhoMinimo = %f \n", parmsPrograma.parmRhoMinimo);
	printf("parmsPrograma.parmTolerancia = %f \n", parmsPrograma.parmTolerancia);
	printf("parmsPrograma.parmFlagImprimeComportamentoFitness = %d\n", parmsPrograma.parmFlagImprimeComportamentoFitness);
	printf("parmsPrograma.parmFlagImprimeTempo = %d \n", parmsPrograma.parmFlagImprimeTempo);
	printf("parmFlagNovaSemente = %d\n", parmsPrograma.parmFlagNovaSemente);
	printf("parmSemente = %d\n", parmsPrograma.parmSemente);

	return 0;
	*/

	// Interface entre os parãmetros e as variáveis locais
	int flagAtingiuTolerancia = 0; // falso.
	float erroAbsolutoNoRho = 1.0F;
	float tolerancia = parmsPrograma.parmTolerancia;
	unsigned short int tipoFitness = parmsPrograma.parmTipoFitness;
	unsigned short int codMaquina = parmsPrograma.parmMaquina;
	unsigned short int tipoPrograma = parmsPrograma.parmSerial_ou_Paralelo;
	unsigned short int flagImprimeComportamentoFitness = parmsPrograma.parmFlagImprimeComportamentoFitness;
	unsigned short int flagImprimeTempo = parmsPrograma.parmFlagImprimeTempo;
	unsigned short int flagSemente = parmsPrograma.parmFlagNovaSemente;
	unsigned long int Semente = parmsPrograma.parmSemente;

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
	unsigned long int				iGeracao = 0;
	float								*Hamiltoniano;

	unsigned long int semente = inicializa_Semente(flagSemente, Semente);

	// Imprime o cabecalho da marcação de tempo
	if (flagImprimeTempo == 1) {
		imprimeTempo(0, 0, 0, 0, 0, 0, 0, 0, &parametrosGA, &parametrosMetodo, &parmsPrograma);
	}

	// GERA HAMILTONIANO
	time_i = clock();
	Hamiltoniano = (float *)malloc(parametrosGA.numGenes * parametrosGA.numGenes*sizeof(float));
	Gera_Matriz_de_Coope_LinhaColuna(Hamiltoniano, parametrosGA.numGenes);
	time_f = clock();
	if (flagImprimeTempo == 1) {
		imprimeTempo(1, 0, 0, 0, 1, 2, time_i, time_f, &parametrosGA, &parametrosMetodo, &parmsPrograma);
	}

	// POPULAÇÃO INICIAL
	time_i = clock();
	GeraPopulacaoInicial_serial(geracao0, &parametrosGA); //, &host_testeGeracao_s);
	//GeraPopInicial_vetsOrtonormais(geracao0, &parametrosGA);
	//printf("\nPOPULACAO INICIAL:\n\n");
	//imprimeGeracao(geracao0,&parametrosGA);
	time_f = clock();
	if (flagImprimeTempo == 1) {
		imprimeTempo(1, 0, 0, 0, 2, 2, time_i, time_f, &parametrosGA, &parametrosMetodo, &parmsPrograma);
	}

	// MATRIZ IDENTIDADE
	time_i = clock();
	float *MatrizIdentidade;
	MatrizIdentidade = (float *)malloc(parametrosGA.numGenes * parametrosGA.numGenes * sizeof(float));
	GeraMatrizIdentidade(MatrizIdentidade, parametrosGA.numGenes);
	time_f = clock();
	if (flagImprimeTempo == 1) {
		imprimeTempo(1, 0, 0, 0, 3, 2, time_i, time_f, &parametrosGA, &parametrosMetodo, &parmsPrograma);
	}

	// "Em C, não existe o tipo booleano (true ou false). O valor zero é interpretado
	// como falso e qualquer valor diferente de zero é considerado verdadeiro."

	// Cabeçalho dos dados de comportamento do fitness
	if (flagImprimeComportamentoFitness == 1) {
		imprimeComportamentoFitness(	0, codMaquina, tipoPrograma, 0, 0, geracao0,
												&parametrosMetodo, &parametrosGA, &parmsPrograma);
	}

	while (iGeracao < parmsPrograma.parmQtdeMaxGeracoes) {
		
		time_i = clock();
		Fitness_Serial(	tipoFitness,
								geracao0,
								&parametrosGA,
								&parametrosMetodo,
								Hamiltoniano,
								parametrosGA.numGenes,
								parametrosMetodo.lambda,
								parametrosMetodo.rho_minimo,
								MatrizIdentidade);
		time_f = clock();
		
		if (flagImprimeTempo == 1) {
			imprimeTempo(1, 0, 0, iGeracao, 4, 2, time_i, time_f, &parametrosGA, &parametrosMetodo, &parmsPrograma);
		}

		if (flagImprimeComportamentoFitness == 1) {
			imprimeComportamentoFitness(	1, codMaquina, tipoPrograma, 0, iGeracao, geracao0,
													&parametrosMetodo, &parametrosGA, &parmsPrograma);
		}

		// Critérios de parada.
		flagAtingiuTolerancia = atingiuCriterioDeParada(tipoFitness, tolerancia, geracao0);
		if (flagAtingiuTolerancia == 1)
			break;

		time_i = clock();
		//testeSelecao_serial(0, geracao0, &parametrosGA, &host_testeGeracao_s);
		Selecao_Por_Torneio_serial(geracao0, geracao1, &parametrosGA); //, &host_testeGeracao_s);
		//testeSelecao_serial(1, geracao1, &parametrosGA, &host_testeGeracao_s);
		time_f = clock();
		if (flagImprimeTempo == 1) {
			imprimeTempo(1, 0, 0, iGeracao, 5, 2, time_i, time_f, &parametrosGA, &parametrosMetodo, &parmsPrograma);
		}
		//testeCrossOver_serial(0, geracao1, &parametrosGA, &host_testeGeracao_s);
		time_i = clock();
		CrossOver2Pontos_serial(geracao1, geracao0, &parametrosGA); //, &host_testeGeracao_s);
		time_f = clock();
		if (flagImprimeTempo == 1) {
			imprimeTempo(1, 0, 0, iGeracao, 6, 2, time_i, time_f, &parametrosGA, &parametrosMetodo, &parmsPrograma);
		}
		//testeCrossOver_serial(1, geracao0, &parametrosGA, &host_testeGeracao_s);

		//testeMutacao(0,geracao0,&parametrosGA, &host_testeGeracao_s);
		time_i = clock();
		Mutacao_Serial(geracao0, &parametrosGA); //, &host_testeGeracao_s);
		time_f = clock();
		if (flagImprimeTempo == 1) {
			imprimeTempo(1, 0, 0, iGeracao, 7, 2, time_i, time_f, &parametrosGA, &parametrosMetodo, &parmsPrograma);
		}
		//testeMutacao(1,geracao0,&parametrosGA, &host_testeGeracao_s);

		iGeracao++;
	}
	
	if (flagAtingiuTolerancia == 0) {
		time_i = clock();
		Fitness_Serial(	tipoFitness,
								geracao0,
								&parametrosGA,
								&parametrosMetodo,
								Hamiltoniano,
								parametrosGA.numGenes,
								parametrosMetodo.lambda,
								parametrosMetodo.rho_minimo,
								MatrizIdentidade);
		time_f = clock();
		if (flagImprimeTempo == 1) {
			imprimeTempo(1, 0, 0, iGeracao, 4, 2, time_i, time_f, &parametrosGA, &parametrosMetodo, &parmsPrograma);
		}

		if (flagImprimeComportamentoFitness == 1) {
			imprimeComportamentoFitness(	1, codMaquina, tipoPrograma, 0, iGeracao, geracao0,
													&parametrosMetodo, &parametrosGA, &parmsPrograma);
		}
	}

	time_pgm_f = clock();

	if (flagImprimeTempo == 1) {
		imprimeTempo(	1, 0, 0, iGeracao, 0, 2, time_pgm_i, time_pgm_f,
						&parametrosGA, &parametrosMetodo, &parmsPrograma);
	}

	unsigned long int tempoTotalProcessamento = time_pgm_f - time_pgm_i;
	gravaEstatistica(	tempoTotalProcessamento,
							semente,
							iGeracao,
							geracao0,
							&parmsPrograma,
							&parametrosGA,
							&parametrosMetodo);

	printf("\n"); printf("Geracao final:");
	imprimeGeracao(geracao0, &parametrosGA);
	
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