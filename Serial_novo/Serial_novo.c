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
	parmsPrograma.parmTipoFitnessEquacao					= atoi(argv[6]);
	parmsPrograma.parmTipoFitnessParalelo					= atoi(argv[7]);
	parmsPrograma.parmTamanhoTorneio							= atoi(argv[8]);
	parmsPrograma.parmProbCrossOver							= (float)atof(argv[9])/100;
	parmsPrograma.parmQtdePontosCorte						= atoi(argv[10]);
	parmsPrograma.parmProbMutacao								= (float)atof(argv[11])/100;
	parmsPrograma.parmIntensidadeMutacao					= (float)atof(argv[12])/10;
	parmsPrograma.parmLambda									= (float)atof(argv[13]);
	parmsPrograma.parmRhoMinimo								= (float)atof(argv[14]);
	parmsPrograma.parmTolerancia								= (float)atof(argv[15]);
	parmsPrograma.parmFlagImprimeComportamentoFitness	= atoi(argv[16]);
	parmsPrograma.parmFlagImprimeTempo						= atoi(argv[17]);
	parmsPrograma.parmFlagNovaSemente						= atoi(argv[18]);
	parmsPrograma.parmSemente									= atoi(argv[19]);
	parmsPrograma.parmTipoCalculoGradRho					= atoi(argv[20]);

	// Interface entre os parãmetros e as variáveis locais
	int flagAtingiuTolerancia = 0; // falso.
	float erroAbsolutoNoRho = 1.0F;
	float tolerancia = parmsPrograma.parmTolerancia;
	unsigned short int parmTipoFitnessEquacao = parmsPrograma.parmTipoFitnessEquacao;
	unsigned short int codMaquina = parmsPrograma.parmMaquina;
	unsigned short int tipoPrograma = parmsPrograma.parmSerial_ou_Paralelo;
	unsigned short int flagImprimeComportamentoFitness = parmsPrograma.parmFlagImprimeComportamentoFitness;
	unsigned short int flagImprimeTempo = parmsPrograma.parmFlagImprimeTempo;
	unsigned short int flagSemente = parmsPrograma.parmFlagNovaSemente;
	unsigned long int Semente = parmsPrograma.parmSemente;
	char TipoCalculoGradRho = parmsPrograma.parmTipoCalculoGradRho;

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
	GeraPopulacaoInicial_serial(geracao0, &parametrosGA);
	time_f = clock();
	if (flagImprimeTempo == 1) {
		imprimeTempo(1, 0, 0, 0, 2, 2, time_i, time_f, &parametrosGA, &parametrosMetodo, &parmsPrograma);
	}

	// MATRIZ IDENTIDADE
	time_i = clock();
	float *MatrizIdentidade;
	if (TipoCalculoGradRho == 0) {
		MatrizIdentidade = (float *)malloc(parametrosGA.numGenes * parametrosGA.numGenes * sizeof(float));
		GeraMatrizIdentidade(MatrizIdentidade, parametrosGA.numGenes);
	}
	else {
		MatrizIdentidade = (float *)malloc(sizeof(float));
	}
	time_f = clock();
	if (flagImprimeTempo == 1) {
		imprimeTempo(1, 0, 0, 0, 3, 2, time_i, time_f, &parametrosGA, &parametrosMetodo, &parmsPrograma);
	}

	// Cabeçalho dos dados de comportamento do fitness
	if (flagImprimeComportamentoFitness == 1) {
		imprimeComportamentoFitness(	0, codMaquina, tipoPrograma, 0, 0, geracao0,
												&parametrosMetodo, &parametrosGA, &parmsPrograma);
	}

	while (iGeracao < parmsPrograma.parmQtdeMaxGeracoes) {
		
		// FITNESS
		time_i = clock();
		Fitness_Serial(	parmTipoFitnessEquacao,
								TipoCalculoGradRho,
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
		flagAtingiuTolerancia = atingiuCriterioDeParada(parmTipoFitnessEquacao, tolerancia, geracao0);
		if (flagAtingiuTolerancia == 1)
			break;

		// SELEÇÃO
		time_i = clock();
		Selecao_Por_Torneio_serial(geracao0, geracao1, &parametrosGA);
		time_f = clock();
		if (flagImprimeTempo == 1) {
			imprimeTempo(1, 0, 0, iGeracao, 5, 2, time_i, time_f, &parametrosGA, &parametrosMetodo, &parmsPrograma);
		}

		// CROSSOVER
		time_i = clock();
		CrossOver2Pontos_serial(geracao1, geracao0, &parametrosGA);
		time_f = clock();
		if (flagImprimeTempo == 1) {
			imprimeTempo(1, 0, 0, iGeracao, 6, 2, time_i, time_f, &parametrosGA, &parametrosMetodo, &parmsPrograma);
		}

		// MUTACAO
		time_i = clock();
		Mutacao_Serial(geracao0, &parametrosGA);
		time_f = clock();
		if (flagImprimeTempo == 1) {
			imprimeTempo(1, 0, 0, iGeracao, 7, 2, time_i, time_f, &parametrosGA, &parametrosMetodo, &parmsPrograma);
		}

		iGeracao++;
	}
	
	if (flagAtingiuTolerancia == 0) {
		time_i = clock();
		Fitness_Serial(	parmTipoFitnessEquacao,
								TipoCalculoGradRho,
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

	printf("\n\n"); printf("Geracao final:");
	imprimeGeracao(geracao0, &parametrosGA);
	
	// ==================================================================
	
	free(geracao0);
	free(geracao1);
	free(Hamiltoniano);
	free(MatrizIdentidade);

	return 0;
}