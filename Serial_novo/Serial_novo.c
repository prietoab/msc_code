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
int main(void) {
//==========================================================================

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
	struct generation				geracao[2];
	struct testeGeracao_s		host_testeGeracao_s;
	unsigned short int			iGeracao;
	float								*Hamiltoniano;

	inicializa_Parametros(&parametrosGA, &parametrosMetodo);
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
	GeraPopulacaoInicial_serial(&geracao[0], &parametrosGA, &host_testeGeracao_s);
	//GeraPopInicial_vetsOrtonormais(&geracao[0], &parametrosGA);
	printf("\nPOPULACAO INICIAL:\n\n");
	imprimeGeracao(&geracao[0],&parametrosGA);
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
	//imprimeComportamentoFitness(0, 0, 0, 0, 0, &geracao[0], &parametrosMetodo);

	for (iGeracao = 0; iGeracao < parametrosGA.numGeracoes; iGeracao++) {
		
		time_i = clock();
		Fitness_Serial(	&geracao[0],
								&parametrosGA,
								Hamiltoniano,
								parametrosGA.numGenes,
								parametrosMetodo.lambda,
								parametrosMetodo.rho_minimo,
								MatrizIdentidade);
		time_f = clock();
		//imprimeTempo(1, 0, 0, iGeracao, 4, 2, time_i, time_f);
		//imprimeComportamentoFitness(1, 0, 0, 0, iGeracao, &geracao[0], &parametrosMetodo);


		time_i = clock();
		//testeSelecao_serial(0, &geracao[0], &parametrosGA, &host_testeGeracao_s);
		Selecao_Por_Torneio_serial(&geracao[0], &geracao[1], &parametrosGA); //, &host_testeGeracao_s);
		//testeSelecao_serial(1, &geracao[1], &parametrosGA, &host_testeGeracao_s);
		time_f = clock();
		//imprimeTempo(1, 0, 0, iGeracao, 5, 2, time_i, time_f);

		testeCrossOver_serial(0, &geracao[1], &parametrosGA, &host_testeGeracao_s);
		time_i = clock();
		//CrossOver1Ponto_serial(&geracao[1], &geracao[0], &parametrosGA); //, &host_testeGeracao_s);
		CrossOver2Pontos_serial(&geracao[1], &geracao[0], &parametrosGA, &host_testeGeracao_s);
		time_f = clock();
		//imprimeTempo(1, 0, 0, iGeracao, 6, 2, time_i, time_f);
		testeCrossOver_serial(1, &geracao[0], &parametrosGA, &host_testeGeracao_s);

		//testeMutacao(0,&geracao[0],&parametrosGA, &host_testeGeracao_s);
		time_i = clock();
		Mutacao_Serial(&geracao[0], &parametrosGA); //, &host_testeGeracao_s);
		time_f = clock();
		//imprimeTempo(1, 0, 0, iGeracao, 7, 2, time_i, time_f);
		//testeMutacao(1,&geracao[0],&parametrosGA, &host_testeGeracao_s);
	}
	
	time_i = clock();
	Fitness_Serial(	&geracao[0],
							&parametrosGA,
							Hamiltoniano,
							parametrosGA.numGenes,
							parametrosMetodo.lambda,
							parametrosMetodo.rho_minimo,
							MatrizIdentidade);
	time_f = clock();
	//imprimeTempo(1, 0, 0, iGeracao, 4, 2, time_i, time_f);
	//imprimeComportamentoFitness(1, 0, 0, 0, iGeracao, &geracao[0], &parametrosMetodo);
	
	time_pgm_f = clock();
	imprimeTempo(1, 0, 0, iGeracao, 0, 2, time_pgm_i, time_pgm_f);

	printf("\n"); printf("Geracao final:");
	imprimeGeracao(&geracao[0], &parametrosGA);

	gravaEstatistica(&geracao[0], &parametrosGA, &parametrosMetodo);
	
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