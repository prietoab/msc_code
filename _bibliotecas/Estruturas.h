#include <math.h>

struct parametrosPrograma {
	unsigned short int parmMaquina;
	unsigned short int parmSerial_ou_Paralelo;
	unsigned short int parmQtdeGenes;
	unsigned long int parmQtdeMaxGeracoes;
	unsigned short int parmQtdeIndividuos; // número de indivíduos por geração
	char parmTipoFitnessEquacao;
	char parmTipoFitnessParalelo;
	char parmTipoCalculoGradRho;
	char parmTamanhoTorneio;
	float	parmProbCrossOver;
	unsigned short int parmQtdePontosCorte;
	float parmProbMutacao;
	float parmIntensidadeMutacao;
	float parmLambda;
	float parmRhoMinimo;
	float parmTolerancia;
	unsigned short int parmFlagImprimeComportamentoFitness;
	unsigned short int parmFlagImprimeTempo;
	unsigned short int parmFlagNovaSemente;
	unsigned long int parmSemente;
};
struct individual {
	float							cociente_Rayleigh;
	float							grad_elevado_ao_quadrado;
	float							rho_menos_rho0_ao_quadrado; // (rho - rho_zero)^2, aquele utilizado no fitness.
	float							Numerador;
	float							CtC;
	float							inverso_de_CtC;
	float							fitness;
	unsigned short int		pontos_de_corte[2]; // fixo
	float							*gene; // alocar constNumGenes
	float							*gradRho; // alocar constNumGenes
};

struct generation {
	float							sumFitness;
	float							FitnessMedio;
	float							rhoMedio;
	float							sumRho;
	float							difRho;
	float							sumGradRho;
	float							gradRhoMedio;
	float							Maior_fitness;
	float							Melhor_cociente_Rayleigh;
	unsigned short int		idxMelhorIndividuo;
	struct individual			*individuo; // alocar constNumIndividuos;
};

struct parametros {
	unsigned short int		numGenes;
	unsigned short int		numIndividuos;
	float							probCrossOver;
	float							probMutacao;
	float							intensidadeMutacao;
	char							tamanho_torneio;
	unsigned short int		qtde_Pontos_de_Corte;
};

struct parametros_Metodo {
	float							lambda;
	float							rho_minimo;
};

//=====================================================================================
float getLambda(unsigned short int N) {
//=====================================================================================
	// log = logaritmo natural, neperiano
	float	lambda = -0.65*(log(0.000001)/pow((1.0001*N - 1.6507),2));

	return lambda;
}

//---------------------------------------------------------------------------------------
void inicializa_Parametros(struct parametrosPrograma *parmsPgm,
									struct parametros *parametrosGA,
									struct parametros_Metodo *parametrosMetodo ) {
//---------------------------------------------------------------------------------------

	parametrosGA->numGenes = parmsPgm->parmQtdeGenes;
	parametrosGA->numIndividuos = parmsPgm->parmQtdeIndividuos;
	parametrosGA->probCrossOver = parmsPgm->parmProbCrossOver;
	parametrosGA->probMutacao = parmsPgm->parmProbMutacao;
	parametrosGA->intensidadeMutacao = parmsPgm->parmIntensidadeMutacao;
	parametrosGA->qtde_Pontos_de_Corte = parmsPgm->parmQtdePontosCorte;
	parametrosGA->tamanho_torneio = parmsPgm->parmTamanhoTorneio;
	//parametrosGA->vlrMaximoGene = 0.1F;

	if (parmsPgm->parmLambda == -1)
		parametrosMetodo->lambda = getLambda(parametrosGA->numGenes);
	else
		parametrosMetodo->lambda = parmsPgm->parmLambda;
	
	parametrosMetodo->rho_minimo = parmsPgm->parmRhoMinimo;


	printf("\n");
	printf("+=====================================================+\n");
	printf("| PARAMETROS DE CHAMADA                                 \n");
	printf("+=====================================================+\n");
	printf("| Maquina : %d \n"														, parmsPgm->parmMaquina);
	printf("| Programa Serial (0) ou Paralelo (1) ? %d \n"				, parmsPgm->parmSerial_ou_Paralelo);
	printf("| Numero de genes (ordem do Hamiltoniano): %d \n"			, parmsPgm->parmQtdeGenes);
	printf("| Numero maximo de geracoes: %d    \n"							, parmsPgm->parmQtdeMaxGeracoes);
	printf("| Numero de individuos na populacao: %d \n"					, parmsPgm->parmQtdeIndividuos);
	printf("| Tipo da equacao do Fitness: %d \n"								, parmsPgm->parmTipoFitnessEquacao);
	printf("|    0: Apenas com (rho - rho_0)^2 \n");
	printf("|    1: Apenas com o (gradiente de rho)^2 \n");
	printf("|    2: Com (rho - rho_0)^2 + (grad(rho))^2 \n");
	printf("|    3: Apenas com grad(rho) \n");
	printf("|    4: Com (rho - rho_0)^2 + grad(rho) \n");
	printf("| Tipo da funcao Fitness no programa paralelo: %d \n"		, parmsPgm->parmTipoFitnessParalelo);
	printf("|    0: 1d, dim3(X,1,1) \n");
	printf("|    1: 2d, dim3(X,Y,1) \n");
	printf("| Tipo do Calculo do Gradiente de rho: %d \n" , parmsPgm->parmTipoCalculoGradRho);
	printf("|    0: usa matriz identidade \n");
	printf("|    1: nao usa matriz identidade (consome menos memoria) \n");
	printf("| Tamanho do torneio para a Selecao: %d \n"					, parmsPgm->parmTamanhoTorneio);
	printf("| Probabilidade de Crossover: %f \n"								, parmsPgm->parmProbCrossOver);
	printf("| Quantidade de pontos de corte para o CrossOver: %d \n"	, parmsPgm->parmQtdePontosCorte);
	printf("|    Observacao: atualmente esse quantidade foi fixada \n");
	printf("|                em 2 e, portanto, esse parametro nao \n");
	printf("|                na execucao do programa. \n");
	printf("| Probabilidade de mutacao: %f \n"								, parmsPgm->parmProbMutacao);
	printf("| Intensidade da mutacao: %f \n"									, parmsPgm->parmIntensidadeMutacao);
	printf("| Lambda do Fitness: %f\n"											, parmsPgm->parmLambda);
	printf("| Rho minimo do Fitness: %f \n"									, parmsPgm->parmRhoMinimo);
	printf("| Tolerancia para parada: %f \n"									, parmsPgm->parmTolerancia);
	printf("| Imprime comportamento do Fitness: %d \n"						, parmsPgm->parmFlagImprimeComportamentoFitness);
	printf("|    0: Nao \n");
	printf("|    1: Sim \n");
	printf("| Imprime marcacoes de tempo: %d \n"								, parmsPgm->parmFlagImprimeTempo);
	printf("|    0: Nao \n");
	printf("|    1: Sim \n");
	printf("| Gera nova semente de numeros pseudo-aleatorios: %d \n"	, parmsPgm->parmFlagNovaSemente);
	printf("|    0: Nao. Usa semente definida no parametro 'Semente'. \n");
	printf("|    1: Sim. Cria nova semente. Parametro 'Semente' ignorado. \n");
	printf("| Semente para geraco de numeros pseudo-aleatorios: %d \n", parmsPgm->parmSemente);
	printf("+=====================================================+\n");
	printf("| PARAMETROS DO ALGORITMO GENETICO \n");
	printf("+=====================================================+\n");
	printf("| Intensidade da mutacao: %f\n"									, parametrosGA->intensidadeMutacao);
	printf("| Numero de genes: %d\n"												, parametrosGA->numGenes);
	printf("| Quantidade de individuos na populacao: %d\n"				, parametrosGA->numIndividuos);
	printf("| Probabilidade de CrossOver: %f\n"								, parametrosGA->probCrossOver);
	printf("| Probabilidade de Mutacao: %f\n"									, parametrosGA->probMutacao);
	printf("| Quantidade de pontos de corte para o CrossOver: %d \n"	, parametrosGA->qtde_Pontos_de_Corte);
	printf("| Tamanho do Torneio para Selecao: %d\n"						, parametrosGA->tamanho_torneio);
	printf("+=====================================================+\n");
	printf("| PARAMETROS DO METODO \n");
	printf("+=====================================================+\n");
	printf("| Lambda do Fitness: %f\n"	, parametrosMetodo->lambda);
	printf("| Rho minimo: %f\n"			, parametrosMetodo->rho_minimo);
	printf("+=====================================================+\n");
	printf("| OUTROS PARAMETROS                                    \n");
	printf("+=====================================================+\n");
	printf("| RAND_MAX ........................: %d                \n", RAND_MAX);
	printf("| CLOCKS_PER_SEC...................: %d                \n", CLOCKS_PER_SEC);
	printf("+=====================================================+\n\n");

}


//=======================================================================================
void ImprimeGenesDoIndividuo(	struct individual *Individuo,
										struct parametros *prt_host_ParametrosGA,
										char flagTabulacao) {
//=======================================================================================
	unsigned int iGene;

	// flagTabulacao
	//		0: insere o caracter de tabulacao ANTES do valor do gene
	//		1: insere o caracter de tabulacao DEPOIS do valor do gene
	//		9: não insere caracter de tabulação

	for (iGene = 0; iGene < prt_host_ParametrosGA->numGenes; iGene++) {
		if (flagTabulacao == 0) printf("\t");
		printf("%1.6f", Individuo->gene[iGene]);
		if (flagTabulacao == 1) printf("\t");
	}
}

//===========================================================================================
void imprimeIndividuo(	unsigned short int idx_Individuo,
								struct individual *individuo,
								struct parametros *ptr_host_parametrosGA,
								char flagCabecalho) {
//===========================================================================================

	unsigned int iGene;
	
	if (flagCabecalho == 0) {
		// Imprime cabeçalho
		//	rho | Grad(rho) | fitness | g1 g2 g3 (...) gn | Grad(rho)^2 | rho^2 | Numerador | CtC | CtC^(-1) 
		// ----------------------------------------------------------------------------------------------------------------------------------------------------
		printf("iIndividuo");
		printf("\tRho");
		printf("\tGrad(rho)");
		printf("\tFitness");
		for (iGene = 0; iGene < ptr_host_parametrosGA->numGenes; iGene++) printf("\tg%d", iGene);
		printf("\tGrad(rho)^2");
		printf("\t(rho - rho_0)^2");
		printf("\tNumerador");
		printf("\tCtC");
		printf("\tCtC^(-1)");
		printf("\n");
	}

	// imprime detalhes
	printf("%d", idx_Individuo);
	printf("\t%1.6f", individuo->cociente_Rayleigh);
	printf("\t%1.6f", -1.0F);
	printf("\t%1.6f", individuo->fitness);
	ImprimeGenesDoIndividuo(individuo, ptr_host_parametrosGA, 0);
	printf("\t%1.6f", individuo->grad_elevado_ao_quadrado);
	printf("\t%1.6f", individuo->rho_menos_rho0_ao_quadrado);
	printf("\t%1.6f", individuo->Numerador);
	printf("\t%1.6f", individuo->CtC);
	printf("\t%1.6f", individuo->inverso_de_CtC);
}


//===========================================================================================
void imprimeGeracao(	struct generation *host_Geracao,
							struct parametros *host_ParametrosGA) {
//===========================================================================================

	unsigned int iIndividuo;
	unsigned int flagCabecalho = 0;
	for (iIndividuo = 0; iIndividuo < host_ParametrosGA->numIndividuos; iIndividuo++) {
		if (iIndividuo > 0) {
			flagCabecalho = 1;
		}
		printf("\n"); imprimeIndividuo(iIndividuo, &host_Geracao->individuo[iIndividuo], host_ParametrosGA, flagCabecalho);
	}
}