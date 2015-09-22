
struct parametrosPrograma {
	unsigned short int parmQtdeGenes;
	unsigned long int parmQtdeMaxGeracoes;
	unsigned short int parmQtdeIndividuos; // número de indivíduos por geração
	char parmTamanhoTorneio;
	float	parmProbCrossOver;
	unsigned short int parmQtdePontosCorte;
	float parmProbMutacao;
	float parmIntensidadeMutacao;
	float parmLambda;
	float parmRhoMinimo;
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
//	float							vlrMaximoGene;
};

struct parametros_Metodo {
//	unsigned int iIndividuo;
	float							lambda;
	float							rho_minimo;
};

// - testes seleçao e crossover - pode ser removido
struct individuo_s {
	unsigned int				*teste_iIndividuo_para_Torneio;
	int							iIndividuo_Vencedor;
	unsigned int				idx_Pai;
	unsigned int				idx_Mae;
	float							p;
	float							f;
	unsigned int				*pontos_de_corte;
	unsigned int				i_Gene_Global; // alocar constNumGenes
	unsigned int				Primeiro_Gene_Individuo_Global;
	unsigned int				Ultimo_Gene_Individuo_Global;
	unsigned int				Ponto_de_Corte_Individuo_Global;
		// para o teste da mutação
	unsigned short int		*gene_L; // alocar constNumGenes;
	int							*gene_termo_do_L; // alocar constNumGenes;
	float							*gene_r; // alocar constNumGenes;
	float							*gene_pAux; // alocar constNumGenes
};

struct testeGeracao_s {
	struct individuo_s		*individuo; // alocar constNumIndividuos
};
// - teste seleçao e crossover - pode ser removido

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

	parametrosMetodo->lambda = parmsPgm->parmLambda;
	parametrosMetodo->rho_minimo = parmsPgm->parmRhoMinimo;

	printf("\n");
	printf("+=====================================================+\n");
	printf("| PARAMETROS DO METODO                                 \n");
	printf("+=====================================================+\n");
	printf("| Fitness                                              \n");
	printf("|    Lambda .......................: %f                \n", parametrosMetodo->lambda);
	printf("|    Rho minimo ...................: %f                \n", parametrosMetodo->rho_minimo);
	printf("+=====================================================+\n");
	printf("| PARAMETROS DO ALGORITMO GENETICO                     \n");
	printf("+=====================================================+\n");
	printf("| Tamanho do cromossomo ...........: %04d              \n", parametrosGA->numGenes);
	printf("| Numero de individuos por geracao : %04d              \n", parametrosGA->numIndividuos);
	printf("| Numero maximo de geracoes .......: %04d              \n", parmsPgm->parmQtdeMaxGeracoes);
	printf("| Seleção? ........................: Sim               \n");
	printf("|    Metodo de selecao ............: Torneio           \n");
	printf("|    Tamanho do Torneio ...........: %d                \n", parametrosGA->tamanho_torneio);
	printf("| Crossover? ......................: Sim               \n");
	printf("|    Probabilidade de Crossover ...: %f                \n", parametrosGA->probCrossOver);
	printf("|    Quantidade de Pontos de Corte.: %d                \n", parametrosGA->qtde_Pontos_de_Corte);
	printf("| Mutacao? ........................: Sim               \n");
	printf("|    Probabilidade de Mutacao .....: %f                \n", parametrosGA->probMutacao);
	printf("|    Intensidade   da Mutacao .....: %f                \n", parametrosGA->intensidadeMutacao);
	printf("| Criterio de parada ..............: Numero maximo     \n");
	printf("|                                    de geracoes       \n");
	printf("+=====================================================+\n");
	printf("| OUTROS PARAMETROS                                    \n");
	printf("+=====================================================+\n");
	printf("| RAND_MAX ........................: %d                \n", RAND_MAX);
	printf("| CLOCKS_PER_SEC...................: %d                \n", CLOCKS_PER_SEC);
	printf("+=====================================================+\n\n");

}

//---------------------------------------------------------------------------------------
void testa_inicializa_Parametros(void) {
//---------------------------------------------------------------------------------------
	//TODO: alterar a inicialização de parâmetros

	/*

	struct parametros parametros_GA;
	struct parametros_Metodo parametros_Metodo;

	unsigned long int globalNumGeracoes = 1;

	inicializa_Parametros(&parametros_GA, &parametros_Metodo);

	printf("\n"); printf("Imprimindo fora da inicializa...");

	printf("\n");
	printf("+=====================================================+\n");
	printf("| PARAMETROS DO METODO                                 \n");
	printf("+=====================================================+\n");
	printf("| Fitness                                              \n");
	printf("|    Lambda .......................: %f                \n", parametros_Metodo.lambda);
	printf("|    Rho minimo ...................: %f                \n", parametros_Metodo.rho_minimo);
	printf("+=====================================================+\n");
	printf("| PARAMETROS DO ALGORITMO GENETICO                     \n");
	printf("+=====================================================+\n");
	printf("| Tamanho do cromossomo ...........: %04d              \n", parametros_GA.numGenes);
	printf("| Numero de individuos por geracao : %04d              \n", parametros_GA.numIndividuos);
	printf("| Numero maximo de geracoes .......: %04d              \n", globalNumGeracoes);
	printf("| Máximo gene (populaão inicial) ..: %f              \n", parametros_GA.vlrMaximoGene);
	printf("| Seleção? ........................: Sim               \n");
	printf("|    Metodo de selecao ............: Torneio           \n");
	printf("|    Tamanho do Torneio ...........: %d                \n", parametros_GA.tamanho_torneio);
	printf("| Crossover? ......................: Sim               \n");
	printf("|    Probabilidade de Crossover ...: %f                \n", parametros_GA.probCrossOver);
	printf("|    Quantidade de Pontos de Corte.: %d                \n", parametros_GA.qtde_Pontos_de_Corte);
	printf("| Mutacao? ........................: Sim               \n");
	printf("|    Probabilidade de Mutacao .....: %f                \n", parametros_GA.probMutacao);
	printf("|    Intensidade   da Mutacao .....: %f                \n", parametros_GA.intensidadeMutacao);
	printf("| Criterio de parada ..............: Numero maximo     \n");
	printf("|                                    de geracoes       \n");
	printf("+=====================================================+\n");
	printf("| OUTROS PARAMETROS                                    \n");
	printf("+=====================================================+\n");
	printf("| RAND_MAX ........................: %d                \n", RAND_MAX);
	printf("| CLOCKS_PER_SEC...................: %d                \n", CLOCKS_PER_SEC);
	printf("+=====================================================+\n\n");

	*/
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

	/*
	printf("\n");
	printf("===> Ponteiro do *Individuo na ImprimeGenesDoIndividuo:\t%p", Individuo);
	printf("\n");
	printf("===> Ponteiro do *prt_host_ParametrosGA na ImprimeGenesDoIndividuo:\t%p", prt_host_ParametrosGA);
	printf("\n");
	*/
	// Endereços dos genes pra teste
	/*
	printf("\n");
	for (iGene = 0; iGene < prt_host_ParametrosGA->numGenes; iGene++) {
		if (flagTabulacao == 0) printf("\t");
		printf("%p", &Individuo->gene[iGene]);
		if (flagTabulacao == 1) printf("\t");
	};
	printf("\n");
	*/

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

//	printf("\n");	printf("*====================================================");
//	printf("\n");  printf(" IMPRIME GERACAO - INICIO");
//	printf("\n");  printf("*====================================================");
//	printf("\n");  printf("ATRIBUTOS DA GERACAO");
//	printf("\n");  printf("------------------------------");
//	printf("\n");  printf("\t\tSoma do fitnes dos individuos: \t%1.6f", host_Geracao->sumFitness);
//	printf("\n");  printf("\t\tFitness medio da geracao ....: \t%1.6f", host_Geracao->FitnessMedio);
//	printf("\n");  printf("\t\tMelhor cociente de Rayleigh .: \t%1.6f", host_Geracao->melhores.Melhor_cociente_Rayleigh);
//	printf("\n");  printf("\t\tPosicao do melhor individuo .: \t%d", host_Geracao->melhores.idxMelhorIndividuo);
//	printf("\n");   
//	printf("\n");  printf("MELHOR INDIVIDUO DA GERACAO");
//	printf("\n");  printf("------------------------------");
//	printf("\n");  imprimeIndividuo(&host_Geracao->melhores.MelhorIndividuo, host_ParametrosGA, 0);
//	printf("\n");   
//	printf("\n");  printf("TODOS OS INDIVIDUOS DA GERACAO");
//	printf("\n");  printf("------------------------------");
	
	unsigned int iIndividuo;
	unsigned int flagCabecalho = 0;
	for (iIndividuo = 0; iIndividuo < host_ParametrosGA->numIndividuos; iIndividuo++) {
		if (iIndividuo > 0) {
			flagCabecalho = 1;
		}
		printf("\n"); imprimeIndividuo(iIndividuo, &host_Geracao->individuo[iIndividuo], host_ParametrosGA, flagCabecalho);
	}

//	printf("\n");	printf("+====================================================");
//	printf("\n");  printf("+IMPRIME GERACAO - FIM");
//	printf("\n");  printf("+====================================================");
}

