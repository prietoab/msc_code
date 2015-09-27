
//---------------------------------------------------------------------------------------
unsigned long int inicializa_Semente(unsigned short int flagNovaSemente, unsigned long int sementeExterna) {
//---------------------------------------------------------------------------------------
   unsigned long int semente;
	
	if (flagNovaSemente == 1) // verdadeiro
		semente = (unsigned long int)time (NULL);
	else
		semente = sementeExterna;

	printf("\nSemente = %d", semente);	
	srand(semente);
	return semente;
}

//---------------------------------------------------------------------------------------
int Randomico_int(unsigned short int range_min, unsigned short int range_max) {
//---------------------------------------------------------------------------------------
   // Generate random numbers in the half-closed interval
   // [range_min, range_max). In other words,
   // range_min <= random number < range_max
	/*if ((range_min > RAND_MAX) || (range_min < 0) || (range_max < 0) || (range_max > RAND_MAX)) {
		printf("===> Erro! Valor do parametro fora do intervalor [0,RAND_MAX]");
		return 9;
	};*/
	
	unsigned short int u;

	if (range_min > RAND_MAX) {
		printf("===> Erro! range_min > RAND_MAX\n");
		return 9;
	}
	else {
		if (range_min < 0) {
			printf("===> Erro! range_min < 0\n");
			return 10;
		}
		else {
			if (range_max < 0) {
				printf("===> Erro! range_max < 0\n");
				return 11;
			}
			else {
				if (range_max > RAND_MAX) {
					printf("===> Erro! range_max > RAND_MAX\n");
					return 12;
				}
				else {
					//u = (float)( (rand() / (RAND_MAX + 1)) * (range_max - range_min) + range_min);
					//u = (float)rand() / (RAND_MAX + 1) * (range_max - range_min) + range_min;
					u = (float)rand() / (RAND_MAX + 1) * (range_max - range_min) + range_min;
					return u;
				}
			}
		}
	};
};

//---------------------------------------------------------------------------------------
void testa_Randomico() {
//---------------------------------------------------------------------------------------

	inicializa_Semente(1, 0);
	
	unsigned short int minimo, maximo;
	unsigned long int iNumero, qtdeNumeros;

	qtdeNumeros = 1000000;
	minimo = 0;
	maximo = 100;
	
	printf("\n"); printf("iNumero\tNumero");
	for (iNumero = 0; iNumero < qtdeNumeros ; iNumero++) {
		printf("\n"); printf("%d\t%d", iNumero, Randomico_int(minimo,maximo));
	}
}

//---------------------------------------------------------------------------------------
unsigned short int imprimeTempo(
							unsigned short int flag_Cabecalho,
							unsigned short int cod_Maquina,
							unsigned short int cod_tipoPrograma,
							unsigned long	int idx_Geracao,
							unsigned short int cod_Funcao,
							unsigned short int cod_TipoMarcacao,
							clock_t clock_Inicio,
							clock_t clock_Fim,
							struct parametros *parametrosGA,
							struct parametros_Metodo *parametrosMetodo,
							struct parametrosPrograma *parmsPrograma) {
//---------------------------------------------------------------------------------------

	if (flag_Cabecalho == 0) {
		printf("\n"); printf("Tempo de processamento");
		printf("\t"); printf("Maquina");						// 0 - Adriano, 1 - FT
		printf("\t"); printf("Serial ou Paralelo?");		// 0 - Serial, 1 - Paralelo
		printf("\t"); printf("Quantidade Geracoes");		// constNumGeracoes
		printf("\t"); printf("Quantidade Individuos");	// constNumIndividuos
		printf("\t"); printf("Quantidade Genes");			// constNumGenes
		printf("\t"); printf("idx_Geracao");
		printf("\t"); printf("Funcao");
		printf("\t"); printf("Tipo Marcacao");				// 0 - Inicio, 1 - Fim, 2 - Duracao
		printf("\t"); printf("Marcacao");					/* Depende do tipo da marcação:
																			0: clock_Inicio
																			1: clock_Fim
																			2: clock_Fim - clock_nicio
																		*/
	}
	else {
		
		char *str_TipoPrograma;
		switch (cod_tipoPrograma) {
			case 0: str_TipoPrograma = "00 - Serial\0";	break;
			case 1: str_TipoPrograma = "01 - Paralelo\0"; break;
			default:	{
					printf("\n");
					printf("Erro\timprimeTempo\tCodigo invalido para o tipo de programa.");
					return 2;
			} break;
		}

		char *str_Funcao;
		switch (cod_Funcao) {
			case  0:	str_Funcao = "00 - Programa\0"; break;
			case  1:	str_Funcao = "01 - Gera Hamiltoniano\0";	break;
			case  2:	str_Funcao = "02 - Gera Populacao Inicial\0"; break;
			case  3: str_Funcao = "03 - Gera Matriz Identidade\0"; break;
			case  4: str_Funcao = "04 - Fitness 1d\0"; break;
			case  5: str_Funcao = "05 - Selecao\0"; break;
			case  6: str_Funcao = "06 - Crossover\0"; break;
			case  7: str_Funcao = "07 - Mutacao\0"; break;
			case  8: str_Funcao = "04 - Fitness 2d\0"; break;
			default: {
				printf("\n");
				printf("Erro\timprimeTempo\tCodigo invalido para a funcao.");
				return 3;
			}; break;
		}


		char *str_tipoMarcacao;
		clock_t valor_Marcacao;

		switch (cod_TipoMarcacao) {
			case 0: {
				str_tipoMarcacao = "00 - Inicio\0";
				valor_Marcacao = clock_Inicio;
			}
			break;
			case 1: {
				str_tipoMarcacao = "01 - Fim\0"; 
				valor_Marcacao = clock_Fim;
			}
			break;
			case 2: {
				str_tipoMarcacao = "02 - Duracao\0"; 
				valor_Marcacao = clock_Fim - clock_Inicio;
			}
			break;
			default: {
				printf("\n");
				printf("Erro\timprimeTempo\tCodigo invalido para o tipo de marcacao.");
				return 4;
			}
			break;

		}
		
		printf("\n"); printf("Tempo de processamento");
		printf("\t"); printf("%d", cod_Maquina);			
		printf("\t"); printf("%s", str_TipoPrograma);	
		printf("\t"); printf("%d", parmsPrograma->parmQtdeMaxGeracoes);	
		printf("\t"); printf("%d", parametrosGA->numIndividuos);	
		printf("\t"); printf("%d", parametrosGA->numGenes);		
		printf("\t"); printf("%d", idx_Geracao);		 															
		printf("\t"); printf("%s", str_Funcao);		 															
		printf("\t"); printf("%s", str_tipoMarcacao);
		printf("\t"); printf("%d", valor_Marcacao);
	} // final do if flag_Cabealho;

	return 0;
}

//---------------------------------------------------------------------------------------
void testa_imprimeTempo() {
//---------------------------------------------------------------------------------------
	//TODO: corrigir parâmetros da função imprimeTempo	
	unsigned short int	int_Cabecalho = 999,
								iMaquina = 9999,
								iTipoPrograma = 999,
								iFuncao = 999,
								iMarcacao = 999;
								
	clock_t	tempo_inicio = 999,
				tempo_fim = 999;

	tempo_inicio = 3;
	tempo_fim = 70;

	int_Cabecalho = 0;
	//imprimeTempo(int_Cabecalho, iMaquina, iTipoPrograma, 0, iFuncao, iMarcacao, clock(), clock());

	int_Cabecalho = 1;
	for (iMaquina = 0; iMaquina < 2; iMaquina++) {
		for (iTipoPrograma = 0; iTipoPrograma < 2; iTipoPrograma++) {
			for (iFuncao = 0; iFuncao < 8; iFuncao++) {
				for (iMarcacao = 0; iMarcacao < 3; iMarcacao++) {
					//imprimeTempo(int_Cabecalho, iMaquina, iTipoPrograma, 0, iFuncao, iMarcacao, tempo_inicio, tempo_fim);
				}
			}
		}
	}
}
//---------------------------------------------------------------------------------------
void criaIndividuosNaGeracao(		struct parametros *parametrosGA,
											struct generation *geracao) {
//---------------------------------------------------------------------------------------
	unsigned short int j, qtdeInd, qtdeGenes;

	qtdeInd = parametrosGA->numIndividuos;
	qtdeGenes = parametrosGA->numGenes;

	geracao->individuo = (struct individual *)malloc(qtdeInd*sizeof(struct individual));
	for (j = 0; j < qtdeInd; j++) {
		geracao->individuo[j].gene = (float *)malloc(qtdeGenes*sizeof(float));
		geracao->individuo[j].gradRho = (float *)malloc(qtdeGenes*sizeof(float));
	}
}