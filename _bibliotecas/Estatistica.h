
// -----------------------------------------------------------------------
float mediaRho(	struct generation *host_geracao,
						struct parametros *host_parametrosGA) {
// -----------------------------------------------------------------------
	float somaRho = 0;
	float mediaRho = 0;
	unsigned short int iIndividuos = 0;

	for (iIndividuos = 0; iIndividuos < host_parametrosGA->numIndividuos; iIndividuos++) {
		somaRho = somaRho + (float)host_geracao->individuo[iIndividuos].cociente_Rayleigh;
	}

	mediaRho = somaRho / host_parametrosGA->numIndividuos;

	return mediaRho;
}

// -----------------------------------------------------------------------
double variancia_rho(struct generation *host_geracao,
							struct parametros *host_parametrosGA) {
// -----------------------------------------------------------------------
	double s2 = 0;
	double media_rho;
	double somatoria = 0;
	double rho_i = 0;
	double rho_menos_media = 0;
	double rho_menos_media_ao_quadrado = 0;

	unsigned short int iIndividuos;

	media_rho = mediaRho(host_geracao, host_parametrosGA);

	for (iIndividuos = 0; iIndividuos < host_parametrosGA->numIndividuos; iIndividuos++) {
		rho_i = (double) host_geracao->individuo[iIndividuos].cociente_Rayleigh;
		rho_menos_media = rho_i - media_rho;
		rho_menos_media_ao_quadrado = rho_menos_media*rho_menos_media;
		somatoria = somatoria + rho_menos_media_ao_quadrado;
	}

	s2 = somatoria / (host_parametrosGA->numIndividuos - 1);

	return s2;

}

// -----------------------------------------------------------------------
double desvio_padrao_rho(	struct generation *host_geracao,
									struct parametros *host_parametrosGA) {
// -----------------------------------------------------------------------
	double s = -1;
	double media = -1;
	double variancia = -1;

	media = mediaRho(host_geracao, host_parametrosGA);
	variancia = variancia_rho(host_geracao, host_parametrosGA);
	s = sqrt(variancia);

	return s;

}

// -----------------------------------------------------------------------
double desvio_da_media_rho(struct generation *host_geracao,
									struct parametros *host_parametrosGA) {
// -----------------------------------------------------------------------
	double dm = -1;
	double raiz_de_n = -1;
	double desvio_padrao = -1;
	
	desvio_padrao = desvio_padrao_rho(host_geracao, host_parametrosGA);
	raiz_de_n = sqrt((double)host_parametrosGA->numIndividuos);
	dm = desvio_padrao / raiz_de_n;

	return dm;
}

// -----------------------------------------------------------------------
void rho_minimo_e_maximo(	float *rho_minimo,
									float *rho_maximo,
									struct generation *host_geracao,
									struct parametros *host_parametrosGA) {
// -----------------------------------------------------------------------

	float minimo = 99999.9F;
	float maximo = -99999.9F;
	float rho_i = -1;

	unsigned short int iIndividuo;

	for (iIndividuo = 0; iIndividuo < host_parametrosGA->numIndividuos ; iIndividuo++) {
		
		rho_i = host_geracao->individuo[iIndividuo].cociente_Rayleigh;
		
		if (rho_i <= minimo)
			minimo = rho_i;

		if (rho_i >= maximo) 
			maximo = rho_i;

	}

	*rho_minimo = minimo;
	*rho_maximo = maximo;
}


// -----------------------------------------------------------------------
float amplitude_total_rho(float *rho_minimo,
									float *rho_maximo,
									struct generation *host_geracao,
									struct parametros *host_parametrosGA) {
// -----------------------------------------------------------------------
	
	rho_minimo_e_maximo(rho_minimo, rho_maximo, host_geracao, host_parametrosGA);
	return *rho_maximo - *rho_minimo;
}

// -----------------------------------------------------------------------
void gravaEstatistica(	unsigned long int tempoTotalProcessamento,
								unsigned long int semente,
								unsigned long int iGeracao,
								struct generation *Geracao,
								struct parametrosPrograma *parmsPrograma,
								struct parametros *parametrosGA,
								struct parametros_Metodo *parametrosMetodo) {
// -----------------------------------------------------------------------

	FILE *arqEstatistica;

	if ( (arqEstatistica = fopen("estatistica.txt", "a")) == NULL) {
		printf("\n"); printf("Arquivo estatistica.txt nao encontrato.");
		printf("\n"); printf("Criando novo arquivo estatistica.txt ...");
		if ( (arqEstatistica = fopen("estatistica.txt", "w+")) == NULL) {
			printf("\n"); printf("Erro na criacao do arquivo estatistica.txt.");
		}
		else {
			printf("\n"); printf("Ok.");
		}
	}

	if (arqEstatistica != NULL) {
		
		float rho_minimo = -1;
		float rho_maximo = -1;
		rho_minimo_e_maximo(&rho_minimo, &rho_maximo, Geracao, parametrosGA);
		
		fprintf(arqEstatistica, "%d\t", semente);
		fprintf(arqEstatistica, "%d\t", parmsPrograma->parmMaquina); // 1
		fprintf(arqEstatistica, "%d\t", parmsPrograma->parmSerial_ou_Paralelo); // 1
		fprintf(arqEstatistica, "%d\t", parametrosGA->numGenes); // 1
		fprintf(arqEstatistica, "%d\t", parametrosGA->numIndividuos); // 2
		fprintf(arqEstatistica, "%d\t", parmsPrograma->parmTipoFitness); // 1
		fprintf(arqEstatistica, "%1.6f\t", parametrosGA->probCrossOver); // 3
		fprintf(arqEstatistica, "%1.6f\t", parametrosGA->probMutacao); // 4
		fprintf(arqEstatistica, "%1.6f\t", parametrosGA->intensidadeMutacao); // 5
		fprintf(arqEstatistica, "%1.6f\t", parametrosMetodo->lambda); // 6
		fprintf(arqEstatistica, "%1.6f\t", parmsPrograma->parmTolerancia);
		fprintf(arqEstatistica, "%1.6f\t", Geracao->gradRhoMedio);
		fprintf(arqEstatistica, "%1.6f\t", Geracao->difRho);
		fprintf(arqEstatistica, "%1.6f\t", Geracao->rhoMedio);
		fprintf(arqEstatistica, "%1.6f\t", parmsPrograma->parmRhoMinimo); // rho_0, limite inferior, (rho-rho_0) do fitness.
		fprintf(arqEstatistica, "%d\t", iGeracao);
		fprintf(arqEstatistica, "%1.6f\t", Geracao->FitnessMedio);
		fprintf(arqEstatistica, "%1.6f\t", Geracao->Maior_fitness);
		fprintf(arqEstatistica, "%1.6f\t", Geracao->Melhor_cociente_Rayleigh);
		fprintf(arqEstatistica, "%3.20f\t", mediaRho(Geracao, parametrosGA)); // 8
		fprintf(arqEstatistica, "%3.20f\t", desvio_padrao_rho(Geracao, parametrosGA)); // 9
		fprintf(arqEstatistica, "%3.20f\t", desvio_da_media_rho(Geracao, parametrosGA));// 10
		fprintf(arqEstatistica, "%3.20f\t", rho_minimo); // esse é o menor rho da última geração
		fprintf(arqEstatistica, "%3.20f\t", rho_maximo);// 12
		fprintf(arqEstatistica, "%3.20f\t", amplitude_total_rho(&rho_minimo, &rho_maximo, Geracao, parametrosGA));// 13
		fprintf(arqEstatistica, "%d\n", tempoTotalProcessamento);
		
		fclose(arqEstatistica);		
	}
}

// -----------------------------------------------------------------------
float erroAbsolutoEntre2(float num1, float num2) {
// -----------------------------------------------------------------------
	return abs(num1 - num2);
}