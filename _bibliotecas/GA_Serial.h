

//=========================================================================================================
void Calcula_Quociente_Rayleigh(	float *VetorC,
											float *Hamiltoniano,
											unsigned short int Ordem_do_Hamiltoniano,
											float *rho) {
//=========================================================================================================

	float *Matriz_Auxiliar, *Numerador, *Denominador;

/*
	rho = [Ct * H * C  ]  / [ Ct * C ] 

	onde
		C = vetor C (candidato à autovetor).
		Ct = transposto do vetor C.
		H = hamiltoniano.g
*/

	Matriz_Auxiliar	= (float *)malloc(Ordem_do_Hamiltoniano * sizeof(float));
	Numerador			= (float *)malloc(sizeof(float));
	Denominador			= (float *)malloc(sizeof(float));

	// Numerador [Ct * H * C  ]
	Multiplica_Matrizes(VetorC, 1, Ordem_do_Hamiltoniano, Hamiltoniano, Ordem_do_Hamiltoniano, Ordem_do_Hamiltoniano, Matriz_Auxiliar);
	Multiplica_Matrizes(Matriz_Auxiliar, 1, Ordem_do_Hamiltoniano, VetorC, Ordem_do_Hamiltoniano, 1, Numerador);

	// Denominador [ Ct * C ]
	Multiplica_Matrizes(VetorC, 1, Ordem_do_Hamiltoniano, VetorC, Ordem_do_Hamiltoniano, 1, Denominador);

	*rho = (float)(*Numerador / *Denominador);

	free(Matriz_Auxiliar);
	free(Numerador);
	free(Denominador);
}


//=========================================================================================================
void gera_Matriz_de_Coope_abs(float *H, unsigned short int ordemH) {
//=========================================================================================================
	unsigned short int iLinha;
	unsigned long int qtdeElementos, iElemento_diagonal, iElemento;

	qtdeElementos = ordemH * ordemH;
	
	for (iElemento = 0; iElemento < qtdeElementos; iElemento++) {
		iLinha = (unsigned short int)(iElemento/ordemH);
		iElemento_diagonal = (1 + ordemH)*iLinha;
		
		if (iElemento == iElemento_diagonal)
			*(H + iElemento) = (float)(2*(iLinha+1) - 1);
		else 
			*(H + iElemento) = 1;
	}
}

//=========================================================================================================
void testa_gera_Matriz_de_Coope_abs(void) {
//=========================================================================================================

	unsigned short int ordem_H = 50;
	unsigned long int qtde_Elementos = ordem_H * ordem_H;
	float *H;

	H = (float *)malloc(qtde_Elementos * sizeof(float));
	
	gera_Matriz_de_Coope_abs(H, ordem_H);

	printf("\n"); printf("H = ");
	imprimeMatriz(H, ordem_H, ordem_H);

}

//=========================================================================================================
void testa_Calcula_Quociente_Rayleigh(void) {
//=========================================================================================================
	
	/* #### TESTE COM HAMILTONIANO ##### */
	unsigned short int i, ordemH = 50;
	unsigned long int qtdeElementos;

	qtdeElementos = ordemH * ordemH;

	float *H, *C, E, E_esperado;
	H = (float *)malloc(qtdeElementos*sizeof(float));
	C = (float *)malloc(ordemH*sizeof(float));

	gera_Matriz_de_Coope_abs(H, ordemH);
	
	// Esse é o primeiro autovetor, com autovalor = 0,296280
	/*
	E_esperado = 0.296280F;
	*(C + 0) = 0.9780472F;
	*(C + 1) = -0.1700842F;
	*(C + 2) = -0.0782391F;
	*(C + 3) = -0.0508047F;
	*(C + 4) = -0.0376151F;
	*(C + 5) = -0.0298623F;
	*(C + 6) = -0.0247593F;
	*(C + 7) = -0.0211458F;
	*(C + 8) = -0.0184527F;
	*(C + 9) = -0.0163681F;
	*(C + 10) = -0.0147067F;
	*(C + 11) = -0.0133514F;
	*(C + 12) = -0.0122249F;
	*(C + 13) = -0.0112737F;
	*(C + 14) = -0.0104598F;
	*(C + 15) = -0.0097555F;
	*(C + 16) = -0.0091401F;
	*(C + 17) = -0.0085977F;
	*(C + 18) = -0.0081161F;
	*(C + 19) = -0.0076856F;
	*(C + 20) = -0.0072985F;
	*(C + 21) = -0.0069484F;
	*(C + 22) = -0.0066305F;
	*(C + 23) = -0.0063403F;
	*(C + 24) = -0.0060745F;
	*(C + 25) = -0.0058301F;
	*(C + 26) = -0.0056045F;
	*(C + 27) = -0.0053958F;
	*(C + 28) = -0.0052021F;
	*(C + 29) = -0.0050218F;
	*(C + 30) = -0.0048536F;
	*(C + 31) = -0.0046962F;
	*(C + 32) = -0.0045488F;
	*(C + 33) = -0.0044103F;
	*(C + 34) = -0.0042801F;
	*(C + 35) = -0.0041573F;
	*(C + 36) = -0.0040413F;
	*(C + 37) = -0.0039316F;
	*(C + 38) = -0.0038278F;
	*(C + 39) = -0.0037292F;
	*(C + 40) = -0.0036357F;
	*(C + 41) = -0.0035467F;
	*(C + 42) = -0.0034619F;
	*(C + 43) = -0.0033811F;
	*(C + 44) = -0.003304F;
	*(C + 45) = -0.0032304F;
	*(C + 46) = -0.0031599F;
	*(C + 47) = -0.0030925F;
	*(C + 48) = -0.0030278F;
	*(C + 49) = -0.0029659F;
	*/

	// Segundo autovetor, com autovalor = 2,337933
	/*
	E_esperado = 2.337933F;
	*(C + 0) = -0.1389437F;
	*(C + 1) = -0.9612598F;
	*(C + 2) = 0.1954439F;
	*(C + 3) = 0.0887042F;
	*(C + 4) = 0.0573714F;
	*(C + 5) = 0.042396F;
	*(C + 6) = 0.0336202F;
	*(C + 7) = 0.0278545F;
	*(C + 8) = 0.0237768F;
	*(C + 9) = 0.0207406F;
	*(C + 10) = 0.018392F;
	*(C + 11) = 0.0165212F;
	*(C + 12) = 0.0149958F;
	*(C + 13) = 0.0137283F;
	*(C + 14) = 0.0126584F;
	*(C + 15) = 0.0117432F;
	*(C + 16) = 0.0109514F;
	*(C + 17) = 0.0102596F;
	*(C + 18) = 0.0096501F;
	*(C + 19) = 0.0091089F;
	*(C + 20) = 0.0086251F;
	*(C + 21) = 0.0081902F;
	*(C + 22) = 0.007797F;
	*(C + 23) = 0.0074399F;
	*(C + 24) = 0.007114F;
	*(C + 25) = 0.0068155F;
	*(C + 26) = 0.006541F;
	*(C + 27) = 0.0062878F;
	*(C + 28) = 0.0060535F;
	*(C + 29) = 0.0058359F;
	*(C + 30) = 0.0056335F;
	*(C + 31) = 0.0054447F;
	*(C + 32) = 0.0052681F;
	*(C + 33) = 0.0051026F;
	*(C + 34) = 0.0049472F;
	*(C + 35) = 0.0048009F;
	*(C + 36) = 0.0046631F;
	*(C + 37) = 0.004533F;
	*(C + 38) = 0.0044099F;
	*(C + 39) = 0.0042933F;
	*(C + 40) = 0.0041827F;
	*(C + 41) = 0.0040777F;
	*(C + 42) = 0.0039779F;
	*(C + 43) = 0.0038828F;
	*(C + 44) = 0.0037921F;
	*(C + 45) = 0.0037056F;
	*(C + 46) = 0.0036229F;
	*(C + 47) = 0.0035439F;
	*(C + 48) = 0.0034682F;
	*(C + 49) = 0.0033957F;
	*/
	
	// Terceiro autovetor, com autovalor = 4,365059
	E_esperado = 4.365059F;
	*(C + 0) = -0.0795874F;
	*(C + 1) = -0.14689F;
	*(C + 2) = -0.9516369F;
	*(C + 3) = 0.2124869F;
	*(C + 4) = 0.0955734F;
	*(C + 5) = 0.0616517F;
	*(C + 6) = 0.0455018F;
	*(C + 7) = 0.0360566F;
	*(C + 8) = 0.0298586F;
	*(C + 9) = 0.0254789F;
	*(C + 10) = 0.0222197F;
	*(C + 11) = 0.0196997F;
	*(C + 12) = 0.0176931F;
	*(C + 13) = 0.0160575F;
	*(C + 14) = 0.0146987F;
	*(C + 15) = 0.013552F;
	*(C + 16) = 0.0125712F;
	*(C + 17) = 0.0117228F;
	*(C + 18) = 0.0109816F;
	*(C + 19) = 0.0103287F;
	*(C + 20) = 0.009749F;
	*(C + 21) = 0.0092309F;
	*(C + 22) = 0.0087651F;
	*(C + 23) = 0.008344F;
	*(C + 24) = 0.0079616F;
	*(C + 25) = 0.0076127F;
	*(C + 26) = 0.007293F;
	*(C + 27) = 0.0069992F;
	*(C + 28) = 0.0067281F;
	*(C + 29) = 0.0064772F;
	*(C + 30) = 0.0062443F;
	*(C + 31) = 0.0060277F;
	*(C + 32) = 0.0058255F;
	*(C + 33) = 0.0056365F;
	*(C + 34) = 0.0054593F;
	*(C + 35) = 0.005293F;
	*(C + 36) = 0.0051365F;
	*(C + 37) = 0.0049889F;
	*(C + 38) = 0.0048496F;
	*(C + 39) = 0.0047179F;
	*(C + 40) = 0.0045932F;
	*(C + 41) = 0.0044748F;
	*(C + 42) = 0.0043625F;
	*(C + 43) = 0.0042556F;
	*(C + 44) = 0.0041538F;
	*(C + 45) = 0.0040568F;
	*(C + 46) = 0.0039642F;
	*(C + 47) = 0.0038758F;
	*(C + 48) = 0.0037912F;
	*(C + 49) = 0.0037102F;

	// E se o vetor foi um mútiplo do autovetor?
	float fator = 13.7983F;
	for (i = 0; i < ordemH; i++) {
		*(C + i) = *(C + i)*fator;
	}

	Calcula_Quociente_Rayleigh(C, H, ordemH, &E);
	printf("\n"); printf("Rho esperado = %f", E_esperado);
	printf("\n"); printf("Rho obtido   = %f", E);

	/* #### TESTE COM MATRIZ SIMPLES #####

	float *A, *w1, *w2;
	
	float e1 = 3.0F;
	float e2 = -1.0F;

	float e1_pgm = 0.0F;
	float e2_pgm = 0.0F;

	unsigned short int ordemMatriz = 2;

	A = (float *)malloc(ordemMatriz*ordemMatriz*sizeof(float));
	w1 = (float *)malloc(ordemMatriz*sizeof(float));
	w2 = (float *)malloc(ordemMatriz*sizeof(float));

	*(A + 0) = +1.0F;
	*(A + 1) = -1.0F;
	*(A + 2) = -4.0F;
	*(A + 3) = +1.0F;

	*(w1 + 0) = +1.0F;
	*(w1 + 1) = -2.0F;

	*(w2 + 0) = +1.0F;
	*(w2 + 1) = +2.0F;

	// Teset para o primeiro autovalor / autovetor
	Calcula_Quociente_Rayleigh(w1, A, ordemMatriz, &e1_pgm);
	printf("\n"); printf("autovalor e1 = %f", e1);
	printf("\n"); printf("autovalor e1_pgm = %f", e1_pgm);

	Calcula_Quociente_Rayleigh(w2, A, ordemMatriz, &e2_pgm);
	printf("\n"); printf("autovalor e1 = %f", e2);
	printf("\n"); printf("autovalor e1_pgm = %f", e2_pgm);
	*/
}

//=========================================================================================================
void Gradiente_de_Rho(float *Hamiltoniano, const unsigned short int Ordem_Hamiltoniano,
							 float *Autovetor_C, float Cociente_Rayleigh, float *Gradiente_de_rho,
							 float *mIdentidade) {
//=========================================================================================================
	
	float *mAuxOxO, *mAuxOx1, CtC;
	unsigned long int	qtdeElementos;
	//unsigned short int iElemento;
	//float *mIdentidade;

	qtdeElementos = Ordem_Hamiltoniano*Ordem_Hamiltoniano;

	// Cria matriz identidade
	
	/*
	mIdentidade = (float *)malloc(qtdeElementos*sizeof(float));
	for (iElemento = 0; iElemento < qtdeElementos; iElemento++) {
		if ( (iElemento + Ordem_Hamiltoniano + 1) % (Ordem_Hamiltoniano + 1) == 0) {
			*(mIdentidade + iElemento) = 1.0F;
		}
		else {
			*(mIdentidade + iElemento) = 0.0F;
		}
	}
	*/

	mAuxOxO = (float *)malloc(qtdeElementos*sizeof(float));
	multiplica_matriz_por_escalar(mIdentidade, Cociente_Rayleigh, qtdeElementos, mAuxOxO);

	// TESTE INÍCIO (pode ser comentado)
	/*
	printf("\n\nMultiplicando a matriz por rho = %f\n", Cociente_Rayleigh);
	for (iElemento = 0; iElemento < qtdeElementos; iElemento++)
		printf("\na[%d] = %f", iElemento, *(mAuxOxO + iElemento));
	*/
	// TESTE FIM

	// TESTE INÍCIO (pode ser comentado)
	/*
	printf("\n\nHamiltoniano\n");
	for (iElemento = 0; iElemento < qtdeElementos; iElemento++)
		printf("\na[%d] = %f", iElemento, *(Hamiltoniano + iElemento));
	*/
	// TESTE FIM

	Subtrai_Matrizes(Hamiltoniano, mAuxOxO, qtdeElementos, mAuxOxO);
	
	// TESTE INÍCIO (pode ser comentado)
	/*
	printf("\n\nHamiltoniano - rho.Indentidade\n");
	for (iElemento = 0; iElemento < qtdeElementos; iElemento++)
		printf("\na[%d] = %f", iElemento, *(mAuxOxO + iElemento));
	*/
	// TESTE FIM

	multiplica_matriz_por_escalar(mAuxOxO, 2, qtdeElementos, mAuxOxO);

	// TESTE INÍCIO (pode ser comentado)
	/*
	printf("\n\nMultiplicando tudo por 2 (vem da equacao)\n");
	for (iElemento = 0; iElemento < qtdeElementos; iElemento++)
		printf("\na[%d] = %f", iElemento, *(mAuxOxO + iElemento));
	*/
	// TESTE FIM

	// TESTE INÍCIO (pode ser comentado)
	/*
	printf("\n\nC (candidato a autovetor)\n");
	for (iElemento = 0; iElemento < Ordem_Hamiltoniano; iElemento++)
		printf("\na[%d] = %f", iElemento, *(Autovetor_C + iElemento));
	*/
	// TESTE FIM

	mAuxOx1 = (float *)malloc(Ordem_Hamiltoniano*sizeof(float));
	Multiplica_Matrizes(mAuxOxO, Ordem_Hamiltoniano, Ordem_Hamiltoniano, Autovetor_C, Ordem_Hamiltoniano, 1, mAuxOx1);

	// TESTE INÍCIO (pode ser comentado)
	/*
	printf("\n\nMultiplicando tudo por C (numerador do gradiente de rho)\n");
	for (iElemento = 0; iElemento < Ordem_Hamiltoniano; iElemento++)
		printf("\na[%d] = %f", iElemento, *(mAuxOx1 + iElemento));
	*/
	// TESTE FIM

	Multiplica_Matrizes(Autovetor_C, 1, Ordem_Hamiltoniano, Autovetor_C, Ordem_Hamiltoniano, 1, &CtC);

	// TESTE INÍCIO (pode ser comentado)
	/*
	printf("\n\nC_transposto x C\n");
	printf("\n%f", CtC);
	*/
	// TESTE FIM

	multiplica_matriz_por_escalar(mAuxOx1, (float)(1/CtC), Ordem_Hamiltoniano, Gradiente_de_rho);

	// TESTE INÍCIO (pode ser comentado)
	/*
	printf("\n\nGradiente de rho na função Gradiente_de_Rho()\n");
	for (iElemento = 0; iElemento < Ordem_Hamiltoniano; iElemento++)
		printf("\na[%d] = %f", iElemento, *(Gradiente_de_rho + iElemento));
	*/
	// TESTE FIM

	//free(mIdentidade);
	//mIdentidade = NULL;
	free(mAuxOxO);
	mAuxOxO = NULL;
	free(mAuxOx1);
	mAuxOx1 = NULL;
}

//=========================================================================================================
void testa_Gradiente_de_Rho(void) {
//=========================================================================================================

	/* #### TESTE COM HAMILTONIANO ##### */
	unsigned short int i, ordemH = 50;
	unsigned long int qtdeElementos;

	qtdeElementos = ordemH * ordemH;

	float *H, *C, *grad_rho, E;
	H = (float *)malloc(qtdeElementos*sizeof(float));
	C = (float *)malloc(ordemH*sizeof(float));
	grad_rho = (float *)malloc(ordemH*sizeof(float));

	// Inicializa o gradiente.

	for (i = 0; i < ordemH; i++) {
		*(grad_rho + i) = 99.0F;
	}

	gera_Matriz_de_Coope_abs(H, ordemH);
	
	// Esse é o primeiro autovetor, com autovalor = 0,296280
	E = 0.296280F;
	*(C + 0) = 0.9780472F;
	*(C + 1) = -0.1700842F;
	*(C + 2) = -0.0782391F;
	*(C + 3) = -0.0508047F;
	*(C + 4) = -0.0376151F;
	*(C + 5) = -0.0298623F;
	*(C + 6) = -0.0247593F;
	*(C + 7) = -0.0211458F;
	*(C + 8) = -0.0184527F;
	*(C + 9) = -0.0163681F;
	*(C + 10) = -0.0147067F;
	*(C + 11) = -0.0133514F;
	*(C + 12) = -0.0122249F;
	*(C + 13) = -0.0112737F;
	*(C + 14) = -0.0104598F;
	*(C + 15) = -0.0097555F;
	*(C + 16) = -0.0091401F;
	*(C + 17) = -0.0085977F;
	*(C + 18) = -0.0081161F;
	*(C + 19) = -0.0076856F;
	*(C + 20) = -0.0072985F;
	*(C + 21) = -0.0069484F;
	*(C + 22) = -0.0066305F;
	*(C + 23) = -0.0063403F;
	*(C + 24) = -0.0060745F;
	*(C + 25) = -0.0058301F;
	*(C + 26) = -0.0056045F;
	*(C + 27) = -0.0053958F;
	*(C + 28) = -0.0052021F;
	*(C + 29) = -0.0050218F;
	*(C + 30) = -0.0048536F;
	*(C + 31) = -0.0046962F;
	*(C + 32) = -0.0045488F;
	*(C + 33) = -0.0044103F;
	*(C + 34) = -0.0042801F;
	*(C + 35) = -0.0041573F;
	*(C + 36) = -0.0040413F;
	*(C + 37) = -0.0039316F;
	*(C + 38) = -0.0038278F;
	*(C + 39) = -0.0037292F;
	*(C + 40) = -0.0036357F;
	*(C + 41) = -0.0035467F;
	*(C + 42) = -0.0034619F;
	*(C + 43) = -0.0033811F;
	*(C + 44) = -0.003304F;
	*(C + 45) = -0.0032304F;
	*(C + 46) = -0.0031599F;
	*(C + 47) = -0.0030925F;
	*(C + 48) = -0.0030278F;
	*(C + 49) = -0.0029659F;
	

	// Segundo autovetor, com autovalor = 2,337933
	/*
	*(C + 0) = -0.1389437F;
	*(C + 1) = -0.9612598F;
	*(C + 2) = 0.1954439F;
	*(C + 3) = 0.0887042F;
	*(C + 4) = 0.0573714F;
	*(C + 5) = 0.042396F;
	*(C + 6) = 0.0336202F;
	*(C + 7) = 0.0278545F;
	*(C + 8) = 0.0237768F;
	*(C + 9) = 0.0207406F;
	*(C + 10) = 0.018392F;
	*(C + 11) = 0.0165212F;
	*(C + 12) = 0.0149958F;
	*(C + 13) = 0.0137283F;
	*(C + 14) = 0.0126584F;
	*(C + 15) = 0.0117432F;
	*(C + 16) = 0.0109514F;
	*(C + 17) = 0.0102596F;
	*(C + 18) = 0.0096501F;
	*(C + 19) = 0.0091089F;
	*(C + 20) = 0.0086251F;
	*(C + 21) = 0.0081902F;
	*(C + 22) = 0.007797F;
	*(C + 23) = 0.0074399F;
	*(C + 24) = 0.007114F;
	*(C + 25) = 0.0068155F;
	*(C + 26) = 0.006541F;
	*(C + 27) = 0.0062878F;
	*(C + 28) = 0.0060535F;
	*(C + 29) = 0.0058359F;
	*(C + 30) = 0.0056335F;
	*(C + 31) = 0.0054447F;
	*(C + 32) = 0.0052681F;
	*(C + 33) = 0.0051026F;
	*(C + 34) = 0.0049472F;
	*(C + 35) = 0.0048009F;
	*(C + 36) = 0.0046631F;
	*(C + 37) = 0.004533F;
	*(C + 38) = 0.0044099F;
	*(C + 39) = 0.0042933F;
	*(C + 40) = 0.0041827F;
	*(C + 41) = 0.0040777F;
	*(C + 42) = 0.0039779F;
	*(C + 43) = 0.0038828F;
	*(C + 44) = 0.0037921F;
	*(C + 45) = 0.0037056F;
	*(C + 46) = 0.0036229F;
	*(C + 47) = 0.0035439F;
	*(C + 48) = 0.0034682F;
	*(C + 49) = 0.0033957F;
	*/
	
	// Terceiro autovetor, com autovalor = 4,365059
	/*
	E = 4.365059F;
	*(C + 0) = -0.0795874F;
	*(C + 1) = -0.14689F;
	*(C + 2) = -0.9516369F;
	*(C + 3) = 0.2124869F;
	*(C + 4) = 0.0955734F;
	*(C + 5) = 0.0616517F;
	*(C + 6) = 0.0455018F;
	*(C + 7) = 0.0360566F;
	*(C + 8) = 0.0298586F;
	*(C + 9) = 0.0254789F;
	*(C + 10) = 0.0222197F;
	*(C + 11) = 0.0196997F;
	*(C + 12) = 0.0176931F;
	*(C + 13) = 0.0160575F;
	*(C + 14) = 0.0146987F;
	*(C + 15) = 0.013552F;
	*(C + 16) = 0.0125712F;
	*(C + 17) = 0.0117228F;
	*(C + 18) = 0.0109816F;
	*(C + 19) = 0.0103287F;
	*(C + 20) = 0.009749F;
	*(C + 21) = 0.0092309F;
	*(C + 22) = 0.0087651F;
	*(C + 23) = 0.008344F;
	*(C + 24) = 0.0079616F;
	*(C + 25) = 0.0076127F;
	*(C + 26) = 0.007293F;
	*(C + 27) = 0.0069992F;
	*(C + 28) = 0.0067281F;
	*(C + 29) = 0.0064772F;
	*(C + 30) = 0.0062443F;
	*(C + 31) = 0.0060277F;
	*(C + 32) = 0.0058255F;
	*(C + 33) = 0.0056365F;
	*(C + 34) = 0.0054593F;
	*(C + 35) = 0.005293F;
	*(C + 36) = 0.0051365F;
	*(C + 37) = 0.0049889F;
	*(C + 38) = 0.0048496F;
	*(C + 39) = 0.0047179F;
	*(C + 40) = 0.0045932F;
	*(C + 41) = 0.0044748F;
	*(C + 42) = 0.0043625F;
	*(C + 43) = 0.0042556F;
	*(C + 44) = 0.0041538F;
	*(C + 45) = 0.0040568F;
	*(C + 46) = 0.0039642F;
	*(C + 47) = 0.0038758F;
	*(C + 48) = 0.0037912F;
	*(C + 49) = 0.0037102F;

	*/

	// E se o vetor foi um mútiplo do autovetor?
	float fator = 13.000000F;
	for (i = 0; i < ordemH; i++) {
		*(C + i) = *(C + i)*fator;
	}

	float *matriz_Identidade;
	matriz_Identidade = (float *)malloc(ordemH*ordemH*sizeof(float));
	GeraMatrizIdentidade(matriz_Identidade, ordemH);
	Gradiente_de_Rho(H, ordemH, C, E, grad_rho, matriz_Identidade);
	free(matriz_Identidade);

	printf("\n"); printf("Gradiente de %f = ", E);
	imprimeMatriz(grad_rho, ordemH, 1);

	/*

	TESTE SIMPLES, COM MATRIZ 2x2

	float *A, *w1, *w2;
	
	float e1 = 3.0F;
	float e2 = -1.0F;

	float e1_pgm = 0.0F;
	float e2_pgm = 0.0F;

	unsigned short int ordemMatriz = 2;

	A = (float *)malloc(ordemMatriz*ordemMatriz*sizeof(float));
	w1 = (float *)malloc(ordemMatriz*sizeof(float));
	w2 = (float *)malloc(ordemMatriz*sizeof(float));

	*(A + 0) = +1.0F;
	*(A + 1) = -1.0F;
	*(A + 2) = -4.0F;
	*(A + 3) = +1.0F;

	*(w1 + 0) = +1.0F;
	*(w1 + 1) = -2.0F;

	*(w2 + 0) = +1.0F;
	*(w2 + 1) = +2.0F;

	float *grad_e1;
	float *grad_e2;

	grad_e1 = (float *)malloc(ordemMatriz*sizeof(float));
	grad_e2 = (float *)malloc(ordemMatriz*sizeof(float));

	unsigned short int i = 0;

	// grad_e1
	for (i = 0; i < ordemMatriz; i++) {
		*(grad_e1 + i) = 99.0F;
	}
	Gradiente_de_Rho(A, ordemMatriz, w1, e1, grad_e1);
	printf("\n"); printf("grad_e1 =");
	imprimeMatriz(grad_e1, ordemMatriz, 1);

	//for (i = 0; i < ordemMatriz; i++) {
	//	printf("\n"); printf("%f", *(grad_e1 + i));
	//}

	// grad_e1
	for (i = 0; i < ordemMatriz; i++) {
		*(grad_e2 + i) = 99.0F;
	}

	Gradiente_de_Rho(A, ordemMatriz, w2, e2, grad_e2);

	printf("\n"); printf("grad_e2 =");
	for (i = 0; i < ordemMatriz; i++) {
		printf("\n"); printf("%f", *(grad_e2 + i));
	}
	*/

}



//=========================================================================================================
void Gera_Matriz_de_Coope_LinhaColuna(float* Matriz, unsigned short int Ordem) {
//=========================================================================================================
	unsigned short int iLinha, iColuna;

	for (iLinha = 0; iLinha < Ordem; iLinha++) {
		for (iColuna = 0; iColuna < Ordem; iColuna++) {
			if (iLinha == iColuna) {
				*(Matriz + Ordem*iLinha + iColuna) = (float)(2*(iLinha+1) - 1);
				//*(Matriz + Ordem*iLinha + iColuna) = (float) -1 * (20 - iLinha);
					// a definição acima gera uma matriz com todos os autovalores negativos.
			}
			else {
				*(Matriz + Ordem*iLinha + iColuna) = 1;
			}
		}
	}
}


//=========================================================================================================
void testa_Gera_Matriz_de_Coope_LinhaColuna(void) {
//=========================================================================================================
	unsigned short int ordem_H = 50;
	unsigned long int qtde_Elementos = ordem_H * ordem_H;
	float *H;

	H = (float *)malloc(qtde_Elementos * sizeof(float));
	
	Gera_Matriz_de_Coope_LinhaColuna(H, ordem_H);

	printf("\n"); printf("H = ");
	imprimeMatriz(H, ordem_H, ordem_H);
}



//=========================================================================================================
void Fitness_Serial(
		unsigned short int tipoFitness,
		struct generation *geracao,
		struct parametros *parametrosGA,
		struct parametros_Metodo *parms_Metodo,
		float *hamiltoniano,
		const unsigned short int ordem_hamiltoniano,
		float lambda,
		float rho_minimo,
		float *matriz_Identidade) {
//=========================================================================================================

	unsigned short int	iIndividuo;
	//unsigned short int		iGene;

	geracao->rhoMedio = -1.0F;
	geracao->sumRho = 0.0F;
	geracao->FitnessMedio = -1.0F;
	geracao->Maior_fitness = -1.0F;
	geracao->sumFitness = 0.0F;
	geracao->Maior_fitness = -1.0F;
	geracao->Melhor_cociente_Rayleigh = -1.0F;

	geracao->sumGradRho = 0.0F;
	geracao->gradRhoMedio = -1.0F;

	for (iIndividuo = 0; iIndividuo < parametrosGA->numIndividuos; iIndividuo++) {
		
		Calcula_Quociente_Rayleigh(geracao->individuo[iIndividuo].gene,
											hamiltoniano,
											ordem_hamiltoniano,
											&geracao->individuo[iIndividuo].cociente_Rayleigh);

		geracao->sumRho = geracao->sumRho + geracao->individuo[iIndividuo].cociente_Rayleigh;
		

		Gradiente_de_Rho(	hamiltoniano,
								ordem_hamiltoniano,
								geracao->individuo[iIndividuo].gene,
								geracao->individuo[iIndividuo].cociente_Rayleigh,
								geracao->individuo[iIndividuo].gradRho,
								matriz_Identidade);


		Multiplica_Matrizes(	geracao->individuo[iIndividuo].gradRho, 1, ordem_hamiltoniano,
									geracao->individuo[iIndividuo].gradRho, ordem_hamiltoniano, 1,
									&geracao->individuo[iIndividuo].grad_elevado_ao_quadrado);
		
		geracao->sumGradRho = geracao->sumGradRho + sqrt(geracao->individuo[iIndividuo].grad_elevado_ao_quadrado);

		geracao->individuo[iIndividuo].rho_menos_rho0_ao_quadrado = pow( (geracao->individuo[iIndividuo].cociente_Rayleigh - rho_minimo) ,2);
		
		switch (tipoFitness) {
			case 0: {
				// Só com (rho - rho_0)^2
				geracao->individuo[iIndividuo].fitness = exp( (-1)*lambda*(geracao->individuo[iIndividuo].rho_menos_rho0_ao_quadrado));
			}
			break;
			case 1: {
				// Fitness só com o grad(rho)
				geracao->individuo[iIndividuo].fitness = exp( (-1)*lambda*(geracao->individuo[iIndividuo].grad_elevado_ao_quadrado));
			}
			break;
			case 2: {
				// Com (rho - rho_0)^2 + grad(rho)
				geracao->individuo[iIndividuo].fitness = exp( (-1)*lambda*(geracao->individuo[iIndividuo].rho_menos_rho0_ao_quadrado + geracao->individuo[iIndividuo].grad_elevado_ao_quadrado));
			}
			break;
			default: {
				geracao->individuo[iIndividuo].fitness = -1;
			}
			break;
		}
				
		geracao->sumFitness = geracao->sumFitness + geracao->individuo[iIndividuo].fitness;
		
		if (geracao->individuo[iIndividuo].fitness > geracao->Maior_fitness) {
			geracao->Maior_fitness = geracao->individuo[iIndividuo].fitness;
			geracao->Melhor_cociente_Rayleigh = geracao->individuo[iIndividuo].cociente_Rayleigh;
			geracao->idxMelhorIndividuo = iIndividuo;
		}
	}
	
	geracao->FitnessMedio = geracao->sumFitness / parametrosGA->numIndividuos;
	geracao->gradRhoMedio = geracao->sumGradRho / parametrosGA->numIndividuos;
	geracao->rhoMedio = geracao->sumRho / parametrosGA->numIndividuos;
	geracao->difRho = abs(geracao->rhoMedio - parms_Metodo->rho_minimo);
	

}

//=========================================================================================================
void testa_Fitness_Serial(void) {
//=========================================================================================================
	//TODO: corrigir parâmetros do Fitness_Serial

	/*
	struct generation geracao;
	struct parametros parametrosGA;
	
	parametrosGA.numGenes = 50;
	parametrosGA.numIndividuos = 10;

	float *H;
	H = (float *)malloc(parametrosGA.numGenes*parametrosGA.numGenes*sizeof(float));
	gera_Matriz_de_Coope_abs(H, parametrosGA.numGenes);

	// Inicializa atributos da geracao;
	
	geracao.FitnessMedio = -1.0F;
	geracao.idxMelhorIndividuo = 99;
	geracao.Maior_fitness = -1.0F;;
	geracao.Melhor_cociente_Rayleigh = -1.0F;
	geracao.sumFitness = -1.0F;
	
	for (unsigned short int iIndividuo = 0; iIndividuo < parametrosGA.numIndividuos; iIndividuo++ ) {
		geracao.individuo[iIndividuo].cociente_Rayleigh = -1.0F;
		geracao.individuo[iIndividuo].fitness = -1.0F;
		geracao.individuo[iIndividuo].grad_elevado_ao_quadrado = -1.0;
		geracao.individuo[iIndividuo].CtC = -1.0F;
		geracao.individuo[iIndividuo].inverso_de_CtC = -1.0F;
		geracao.individuo[iIndividuo].Numerador = -1.0F;
		geracao.individuo[iIndividuo].rho_menos_rho0_ao_quadrado = -1.0F;
		for (unsigned short int iGene = 0; iGene < parametrosGA.numGenes; iGene++ ) {
			geracao.individuo[iIndividuo].gene[iGene] = -1.0F;			
		}
	}

	printf("\n"); printf("Geracao sem autovalores = ");
	imprimeGeracao(&geracao, &parametrosGA);
	
	// ------------------------------
	//	INICIALIZAÇÃO DOS VETORES COM OS GENES DOS AUTOVETORES
	//   -------------------------------                         

	// Indivíduo 0
	geracao.individuo[0].gene[0] = 0.9780472F;
	geracao.individuo[0].gene[1] = -0.1700842F;
	geracao.individuo[0].gene[2] = -0.0782391F;
	geracao.individuo[0].gene[3] = -0.0508047F;
	geracao.individuo[0].gene[4] = -0.0376151F;
	geracao.individuo[0].gene[5] = -0.0298623F;
	geracao.individuo[0].gene[6] = -0.0247593F;
	geracao.individuo[0].gene[7] = -0.0211458F;
	geracao.individuo[0].gene[8] = -0.0184527F;
	geracao.individuo[0].gene[9] = -0.0163681F;
	geracao.individuo[0].gene[10] = -0.0147067F;
	geracao.individuo[0].gene[11] = -0.0133514F;
	geracao.individuo[0].gene[12] = -0.0122249F;
	geracao.individuo[0].gene[13] = -0.0112737F;
	geracao.individuo[0].gene[14] = -0.0104598F;
	geracao.individuo[0].gene[15] = -0.0097555F;
	geracao.individuo[0].gene[16] = -0.0091401F;
	geracao.individuo[0].gene[17] = -0.0085977F;
	geracao.individuo[0].gene[18] = -0.0081161F;
	geracao.individuo[0].gene[19] = -0.0076856F;
	geracao.individuo[0].gene[20] = -0.0072985F;
	geracao.individuo[0].gene[21] = -0.0069484F;
	geracao.individuo[0].gene[22] = -0.0066305F;
	geracao.individuo[0].gene[23] = -0.0063403F;
	geracao.individuo[0].gene[24] = -0.0060745F;
	geracao.individuo[0].gene[25] = -0.0058301F;
	geracao.individuo[0].gene[26] = -0.0056045F;
	geracao.individuo[0].gene[27] = -0.0053958F;
	geracao.individuo[0].gene[28] = -0.0052021F;
	geracao.individuo[0].gene[29] = -0.0050218F;
	geracao.individuo[0].gene[30] = -0.0048536F;
	geracao.individuo[0].gene[31] = -0.0046962F;
	geracao.individuo[0].gene[32] = -0.0045488F;
	geracao.individuo[0].gene[33] = -0.0044103F;
	geracao.individuo[0].gene[34] = -0.0042801F;
	geracao.individuo[0].gene[35] = -0.0041573F;
	geracao.individuo[0].gene[36] = -0.0040413F;
	geracao.individuo[0].gene[37] = -0.0039316F;
	geracao.individuo[0].gene[38] = -0.0038278F;
	geracao.individuo[0].gene[39] = -0.0037292F;
	geracao.individuo[0].gene[40] = -0.0036357F;
	geracao.individuo[0].gene[41] = -0.0035467F;
	geracao.individuo[0].gene[42] = -0.0034619F;
	geracao.individuo[0].gene[43] = -0.0033811F;
	geracao.individuo[0].gene[44] = -0.003304F;
	geracao.individuo[0].gene[45] = -0.0032304F;
	geracao.individuo[0].gene[46] = -0.0031599F;
	geracao.individuo[0].gene[47] = -0.0030925F;
	geracao.individuo[0].gene[48] = -0.0030278F;
	geracao.individuo[0].gene[49] = -0.0029659F;

	// Indivíduo 1
	geracao.individuo[1].gene[0] = -0.1389437F;
	geracao.individuo[1].gene[1] = -0.9612598F;
	geracao.individuo[1].gene[2] = 0.1954439F;
	geracao.individuo[1].gene[3] = 0.0887042F;
	geracao.individuo[1].gene[4] = 0.0573714F;
	geracao.individuo[1].gene[5] = 0.042396F;
	geracao.individuo[1].gene[6] = 0.0336202F;
	geracao.individuo[1].gene[7] = 0.0278545F;
	geracao.individuo[1].gene[8] = 0.0237768F;
	geracao.individuo[1].gene[9] = 0.0207406F;
	geracao.individuo[1].gene[10] = 0.018392F;
	geracao.individuo[1].gene[11] = 0.0165212F;
	geracao.individuo[1].gene[12] = 0.0149958F;
	geracao.individuo[1].gene[13] = 0.0137283F;
	geracao.individuo[1].gene[14] = 0.0126584F;
	geracao.individuo[1].gene[15] = 0.0117432F;
	geracao.individuo[1].gene[16] = 0.0109514F;
	geracao.individuo[1].gene[17] = 0.0102596F;
	geracao.individuo[1].gene[18] = 0.0096501F;
	geracao.individuo[1].gene[19] = 0.0091089F;
	geracao.individuo[1].gene[20] = 0.0086251F;
	geracao.individuo[1].gene[21] = 0.0081902F;
	geracao.individuo[1].gene[22] = 0.007797F;
	geracao.individuo[1].gene[23] = 0.0074399F;
	geracao.individuo[1].gene[24] = 0.007114F;
	geracao.individuo[1].gene[25] = 0.0068155F;
	geracao.individuo[1].gene[26] = 0.006541F;
	geracao.individuo[1].gene[27] = 0.0062878F;
	geracao.individuo[1].gene[28] = 0.0060535F;
	geracao.individuo[1].gene[29] = 0.0058359F;
	geracao.individuo[1].gene[30] = 0.0056335F;
	geracao.individuo[1].gene[31] = 0.0054447F;
	geracao.individuo[1].gene[32] = 0.0052681F;
	geracao.individuo[1].gene[33] = 0.0051026F;
	geracao.individuo[1].gene[34] = 0.0049472F;
	geracao.individuo[1].gene[35] = 0.0048009F;
	geracao.individuo[1].gene[36] = 0.0046631F;
	geracao.individuo[1].gene[37] = 0.004533F;
	geracao.individuo[1].gene[38] = 0.0044099F;
	geracao.individuo[1].gene[39] = 0.0042933F;
	geracao.individuo[1].gene[40] = 0.0041827F;
	geracao.individuo[1].gene[41] = 0.0040777F;
	geracao.individuo[1].gene[42] = 0.0039779F;
	geracao.individuo[1].gene[43] = 0.0038828F;
	geracao.individuo[1].gene[44] = 0.0037921F;
	geracao.individuo[1].gene[45] = 0.0037056F;
	geracao.individuo[1].gene[46] = 0.0036229F;
	geracao.individuo[1].gene[47] = 0.0035439F;
	geracao.individuo[1].gene[48] = 0.0034682F;
	geracao.individuo[1].gene[49] = 0.0033957F;

	// Indivíduo 2
	geracao.individuo[2].gene[0] = -0.0795874F;
	geracao.individuo[2].gene[1] = -0.14689F;
	geracao.individuo[2].gene[2] = -0.9516369F;
	geracao.individuo[2].gene[3] = 0.2124869F;
	geracao.individuo[2].gene[4] = 0.0955734F;
	geracao.individuo[2].gene[5] = 0.0616517F;
	geracao.individuo[2].gene[6] = 0.0455018F;
	geracao.individuo[2].gene[7] = 0.0360566F;
	geracao.individuo[2].gene[8] = 0.0298586F;
	geracao.individuo[2].gene[9] = 0.0254789F;
	geracao.individuo[2].gene[10] = 0.0222197F;
	geracao.individuo[2].gene[11] = 0.0196997F;
	geracao.individuo[2].gene[12] = 0.0176931F;
	geracao.individuo[2].gene[13] = 0.0160575F;
	geracao.individuo[2].gene[14] = 0.0146987F;
	geracao.individuo[2].gene[15] = 0.013552F;
	geracao.individuo[2].gene[16] = 0.0125712F;
	geracao.individuo[2].gene[17] = 0.0117228F;
	geracao.individuo[2].gene[18] = 0.0109816F;
	geracao.individuo[2].gene[19] = 0.0103287F;
	geracao.individuo[2].gene[20] = 0.009749F;
	geracao.individuo[2].gene[21] = 0.0092309F;
	geracao.individuo[2].gene[22] = 0.0087651F;
	geracao.individuo[2].gene[23] = 0.008344F;
	geracao.individuo[2].gene[24] = 0.0079616F;
	geracao.individuo[2].gene[25] = 0.0076127F;
	geracao.individuo[2].gene[26] = 0.007293F;
	geracao.individuo[2].gene[27] = 0.0069992F;
	geracao.individuo[2].gene[28] = 0.0067281F;
	geracao.individuo[2].gene[29] = 0.0064772F;
	geracao.individuo[2].gene[30] = 0.0062443F;
	geracao.individuo[2].gene[31] = 0.0060277F;
	geracao.individuo[2].gene[32] = 0.0058255F;
	geracao.individuo[2].gene[33] = 0.0056365F;
	geracao.individuo[2].gene[34] = 0.0054593F;
	geracao.individuo[2].gene[35] = 0.005293F;
	geracao.individuo[2].gene[36] = 0.0051365F;
	geracao.individuo[2].gene[37] = 0.0049889F;
	geracao.individuo[2].gene[38] = 0.0048496F;
	geracao.individuo[2].gene[39] = 0.0047179F;
	geracao.individuo[2].gene[40] = 0.0045932F;
	geracao.individuo[2].gene[41] = 0.0044748F;
	geracao.individuo[2].gene[42] = 0.0043625F;
	geracao.individuo[2].gene[43] = 0.0042556F;
	geracao.individuo[2].gene[44] = 0.0041538F;
	geracao.individuo[2].gene[45] = 0.0040568F;
	geracao.individuo[2].gene[46] = 0.0039642F;
	geracao.individuo[2].gene[47] = 0.0038758F;
	geracao.individuo[2].gene[48] = 0.0037912F;
	geracao.individuo[2].gene[49] = 0.0037102F;

	float *matriz_Identidade;
	matriz_Identidade = (float *)malloc(parametrosGA.numGenes*parametrosGA.numGenes*sizeof(float));
	GeraMatrizIdentidade(matriz_Identidade, parametrosGA.numGenes);
	Fitness_Serial(&geracao, &parametrosGA, H, parametrosGA.numGenes, 0.01F, 0.1F, matriz_Identidade);
	free(matriz_Identidade);

	printf("\n"); printf("Geracao depois = ");

	imprimeGeracao(&geracao, &parametrosGA);

	*/
}

// =========================================================================================================
void SelecaoViaRoleta_serial(	struct generation *host_PopulacaoAntes,
										struct generation *host_PopulacaoDepois,
										struct parametros *host_parametrosGA) {
// =========================================================================================================

	unsigned short int iIndividuo;
	
	float fatorRoleta, vlrRoleta, sumRoleta, iSelecionado;

	for (iIndividuo = 0; iIndividuo < host_parametrosGA->numIndividuos; iIndividuo++) {
		
		//printf("\n"); printf("00 - SelecaoViaRoleta_serial - individuo = %d", iIndividuo);
		
		fatorRoleta = (float)rand() / (float)RAND_MAX;
		//printf("\n"); printf("01 - SelecaoViaRoleta_serial - fatorRoleta = %f", fatorRoleta);
		
		//printf("\n"); printf("02 - SelecaoViaRoleta_serial - host_PopulacaoAntes->sumFitness = %f", host_PopulacaoAntes->sumFitness);

		vlrRoleta =  fatorRoleta * host_PopulacaoAntes->sumFitness;
		//printf("\n"); printf("03 - SelecaoViaRoleta_serial - vlrRoleta = %f", vlrRoleta);
		
		iSelecionado = -1.0F;
		
		sumRoleta = 0;
		//printf("\n"); printf("04 - SelecaoViaRoleta_serial - sumRoleta = %f", sumRoleta);
		do {
			iSelecionado = iSelecionado + 1;
			//printf("\n"); printf("05 - SelecaoViaRoleta_serial - iSelecionado = %d", (unsigned short int)iSelecionado);

			//printf("\n"); printf("06 - SelecaoViaRoleta_serial - individuo[iSelecionado].fitness = %f", host_PopulacaoAntes->individuo[(unsigned short int)iSelecionado].fitness);

			sumRoleta = sumRoleta + host_PopulacaoAntes->individuo[(unsigned short int)iSelecionado].fitness;
			//printf("\n"); printf("07 - SelecaoViaRoleta_serial - sumRoleta = %f", sumRoleta);

		} while (sumRoleta <= vlrRoleta);
		
		//printf("\n"); printf("08 - SelecaoViaRoleta_serial - iSelecionado = %d", (unsigned short int)iSelecionado);
		// Grava o individuo selecionado (iSelecionado) na próxima geração na posicao iIndividuo.
		host_PopulacaoDepois->individuo[iIndividuo] = host_PopulacaoAntes->individuo[(unsigned short int)iSelecionado];
	};
};


//=========================================================================================================
void Selecao_Por_Torneio_serial(	struct generation *PopulacaoAntes,
											struct generation *PopulacaoDepois,
											struct parametros *parametrosGA) {
											//struct testeGeracao_s *h_testeGeracao_s) {
//=========================================================================================================

	unsigned short int iIndividuo;
	char iTamanhoDoTorneio;
	unsigned short int iIndividuo_para_Torneio;
	struct individual *Individuo;
	float melhorFitness;
	
	for (iIndividuo = 0; iIndividuo < parametrosGA->numIndividuos; iIndividuo++) {
		
		melhorFitness = -1.0F;
		//h_testeGeracao_s->individuo[iIndividuo].iIndividuo_Vencedor = 9999; // apenas para teste

		for (iTamanhoDoTorneio = 0; iTamanhoDoTorneio < parametrosGA->tamanho_torneio; iTamanhoDoTorneio++) {
			
			iIndividuo_para_Torneio = Randomico_int(0,parametrosGA->numIndividuos);
			//printf("\n===> SELECAO -- iIndividuo[%d,%d] para o torneio %d", iTamanhoDoTorneio, iIndividuo_para_Torneio, iIndividuo);

			Individuo = &PopulacaoAntes->individuo[iIndividuo_para_Torneio];
			//h_testeGeracao_s->individuo[iIndividuo].teste_iIndividuo_para_Torneio[iTamanhoDoTorneio] = iIndividuo_para_Torneio;
			//printf("\n===> SELECAO -- Individuo.Fitness = [%f]", Individuo->fitness);
			//printf("\n===> SELECAO -- melhorFitness    = [%f]", melhorFitness);

			if (Individuo->fitness > melhorFitness) {
				melhorFitness = Individuo->fitness;
				PopulacaoDepois->individuo[iIndividuo] = *Individuo;
				//h_testeGeracao_s->individuo[iIndividuo].iIndividuo_Vencedor = iIndividuo_para_Torneio;
			};

			//printf("\n===> SELECAO -- melhorFitness    = [%f]", melhorFitness);
		};
	};
};

//=========================================================================================================
void CrossOver1Ponto_serial(	struct generation *g0,
										struct generation *g1,
										struct parametros *parametrosGA ) {
										//struct testeGeracao_s *host_TesteGeracao_s) {
//=========================================================================================================
	unsigned short int iIndividuo, idx_Pai, idx_Mae;
	unsigned long int iGene, ponto_de_corte;

	struct individual *Pai, *Mae;
	float f, p_aux;	// p = probabilidade auxiliar para comparacao com a probabilidade do crossover

	for (iIndividuo = 0; iIndividuo < parametrosGA->numIndividuos; iIndividuo = iIndividuo + 1) {	
		// iIndividuo = iIndividuo + 2 limita numIndividuos como um número par.

		// Coloquei a obtenção do Pai e da Mãe aqui em cima porque eles são utilizados mesmo se 
		// o Crossover não acontece (pra apenas repassar os indivíduos).

		idx_Pai = Randomico_int(0, parametrosGA->numIndividuos - 1);
		//idx_Pai = iIndividuo;
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("idx_Pai =\t%d", idx_Pai);
		
		//host_TesteGeracao_s->individuo[iIndividuo].idx_Pai = idx_Pai;
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("host_TesteGeracao_s->individuo[iIndividuo].idx_Pai =\t%d", host_TesteGeracao_s->individuo[iIndividuo].idx_Pai);
		
		idx_Mae = Randomico_int(0, parametrosGA->numIndividuos - 1);
		//printf("\nCrossOver1Ponto_serial\t"); printf("idx_Mae =\t%d", idx_Mae);
		//host_TesteGeracao_s->individuo[iIndividuo].idx_Mae = idx_Mae;
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("host_TesteGeracao_s->individuo[iIndividuo].idx_Mae =\t%d", host_TesteGeracao_s->individuo[iIndividuo].idx_Mae);

		Pai = &g0->individuo[idx_Pai];
		Mae = &g0->individuo[idx_Mae];

		p_aux = (float)((float)rand() / (float)RAND_MAX);
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("p_aux =\t%f", p_aux);
		//host_TesteGeracao_s->individuo[iIndividuo].p = p_aux;
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("host_TesteGeracao_s->individuo[iIndividuo].p =\t%f", host_TesteGeracao_s->individuo[iIndividuo].p);

		f = -1.0F;
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("f =\t%f", f);
		//host_TesteGeracao_s->individuo[iIndividuo].f = f;
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("host_TesteGeracao_s->individuo[iIndividuo].f =\t%f", host_TesteGeracao_s->individuo[iIndividuo].f);

		ponto_de_corte = Randomico_int(0, (parametrosGA->numGenes-1));
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("ponto_de_corte =\t%d", ponto_de_corte);
		//host_TesteGeracao_s->individuo[iIndividuo].pontos_de_corte[0] = ponto_de_corte;
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("host_TesteGeracao_s->individuo[iIndividuo].pontos_de_corte[0] =\t%d", host_TesteGeracao_s->individuo[iIndividuo].pontos_de_corte[0]);
		//host_TesteGeracao_s->individuo[iIndividuo].pontos_de_corte[1] = 0;
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("host_TesteGeracao_s->individuo[iIndividuo].pontos_de_corte[1] =\t%d", host_TesteGeracao_s->individuo[iIndividuo].pontos_de_corte[1]);
	
		if (p_aux <= parametrosGA->probCrossOver) {
			
			f = (float)((float)rand()/(float)RAND_MAX);
			//host_TesteGeracao_s->individuo[iIndividuo].f = f;			

			for (iGene = 0; iGene < ponto_de_corte; iGene++) {
				g1->individuo[iIndividuo    ].gene[iGene] = Pai->gene[iGene];				
			}
			
			// na segunda parte há alteração, mas apenas entre os dois pontos de corte.
			for (iGene = ponto_de_corte; iGene < parametrosGA->numGenes; iGene++) {				
				g1->individuo[iIndividuo].gene[iGene] =       f*Pai->gene[iGene] + (1 - f)*Mae->gene[iGene];				
			}
		}
		else {
			// Se não acontece o crossover, o mesmo indivíduo da 
			// geração anterior é passado pra próxima.
			for (iGene = 0; iGene < parametrosGA->numGenes; iGene++) {
				g1->individuo[iIndividuo].gene[iGene] =	g0->individuo[iIndividuo].gene[iGene];			
			}
		}
		
	}
}

//=========================================================================================================
void CrossOver2Pontos_serial(	struct generation *g0,
										struct generation *g1,
										struct parametros *parametrosGA) {
										// struct testeGeracao_s *host_TesteGeracao_s) {
//=========================================================================================================
	unsigned short int iIndividuo, idx_Pai, idx_Mae;
	unsigned long int iGene;

	struct individual *Pai, *Mae;
	float f, p_aux;	// p = probabilidade auxiliar para comparacao com a probabilidade do crossover

	for (iIndividuo = 0; iIndividuo < parametrosGA->numIndividuos; iIndividuo = iIndividuo + 1) {	
		// iIndividuo = iIndividuo + 2 limita numIndividuos como um número par.

		// Coloquei a obtenção do Pai e da Mãe aqui em cima porque eles são utilizados mesmo se 
		// o Crossover não acontece (pra apenas repassar os indivíduos).

		idx_Pai = Randomico_int(0, parametrosGA->numIndividuos - 1);
		//idx_Pai = iIndividuo;
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("idx_Pai =\t%d", idx_Pai);
		//host_TesteGeracao_s->individuo[iIndividuo].idx_Pai = idx_Pai;
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("host_TesteGeracao_s->individuo[iIndividuo].idx_Pai =\t%d", host_TesteGeracao_s->individuo[iIndividuo].idx_Pai);
		
		idx_Mae = Randomico_int(0, parametrosGA->numIndividuos - 1);
		//printf("\nCrossOver1Ponto_serial\t"); printf("idx_Mae =\t%d", idx_Mae);
		//host_TesteGeracao_s->individuo[iIndividuo].idx_Mae = idx_Mae;
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("host_TesteGeracao_s->individuo[iIndividuo].idx_Mae =\t%d", host_TesteGeracao_s->individuo[iIndividuo].idx_Mae);

		Pai = &g0->individuo[idx_Pai];
		Mae = &g0->individuo[idx_Mae];

		p_aux = (float)((float)rand() / (float)RAND_MAX);
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("p_aux =\t%f", p_aux);
		//host_TesteGeracao_s->individuo[iIndividuo].p = p_aux;
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("host_TesteGeracao_s->individuo[iIndividuo].p =\t%f", host_TesteGeracao_s->individuo[iIndividuo].p);

		f = -1.0F;
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("f =\t%f", f);
		//host_TesteGeracao_s->individuo[iIndividuo].f = f;
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("host_TesteGeracao_s->individuo[iIndividuo].f =\t%f", host_TesteGeracao_s->individuo[iIndividuo].f);

		// Com a utilização de 'g1->individuo[iIndividuo].pontos_de_corte[0]'
		// no segundo Randomico eu garanto que ponto1 >= ponto2.
		// Assim, não preciso de código para ordenação.
		g1->individuo[iIndividuo].pontos_de_corte[0] = Randomico_int(0, (parametrosGA->numGenes-1));
		g1->individuo[iIndividuo].pontos_de_corte[1] = Randomico_int(g1->individuo[iIndividuo].pontos_de_corte[0], (parametrosGA->numGenes-1));

		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("ponto_de_corte =\t%d", ponto_de_corte);
		//host_TesteGeracao_s->individuo[iIndividuo].pontos_de_corte[0] = g1->individuo[iIndividuo].pontos_de_corte[0];
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("host_TesteGeracao_s->individuo[iIndividuo].pontos_de_corte[0] =\t%d", host_TesteGeracao_s->individuo[iIndividuo].pontos_de_corte[0]);
		//host_TesteGeracao_s->individuo[iIndividuo].pontos_de_corte[1] = g1->individuo[iIndividuo].pontos_de_corte[1];
		// TESTE - printf("\nCrossOver1Ponto_serial\t"); printf("host_TesteGeracao_s->individuo[iIndividuo].pontos_de_corte[1] =\t%d", host_TesteGeracao_s->individuo[iIndividuo].pontos_de_corte[1]);
	
		if (p_aux <= parametrosGA->probCrossOver) {
			
			f = (float)((float)rand()/(float)RAND_MAX);
			//host_TesteGeracao_s->individuo[iIndividuo].f = f;			

			for (iGene = 0; iGene < g1->individuo[iIndividuo].pontos_de_corte[0]; iGene++) {
				g1->individuo[iIndividuo    ].gene[iGene] = Pai->gene[iGene];				
			}
			
			// na segunda parte há alteração, mas apenas entre os dois pontos de corte.
			for (iGene = g1->individuo[iIndividuo].pontos_de_corte[0]; iGene < g1->individuo[iIndividuo].pontos_de_corte[1]; iGene++) {				
				g1->individuo[iIndividuo].gene[iGene] =       f*Pai->gene[iGene] + (1 - f)*Mae->gene[iGene];				
			}

			for (iGene = g1->individuo[iIndividuo].pontos_de_corte[1]; iGene < parametrosGA->numGenes; iGene++) {
				g1->individuo[iIndividuo    ].gene[iGene] = Pai->gene[iGene];				
			}
		}
		else {
			// Se não acontece o crossover, o mesmo indivíduo da 
			// geração anterior é passado pra próxima.
			for (iGene = 0; iGene < parametrosGA->numGenes; iGene++) {
				g1->individuo[iIndividuo].gene[iGene] =	g0->individuo[iIndividuo].gene[iGene];			
			}
		}
	}
}

//---------------------------------------------------------------------------------------
void Mutacao_Serial(	struct generation *geracao,
							struct parametros *parametrosGA) {
							//struct testeGeracao_s *host_testeGeracao) {
//---------------------------------------------------------------------------------------

	/*
	* A mutação seguiu a equação 10 do artigo "Diagonalization of a real-symmetric
	* Hamiltonian by genetic algorithm".
	*/

	unsigned short int iIndividuo;
	unsigned long iGene;
	unsigned short int L;
	int termo_do_L;
	float r, pAux;
	float intensidadeMutacao = 10*geracao->Maior_fitness;

	for (iIndividuo = 0; iIndividuo < parametrosGA->numIndividuos; iIndividuo++) {
		
		for (iGene = 0; iGene < parametrosGA->numGenes; iGene++) {
	
			pAux = (float)((float)rand() / (float)RAND_MAX);
			//host_testeGeracao->individuo[iIndividuo].gene_pAux[iGene] = pAux;
	
			if (pAux <= parametrosGA->probMutacao) {
			//	printf("# OCORRERA MUTACAO NESSE GENE#\n");	
				L = Randomico_int(0, 10);
				//host_testeGeracao->individuo[iIndividuo].gene_L[iGene] = L;
					// No artigo esse L é definido como um número
					// randômico inteiro. Nada é dito sobre os limites
					// inferior ou superior de L. Decidi arbitrariamente
					// que L = (0, 10].
				r = (float)((float)rand() / (float)RAND_MAX) + 0.000002;
				//host_testeGeracao->individuo[iIndividuo].gene_r[iGene] = r;
					// r é um número randômico entre zero e 1. Como
					// rand() = (0, RAND_MAX), a expressão acima funciona
					// quase perfeitamente. Quase pois nunca retornará
					// o zero.
				termo_do_L = (int)pow((float)(-1),(int)L);
				//host_testeGeracao->individuo[iIndividuo].gene_termo_do_L[iGene] = termo_do_L;
					// é o (-1) elevado ao L. Só pode ser + ou - 1,
					// portanto, um inteiro.

				geracao->individuo[iIndividuo].gene[iGene] =
						geracao->individuo[iIndividuo].gene[iGene] +
						(float)(termo_do_L * r * intensidadeMutacao);
						//(float)(termo_do_L * r * parametrosGA->intensidadeMutacao);
								// parametrosGA->intensidadeMutacao é o 
								// DELTA da equação 10 do artigo.
			}
		}
	}
}

//---------------------------------------------------------------------------------------
void GeraPopulacaoInicial_serial(
						struct generation *GeracaoInicial,
						struct parametros *parametrosGA ) {
						//struct testeGeracao_s *host_testeGeracao) {
//---------------------------------------------------------------------------------------

	unsigned short int iIndividuo, qtdePontosCorte, iPontoCorte;
	unsigned long int iGene;

	GeracaoInicial->FitnessMedio = -1.0F;
	GeracaoInicial->idxMelhorIndividuo = 99;
	GeracaoInicial->Maior_fitness = -1.0F;
	GeracaoInicial->Melhor_cociente_Rayleigh = -1.0F;
	GeracaoInicial->sumFitness = -1.0F;	
	
	unsigned short int L;
	short int sinal;
	
	for (iIndividuo = 0; iIndividuo < parametrosGA->numIndividuos; iIndividuo++) {		
		
		GeracaoInicial->individuo[iIndividuo].cociente_Rayleigh = -1.0F;
		GeracaoInicial->individuo[iIndividuo].CtC = -1.0F;
		GeracaoInicial->individuo[iIndividuo].fitness = -1.0F;
		GeracaoInicial->individuo[iIndividuo].grad_elevado_ao_quadrado = -1.0F;
		GeracaoInicial->individuo[iIndividuo].inverso_de_CtC = -1.0F;
		GeracaoInicial->individuo[iIndividuo].Numerador = -1.0F;
		GeracaoInicial->individuo[iIndividuo].rho_menos_rho0_ao_quadrado = -1.0F;
		
		qtdePontosCorte = parametrosGA->qtde_Pontos_de_Corte;
		for (iPontoCorte = 0; iPontoCorte < qtdePontosCorte; iPontoCorte++) {
			GeracaoInicial->individuo[iIndividuo].pontos_de_corte[iPontoCorte] = 11111;
		}

		for (iGene = 0; iGene < parametrosGA->numGenes; iGene++) {			
		
			//host_testeGeracao->individuo[iIndividuo].gene_L[iGene] = 11111;
			//host_testeGeracao->individuo[iIndividuo].gene_pAux[iGene] = 1.0F;
			//host_testeGeracao->individuo[iIndividuo].gene_r[iGene] = 1.0F;
			//host_testeGeracao->individuo[iIndividuo].gene_termo_do_L[iGene] = 11111;
		
			L = Randomico_int(0, 2);
			sinal = (short int)pow(-1.0F, (float)L);
			GeracaoInicial->individuo[iIndividuo].gene[iGene] = (float) sinal * rand() / RAND_MAX;
		}
	}
}


//---------------------------------------------------------------------------------------
void GeraPopInicial_vetsOrtonormais(
						struct generation *GeracaoInicial,
						struct parametros *parametrosGA) {
						//struct testeGeracao_s *host_testeGeracao) {
//---------------------------------------------------------------------------------------

	unsigned short int iIndividuo, qtdePontosCorte, iPontoCorte;
	unsigned long int iGene;

	GeracaoInicial->FitnessMedio = -1.0F;
	GeracaoInicial->idxMelhorIndividuo = 99;
	GeracaoInicial->Maior_fitness = -1.0F;
	GeracaoInicial->Melhor_cociente_Rayleigh = -1.0F;
	GeracaoInicial->sumFitness = -1.0F;	
	
	unsigned short int L;
	short int sinal;
	
	// Parâmetros específicos dos indivíduos.
	for (iIndividuo = 0; iIndividuo < parametrosGA->numIndividuos; iIndividuo++) {		
		
		GeracaoInicial->individuo[iIndividuo].cociente_Rayleigh = -1.0F;
		GeracaoInicial->individuo[iIndividuo].CtC = -1.0F;
		GeracaoInicial->individuo[iIndividuo].fitness = -1.0F;
		GeracaoInicial->individuo[iIndividuo].grad_elevado_ao_quadrado = -1.0F;
		GeracaoInicial->individuo[iIndividuo].inverso_de_CtC = -1.0F;
		GeracaoInicial->individuo[iIndividuo].Numerador = -1.0F;
		GeracaoInicial->individuo[iIndividuo].rho_menos_rho0_ao_quadrado = -1.0F;
		
		qtdePontosCorte = parametrosGA->qtde_Pontos_de_Corte;
		for (iPontoCorte = 0; iPontoCorte < qtdePontosCorte; iPontoCorte++) {
			GeracaoInicial->individuo[iIndividuo].pontos_de_corte[iPontoCorte] = 11111;
		}

		float flutuacao = 0;

		for (iGene = 0; iGene < parametrosGA->numIndividuos; iGene++) {			
			
			L = Randomico_int(0, 2);
			sinal = (short int)pow(-1.0F, (float)L);

			//flutuacao = (float) (sinal * rand()/RAND_MAX)*(float)0.01F;
			flutuacao = (float) sinal * 0.00001F;

			if (iGene == iIndividuo)
				GeracaoInicial->individuo[iIndividuo].gene[iGene] = 1.0F;
			else
				GeracaoInicial->individuo[iIndividuo].gene[iGene] = flutuacao;
	
		}

	}

	
	
}

//---------------------------------------------------------------------------------------
void teste_GeraPopInicial_vetsOrtonormais(void) {
//---------------------------------------------------------------------------------------

	struct generation geracaoInicial;
	struct parametros parametrosGA;

	parametrosGA.numIndividuos = 10;
	parametrosGA.numGenes = 10;

	GeraPopInicial_vetsOrtonormais(
						&geracaoInicial,
						&parametrosGA);

	imprimeGeracao(&geracaoInicial, &parametrosGA);

}
//---------------------------------------------------------------------------------------
void teste_GeraPopulacaoInicial_serial(void) {
//---------------------------------------------------------------------------------------
	//TODO: deve ser alterada para inicialização manua dos parâmetros

	struct generation				geracaoTeste;
	struct parametros				parametrosGA;
	//struct parametros_Metodo	parametrosMetodo;
	//struct testeGeracao_s		host_testeGeracao_s;
	
	//inicializa_Parametros(&parametrosGA, &parametrosMetodo);
	GeraPopulacaoInicial_serial(&geracaoTeste, &parametrosGA); // , &host_testeGeracao_s);

	imprimeGeracao(&geracaoTeste, &parametrosGA);

}

//===========================================================================================
void testeSelecao_serial(	unsigned int flag_Momento,
									struct generation *host_Geracao,
									struct parametros *host_ParametrosGA,
									struct testeGeracao_s *host_testeSelecao) {
//===========================================================================================

	if (flag_Momento == 0) {
		// antes da seleção atuar
		// apenas imprimo a geracao.
		printf("\n");	printf("Teste Selecao - Populacao Antes");
		printf("\n");	printf("-----------------------------------------------------------");
		imprimeGeracao(host_Geracao, host_ParametrosGA);
	}
	else if (flag_Momento == 1) {
		unsigned int iIndividuo;
		int iTorneio;
		printf("\n");	printf("Teste Selecao - Populacao Depois");
		printf("\n");	printf("-----------------------------------------------------------");
		printf("\n");	printf("#");
		printf("\t");	printf("Individuo 0");
		printf("\t");	printf("Individuo 1");
		printf("\t");	printf("Vencedor");
		for (iIndividuo = 0; iIndividuo < host_ParametrosGA->numIndividuos; iIndividuo++) {
			printf("\n");	printf("%d", iIndividuo);
			for (iTorneio = 0; iTorneio < host_ParametrosGA->tamanho_torneio; iTorneio++) {
				printf("\t");	printf("%d", host_testeSelecao->individuo[iIndividuo].teste_iIndividuo_para_Torneio[iTorneio]);
			}
			printf("\t");	printf("%d", host_testeSelecao->individuo[iIndividuo].iIndividuo_Vencedor);
		}
	}
}


//===========================================================================================
void testeCrossOver_serial(	unsigned int flag_Momento,
										struct generation *host_Geracao,
										struct parametros *host_ParametrosGA,
										struct testeGeracao_s *host_testeGeracao) {
//===========================================================================================
	
	if (flag_Momento == 2) { // antes da aplicação da CrossOver
		printf("\n");	printf("Teste Crossover - Proxima Populacao Antes de receber os novos genes");
		printf("\n");	printf("-----------------------------------------------------------");

		imprimeGeracao(host_Geracao, host_ParametrosGA);
	}

	if (flag_Momento == 0) { // antes da aplicação da CrossOver
		printf("\n");	printf("Teste Crossover - Populacao Antes");
		printf("\n");	printf("-----------------------------------------------------------");

		imprimeGeracao(host_Geracao, host_ParametrosGA);
	}
	else if (flag_Momento == 1) { // depois da aplicação da CrossOver
		
		printf("\n");	printf("Teste Crossover - Termos auxiliares");
		printf("\n");	printf("-----------------------------------------------------------");

		unsigned int iIndividuo;

		printf("\n");	printf("iIndividuo");
		printf("\t");	printf("idx_Pai");
		printf("\t");	printf("idx_Mae");
		printf("\t");	printf("p");
		printf("\t");	printf("f");
		printf("\t");	printf("pontos_de_corte[0]");
		printf("\t");	printf("pontos_de_corte[1]");

		for (iIndividuo = 0; iIndividuo < host_ParametrosGA->numIndividuos; iIndividuo++) {
			printf("\n");	printf("%d", iIndividuo);
			printf("\t");	printf("%d", host_testeGeracao->individuo[iIndividuo].idx_Pai);
			printf("\t");	printf("%d", host_testeGeracao->individuo[iIndividuo].idx_Mae);
			printf("\t");	printf("%f", host_testeGeracao->individuo[iIndividuo].p);
			printf("\t");	printf("%f", host_testeGeracao->individuo[iIndividuo].f);
			printf("\t");	printf("%d", host_testeGeracao->individuo[iIndividuo].pontos_de_corte[0]);
			printf("\t");	printf("%d", host_testeGeracao->individuo[iIndividuo].pontos_de_corte[1]);
		}

		printf("\n");	printf("Teste Crossover - Populacao Depois");
		printf("\n");	printf("-----------------------------------------------------------");
		imprimeGeracao(host_Geracao, host_ParametrosGA);
	}
}

//===========================================================================================
void testeMutacao(int flag_Momento,
						struct generation *host_Geracao,
						struct parametros *host_parametros,
						struct testeGeracao_s *host_testeGeracao) {
//===========================================================================================
	
	unsigned int iIndividuo, iGene;
	
	//struct generation host_Geracao;
	//cudaErrorCheck(cudaMemcpy(&host_Geracao, dev_Geracao, sizeof(struct generation), cudaMemcpyDeviceToHost));

	if (flag_Momento == 0) { // antes da Mutacao
		printf("\n");	printf("===========================================================");
		printf("\n");	printf("TESTE MUTACAO: antes da devMutacao");
		printf("\n");	printf("===========================================================");
		for (iIndividuo = 0; iIndividuo < host_parametros->numIndividuos; iIndividuo++) {
			
			printf("\n");	printf("Individuo %d", iIndividuo);
			printf("\n");	printf("-----------------------------------------------------------");
			printf("\n");	printf("iGene");
			printf("\t");	printf("Gene");
			//printf("\t");	printf("pAux");
			//printf("\t");	printf("L");
			//printf("\t");	printf("r");
			//printf("\t");	printf("Termo do L");

			for (iGene = 0; iGene < host_parametros->numGenes; iGene++) {
				printf("\n");	printf("%d", iGene);
				printf("\t");	printf("%f", host_Geracao->individuo[iIndividuo].gene[iGene]);
				//printf("\t");	printf("%f", host_Geracao->individuo[iIndividuo].gene_pAux[iGene]);
				//printf("\t");	printf("%f", host_Geracao->individuo[iIndividuo].gene_L[iGene]);
				//printf("\t");	printf("%f", host_Geracao->individuo[iIndividuo].gene_r[iGene]);
				//printf("\t");	printf("%f", host_Geracao->individuo[iIndividuo].gene_termo_do_L[iGene]);
			}
		}
	} 
	else if (flag_Momento == 1) { // depois da mutacao
		printf("\n");	printf("===========================================================");
		printf("\n");	printf("TESTE MUTACAO: depois da devMutacao");
		printf("\n");	printf("===========================================================");
		for (iIndividuo = 0; iIndividuo < host_parametros->numIndividuos; iIndividuo++) {
			
			printf("\n");	printf("Individuo %d", iIndividuo);
			printf("\n");	printf("-----------------------------------------------------------");
			printf("\n");	printf("iGene");
			printf("\t");	printf("Gene");
			printf("\t");	printf("pAux");
			printf("\t");	printf("L");
			printf("\t");	printf("r");
			printf("\t");	printf("Termo do L");

			for (iGene = 0; iGene < host_parametros->numGenes; iGene++) {
				printf("\n");	printf("%d", iGene);
				printf("\t");	printf("%f", host_Geracao->individuo[iIndividuo].gene[iGene]);
				printf("\t");	printf("%f", host_testeGeracao->individuo[iIndividuo].gene_pAux[iGene]);
				printf("\t");	printf("%d", host_testeGeracao->individuo[iIndividuo].gene_L[iGene]);
				printf("\t");	printf("%f", host_testeGeracao->individuo[iIndividuo].gene_r[iGene]);
				printf("\t");	printf("%d", host_testeGeracao->individuo[iIndividuo].gene_termo_do_L[iGene]);
			}
		}
	}
}

//---------------------------------------------------------------------------------------
unsigned short int imprimeComportamentoFitness(
							unsigned short int flag_Cabecalho,
							unsigned short int cod_Maquina,
							unsigned short int cod_tipoPrograma,
							unsigned short int cod_Tipo_Fitness,
							unsigned long	int parm_iGeracao,
							struct generation *host_Geracao,
							struct parametros_Metodo *host_parametros_Metodo,
							struct parametros *host_parametrosGA,
							struct parametrosPrograma *host_parametrosPrograma) {
//---------------------------------------------------------------------------------------

	if (flag_Cabecalho == 0) {
		printf("\n"); printf("Comportamento do Fitness");
		printf("\t"); printf("Maquina");						// 0 - Adriano, 1 - FT
		printf("\t"); printf("Serial ou Paralelo?");		// 0 - Serial, 1 - Paralelo
		printf("\t"); printf("Quantidade Geracoes");		// constNumGeracoes
		printf("\t"); printf("Quantidade Individuos");	// constNumIndividuos
		printf("\t"); printf("Quantidade Genes");			// constNumGenes
		printf("\t"); printf("Tipo Fitness");
		printf("\t"); printf("Lambda");
		printf("\t"); printf("Rho minimo");
		printf("\t"); printf("Rho medio");
		printf("\t"); printf("Diferenca Rho");
		printf("\t"); printf("iGeracao");				// parm_iGeracao
		printf("\t"); printf("Fitness Medio");
		printf("\t"); printf("Maior Fitness");
		printf("\t"); printf("Melhor rho");
		printf("\t"); printf("Grad rho medio");
		printf("\t"); printf("Posicao Melhor Individuo");
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

		printf("\n"); printf("Comportamento do Fitness");
		printf("\t"); printf("%d", cod_Maquina);			
		printf("\t"); printf("%s", str_TipoPrograma);	
		printf("\t"); printf("%d", host_parametrosPrograma->parmQtdeMaxGeracoes);	
		printf("\t"); printf("%d", host_parametrosGA->numIndividuos);	
		printf("\t"); printf("%d", host_parametrosGA->numGenes);		
		printf("\t"); printf("%d", cod_Tipo_Fitness);
		printf("\t"); printf("%f", host_parametros_Metodo->lambda);
		printf("\t"); printf("%f", host_parametros_Metodo->rho_minimo);
		printf("\t"); printf("%f", host_Geracao->rhoMedio);
		printf("\t"); printf("%f", host_Geracao->difRho);
		printf("\t"); printf("%d", parm_iGeracao);
		printf("\t"); printf("%f", host_Geracao->FitnessMedio);
		printf("\t"); printf("%f", host_Geracao->Maior_fitness);
		printf("\t"); printf("%f", host_Geracao->Melhor_cociente_Rayleigh);
		printf("\t"); printf("%f", host_Geracao->gradRhoMedio);
		printf("\t"); printf("%d", host_Geracao->idxMelhorIndividuo);

	} // final do if flag_Cabealho;

	return 0;
}

//------------------------------------------------------------------------------------
unsigned short int atingiuCriterioDeParada(
							unsigned short int tipoFitness,
							float tolerancia,
							struct generation *geracao) {
//------------------------------------------------------------------------------------

	unsigned short int flagAtingiuTolerancia = 0; // falso

	if (tipoFitness == 0 || tipoFitness == 2) {
		if ( (geracao->gradRhoMedio <= tolerancia) || (geracao->difRho <= tolerancia) ) {
			flagAtingiuTolerancia = 1; // verdadeiro
		}
	}
	else {
		if (tipoFitness == 1) {
			if (geracao->gradRhoMedio <= tolerancia) {	
				flagAtingiuTolerancia = 1; // verdadeiro
			}
		}
	}

	return flagAtingiuTolerancia;
}