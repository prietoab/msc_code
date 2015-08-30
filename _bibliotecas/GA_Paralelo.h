// Cálculo do número de threads necessárias
/*
Maximo de threads por bloco ........: 	 512
Dimensoes	maximas para um bloco ....: 	 (512, 512, 64) 
Dimensoes maximas para o grid ......: 	 (65535, 65535, 1) 
*/

/*
----------------------------------------------------------
Mapeamento SciLab
----------------------------------------------------------
Device 0: "GeForce GT 130M"
  CUDA Driver Version:                           6.50
  CUDA Runtime Version:                          4.0
  CUDA Capability Major revision number:         1
  CUDA Capability Minor revision number:         1
  Total amount of global memory:                 536870912 bytes
  Number of multiprocessors:                     4
  Total amount of constant memory:               65536 bytes
  Total amount of shared memory per block:       16384 bytes
  Total number of registers available per block: 8192
  Warp size:                                     32
  Maximum number of threads per block:           512
  Maximum sizes of each dimension of a block:    512 x 512 x 64
  Maximum sizes of each dimension of a grid:     65535 x 65535 x 1
  Maximum memory pitch:                          2147483647 bytes
  Texture alignment:                             256 bytes
  Clock rate:                                    1.50 GHz
  Concurrent copy and execution:                 Yes
  Run time limit on kernels:                     Yes
  Integrated:                                    No
  Support host page-locked memory mapping:       Yes
  Compute mode:                                  Default (multiple host threads can use this device simultaneously)
  Concurrent kernel execution:                   No
  Device has ECC support enabled:                No
*/

//-------------------------------------------------------------------------------------------
__global__ void dev_GradRho_CalculaInversoCtC(struct generation *dev_Geracao,
															 unsigned int host_Individuo) {
//-------------------------------------------------------------------------------------------

	// Essa função global executa, na verdade, de maneira serial.
	// Motivo: perda de informação da matriz auxiliar Ox1.

	unsigned int iThread = blockIdx.x * blockDim.x + threadIdx.x;

	if (iThread == 0) {
		dev_Geracao->individuo[host_Individuo].inverso_de_CtC = (float) 1.0F / dev_Geracao->individuo[host_Individuo].CtC;
	}
}


//-------------------------------------------------------------------------------------------
__global__ void dev_CalculaFitness(	struct generation *devGeracao,
												struct parametros *devParametrosGA,
												struct parametros_Metodo *dev_ParametrosMetodo) {
//-------------------------------------------------------------------------------------------
	unsigned int iIndividuo = threadIdx.x + blockIdx.x*blockDim.x;

	if (iIndividuo < devParametrosGA->numIndividuos) {
		float rho_menos_rho0 = devGeracao->individuo[iIndividuo].cociente_Rayleigh - dev_ParametrosMetodo->rho_minimo;
		devGeracao->individuo[iIndividuo].rho_menos_rho0_ao_quadrado = pow( rho_menos_rho0, 2);
		devGeracao->individuo[iIndividuo].fitness = exp( (-1)*dev_ParametrosMetodo->lambda*(devGeracao->individuo[iIndividuo].rho_menos_rho0_ao_quadrado + devGeracao->individuo[iIndividuo].grad_elevado_ao_quadrado));
	}
}

/*
//===========================================================================================
__global__ void dev_CalculaRho(struct generation *dev_Geracao,
										 struct parametros *dev_ParametrosGA) {
//===========================================================================================
	
	unsigned int iIndividuo = threadIdx.x + blockIdx.x*blockDim.x;
	
	if (iIndividuo < dev_ParametrosGA->numIndividuos) {
		dev_Geracao->individuo[iIndividuo].cociente_Rayleigh = (dev_Geracao->individuo[iIndividuo].Numerador / dev_Geracao->individuo[iIndividuo].CtC);
	}
}
*/

//===========================================================================================
__global__ void dev_CalculaRho_Serial_naGPU(	unsigned int host_iIndividuo,
															struct generation *dev_Geracao) {
//===========================================================================================
	
	unsigned int iThread = threadIdx.x + blockIdx.x*blockDim.x;
	
	if (iThread == 0) {
		dev_Geracao->individuo[host_iIndividuo].cociente_Rayleigh = (dev_Geracao->individuo[host_iIndividuo].Numerador / dev_Geracao->individuo[host_iIndividuo].CtC);
	}
}

//===========================================================================================
void dev_CalculaRho(
				unsigned int host_iIndividuo,
				float *devHamiltoniano,
				struct generation *devGeracao,
				struct parametros *host_ParametrosGA,
				struct parametros *dev_ParametrosGA,
				struct parametros_Metodo *dev_ParametrosMetodo,
				float *dev_Matriz_Auxiliar_OxO,
				float *dev_Matriz_Auxiliar_Ox1,
				float *dev_Matriz_Identidade) {
	
//===========================================================================================

	/*-----------------------------------------------
	** Cálculo do cociente de Rayleigh - Início
	**-----------------------------------------------*/

	//float Numerador = 0,
	//		CtC = 0;

	//rho = [Ct * H * C  ]  / [ Ct * C ] 
	//
	//	onde
	//		C = vetor C (candidato à autovetor).
	//		Ct = transposto do vetor C.
	//		H = hamiltoniano.g

	float *host_Matriz_Auxiliar_OxO, *host_Matriz_Auxiliar_Ox1;

	host_Matriz_Auxiliar_Ox1 = (float *)malloc(host_ParametrosGA->numGenes*sizeof(float));
	host_Matriz_Auxiliar_OxO = (float *)malloc(host_ParametrosGA->numGenes*host_ParametrosGA->numGenes*sizeof(float));

// Numerador [Ct * H * C]

	printf("\n"); printf("Individuo %d\tPASSO 2\tdev_CalculaRho\tCalculo Ct * H", host_iIndividuo);
	printf("\n"); printf("Genes do Ct (Ox1 antes ):");
	cudaMemcpy(&host_Matriz_Auxiliar_Ox1, &devGeracao->individuo[host_iIndividuo].gene[0], host_ParametrosGA->numGenes*sizeof(float), cudaMemcpyDeviceToHost);
	imprimeMatriz(&host_Matriz_Auxiliar_Ox1, 1, host_ParametrosGA->numGenes);

	cudaErrorCheck( cudaMemcpy(&host_Geracao, devGeracao0, sizeof(struct generation), cudaMemcpyDeviceToHost) );
	imprimeGeracao(&host_Geracao, hostParametrosGA);

	//Ct * H                       devGeracao->iIndividuo
	MatMul(	&devGeracao->individuo[host_iIndividuo].gene[0], // Ct
				1,  // número de linhas de Ct
				host_ParametrosGA->numGenes, // número de colunas de Ct (ordem do hamiltoniano)
				devHamiltoniano, // Hamiltoniano H
				host_ParametrosGA->numGenes, // número de linhas de H
				host_ParametrosGA->numGenes, // número de colunas de H
				dev_Matriz_Auxiliar_Ox1, // resultado da multiplicação
				host_ParametrosGA
	);

	printf("\n"); printf("Individuo %d\tPASSO 3\tdev_CalculaRho\tCalculo Ct", host_iIndividuo);
	printf("\n"); printf("Genes do Ct:");
	cudaMemcpy(&host_Matriz_Auxiliar_Ox1, &devGeracao->individuo[host_iIndividuo].gene[0], host_ParametrosGA->numGenes*sizeof(float), cudaMemcpyDeviceToHost);
	imprimeMatriz(&host_Matriz_Auxiliar_Ox1, 1, host_ParametrosGA->numGenes);

	// Pelos testes que fiz nesse passo, os dois comandos abaixo são desnecessários.
	//cudaErrorCheck( cudaDeviceSynchronize() );
	//cudaErrorCheck( cudaPeekAtLastError() );

	// == TESTE PASSO 1 - inicio

	/*
	printf("\n"); printf("==================================================");
	printf("\n"); printf("Passo 1: Ct * H");
	printf("\n"); printf("==================================================");

	float *host_Matriz_Auxiliar_Ox1;
	host_Matriz_Auxiliar_Ox1 = (float *)malloc(host_ParametrosGA->numGenes*sizeof(float));
	
	cudaErrorCheck(cudaMemcpy(host_Matriz_Auxiliar_Ox1, dev_Matriz_Auxiliar_Ox1, host_ParametrosGA->numGenes*sizeof(float), cudaMemcpyDeviceToHost));
	printf("\n"); printf("host_Matriz_Auxiliar_Ox1 = [");
	printf("\n"); imprimeMatriz(host_Matriz_Auxiliar_Ox1, 1, host_ParametrosGA->numGenes);
	printf("\n"); printf("] ");
	*/
	//host_Matriz_Auxiliar_Ox1 = NULL; free(host_Matriz_Auxiliar_Ox1);
	
	// == TESTE PASSO 1 - FIM

	
	// (Ct * H) * C

	
	MatMul(	dev_Matriz_Auxiliar_Ox1, // Ct * H
				1, // número de linhas de Ct * H
				host_ParametrosGA->numGenes, // número de colunas de Ct * H (ordem do hamiltoniano)
				&devGeracao->individuo[host_iIndividuo].gene[0], // vetor C
				host_ParametrosGA->numGenes, // número de linhas de C
				1, // número de colunas de C
				&devGeracao->individuo[host_iIndividuo].Numerador,
				host_ParametrosGA
	);
	
	/*
	printf("\n"); printf("==================================================");
	printf("\n"); printf("Passo 2: (Ct * H)*C");
	printf("\n"); printf("==================================================");

	float host_Numerador;
	cudaErrorCheck( cudaMemcpy(&host_Numerador, &devGeracao->individuo[host_iIndividuo].Numerador, sizeof(float), cudaMemcpyDeviceToHost) );
	printf("\n"); printf("(Ct * H)*C [Numerador] = %f", host_Numerador);
	*/

	/*
	printf("\n"); printf("==================================================");
	printf("\n"); printf("Passo 3: Denominador [ Ct * C ]");
	printf("\n"); printf("==================================================");
	*/

	// Denominador [ Ct * C ]
	/*
		**** Possível ponto de melhoria *****
		Como o acesso à memória dos elementos de um vetor é muito mais
		rápido quando eles são acessados contíguamente, apenas definir
		C transposto como uma mudança de índices i->j e j->i pode levar
		a um tempo excessivo. Talvez haja ganho considerável se, utilizando
		cuBLAS:

		1) Transpor de fato a matriz C antes de usá-la. Contra: mais uso
			de memória global
		2) Configurar a função de multiplicação do BLAS para não transpor
			a matriz (questão do acesso do FORTRAN). Não tenho certeza se
			dá certo.
	*/
	
	
	MatMul(	&devGeracao->individuo[host_iIndividuo].gene[0],
				1, 
				host_ParametrosGA->numGenes,
				&devGeracao->individuo[host_iIndividuo].gene[0],
				host_ParametrosGA->numGenes,
				1,
				&devGeracao->individuo[host_iIndividuo].CtC,
				host_ParametrosGA
	);

	/*
	float host_Denominador;
	cudaErrorCheck( cudaMemcpy(&host_Denominador, &devGeracao->individuo[host_iIndividuo].CtC, sizeof(float), cudaMemcpyDeviceToHost) );
	printf("\n"); printf("Denominador [ Ct * C ] = %f", host_Denominador);
	*/
	
	dev_CalculaRho_Serial_naGPU<<<1,32>>>(host_iIndividuo, devGeracao);

	free(host_Matriz_Auxiliar_Ox1);
	frre(host_Matriz_Auxiliar_OxO);

}

//===========================================================================================
void GradRho(	unsigned int host_iIndividuo,
					float *devHamiltoniano,
					struct generation *devGeracao,
					struct parametros *host_ParametrosGA,
					struct parametros *dev_ParametrosGA,
					struct parametros_Metodo *dev_ParametrosMetodo,
					float *dev_Matriz_Auxiliar_OxO,
					float *dev_Matriz_Auxiliar_Ox1,
					float *dev_Matriz_Identidade) {
//===========================================================================================

	unsigned int numLinhas, numColunas;
	//
	// IMPORTANTE: essa função é executado para cada indivíduo.
	//
	//-----------------------------------------------
	// Cálculo do Gradiente de Rho - Início
	//-----------------------------------------------
	// 2*[H - rho*I]*C / Ct*C


	//printf("\n"); printf("Inicio do calculo do Gradiente de Rho");
	
	// rho*I -----------------------------------------------------------------
	
	/*
	printf("\n"); printf("==================================================");
	printf("\n"); printf("Passo 5: rho*I");
	printf("\n"); printf("==================================================");	
	*/
	
	numLinhas = host_ParametrosGA->numGenes;
	numColunas = host_ParametrosGA->numGenes;

	multiplica_matriz_por_escalar(
						dev_Matriz_Identidade,
						numLinhas,
						numColunas,
						1, // o escalar está na GPU
						&devGeracao->individuo[host_iIndividuo].cociente_Rayleigh,
						-99, // coloquei -99 como uma valor arbitrário, já que utilizarei um escalar na gpu
						dev_Matriz_Auxiliar_OxO);

	// H - rho*I ----------------------------------------------------------------------

	/*
	printf("\n"); printf("==================================================");
	printf("\n"); printf("Passo 6: H - rho*I");
	printf("\n"); printf("==================================================");
	*/

	unsigned int numBlocosLinha = (numLinhas + nThreadsPorBloco)/nThreadsPorBloco;
	unsigned int numBlocosColuna = (numColunas + nThreadsPorBloco)/nThreadsPorBloco;

	dim3 numThreadsI(nThreadsPorBloco,nThreadsPorBloco);
	dim3 numBlocosI(numBlocosLinha,numBlocosColuna);

	dev_Subtrai_Matrizes<<<numBlocosI,numThreadsI>>>(
				devHamiltoniano,
				dev_Matriz_Auxiliar_OxO,
				numLinhas,
				numColunas,
				dev_Matriz_Auxiliar_OxO);

	/*
	float *host_Matriz_Auxiliar_OxO;
	unsigned int qtdeElementos = numLinhas*numColunas;
	host_Matriz_Auxiliar_OxO = (float *)malloc(qtdeElementos*sizeof(float));
	cudaErrorCheck(cudaMemcpy(host_Matriz_Auxiliar_OxO, dev_Matriz_Auxiliar_OxO, qtdeElementos*sizeof(float), cudaMemcpyDeviceToHost));
	printf("\n"); printf("H - rho*I = [");
	printf("\n"); imprimeMatriz(host_Matriz_Auxiliar_OxO, numLinhas, numColunas);	
	printf("\n"); printf("]");
	host_Matriz_Auxiliar_OxO = NULL; free(host_Matriz_Auxiliar_OxO);
	*/

	// 2*[H - rho*I] ---------------------------------------------------------------------

	/*
	printf("\n"); printf("==================================================");
	printf("\n"); printf("Passo 7: 2*[H - rho*I]");
	printf("\n"); printf("==================================================");
	*/
	
	numLinhas = host_ParametrosGA->numGenes;
	numColunas = host_ParametrosGA->numGenes;
	
	float *ptrDummy;
	ptrDummy = (float *)malloc(sizeof(float));

	multiplica_matriz_por_escalar(
			dev_Matriz_Auxiliar_OxO,
			numLinhas,
			numColunas,
			0, // 0: o escalar está na CPU
			ptrDummy, // coloquei um valor arbitrário, já que utilizarei o 2 como escalar na cpu
			2, 
			dev_Matriz_Auxiliar_OxO);
	
	/*
	host_Matriz_Auxiliar_OxO = (float *)malloc(qtdeElementos*sizeof(float));
	cudaErrorCheck(cudaMemcpy(host_Matriz_Auxiliar_OxO, dev_Matriz_Auxiliar_OxO, qtdeElementos*sizeof(float), cudaMemcpyDeviceToHost));
	printf("\n"); printf("2*[H - rho*I] = [");
	printf("\n"); imprimeMatriz(host_Matriz_Auxiliar_OxO, numLinhas, numColunas);	
	printf("\n"); printf("]");
	host_Matriz_Auxiliar_OxO = NULL; free(host_Matriz_Auxiliar_OxO);
	ptrDummy = NULL; free(ptrDummy);

	*/

	// 2*[H - rho*I]*C

	/*
	printf("\n"); printf("==================================================");
	printf("\n"); printf("Passo 8: 2*[H - rho*I]*C");
	printf("\n"); printf("==================================================");
	*/

	MatMul(	dev_Matriz_Auxiliar_OxO, // mAuxOxO,
				host_ParametrosGA->numGenes, // ordem do hamiltoniano,
				host_ParametrosGA->numGenes,
				devGeracao->individuo[host_iIndividuo].gene, // Autovetor_C,
				host_ParametrosGA->numGenes, // Ordem_Hamiltoniano,
				1,
				dev_Matriz_Auxiliar_Ox1,
				host_ParametrosGA
	);

	/*
	float *host_Matriz_Auxiliar_Ox1;
	numLinhas = host_ParametrosGA->numGenes;
	numColunas = 1;
	qtdeElementos = numLinhas * numColunas;
	host_Matriz_Auxiliar_Ox1 = (float *)malloc(qtdeElementos*sizeof(float));
	cudaErrorCheck(cudaMemcpy(host_Matriz_Auxiliar_Ox1, dev_Matriz_Auxiliar_Ox1, qtdeElementos*sizeof(float), cudaMemcpyDeviceToHost));
	printf("\n"); printf("2*(H - rho*I)*C = [");
	printf("\n"); imprimeMatriz(host_Matriz_Auxiliar_Ox1, numLinhas, numColunas);	
	printf("\n"); printf("]");
	host_Matriz_Auxiliar_Ox1 = NULL; free(host_Matriz_Auxiliar_Ox1);
	*/

	// ==========================================================
	// PASSO 9 - Calculo do inverso de CtC
	// ==========================================================
	// OBS.: O teste do passo 9 foi feito imprimindo toda a geração antes de depois
	// do GradRho.
	
	dev_GradRho_CalculaInversoCtC<<<1,16>>>(devGeracao, host_iIndividuo);

	// passo 10

	/*
	printf("\n"); printf("============================================================");
	printf("\n"); printf("Passo 10: Gradiente de Rho (2*[H - rho*I]*C * (Ct*C)^-1) ");
	printf("\n"); printf("============================================================");
	*/

	numLinhas = host_ParametrosGA->numGenes;
	numColunas = 1;

	multiplica_matriz_por_escalar(
			dev_Matriz_Auxiliar_Ox1,
			numLinhas,
			numColunas,
			1, // 1: o escalar está na GPU
			&devGeracao->individuo[host_iIndividuo].inverso_de_CtC,
			-99, // arbitrário pois o valor escalar está na gpu. 
			dev_Matriz_Auxiliar_Ox1);
	
	// float *host_Matriz_Auxiliar_Ox1; - já foi declarada no escopo
	/*
	qtdeElementos = numLinhas * numColunas;
	host_Matriz_Auxiliar_Ox1 = (float *)malloc(qtdeElementos*sizeof(float));
	cudaErrorCheck(cudaMemcpy(host_Matriz_Auxiliar_Ox1, dev_Matriz_Auxiliar_Ox1, qtdeElementos*sizeof(float), cudaMemcpyDeviceToHost));
	printf("\n"); printf("Gradiente do individuo [%d]: ", host_iIndividuo);
	printf("\n"); imprimeMatriz(host_Matriz_Auxiliar_Ox1, numLinhas, numColunas);	
	printf("\n"); printf("]");
	host_Matriz_Auxiliar_Ox1 = NULL; free(host_Matriz_Auxiliar_Ox1);
	*/

	//printf("\n"); printf("============================================================");
	//printf("\n"); printf("Passo 11: modulo quadrado do gradiente de rho;
	//printf("\n"); printf("============================================================");

	MatMul(	dev_Matriz_Auxiliar_Ox1,
				1,
				host_ParametrosGA->numGenes,
				dev_Matriz_Auxiliar_Ox1,
				host_ParametrosGA->numGenes,
				1,
				&devGeracao->individuo[host_iIndividuo].grad_elevado_ao_quadrado,
				host_ParametrosGA
	);
}

//*******************************************************************************************************
__global__ void devFitness(float *devHamiltoniano,
									struct generation *devGeracao,
									struct parametros *devParametrosGA,
									struct parametros_Metodo *devParametrosMetodo,
									float *dev_Matriz_Auxiliar_OxO,
									float *dev_Matriz_Auxiliar_Ox1,
									float *dev_Matriz_Identidade) {
//*******************************************************************************************************

	unsigned int iIndividuo = threadIdx.x + blockIdx.x*blockDim.x;
	
	if (iIndividuo < devParametrosGA->numIndividuos) {

		float	grad_elevado_ao_quadrado,
				rhos_ao_quadrado,
				cociente_Rayleigh,
				fitness,
				Numerador,
				CtC,
				local_rho_minimo = devParametrosMetodo->rho_minimo,
				local_lambda = devParametrosMetodo->lambda;

		unsigned short int local_numGenes = devParametrosGA->numGenes;


		/*-----------------------------------------------
		** Cálculo do cociente de Rayleigh - Início
		**-----------------------------------------------
		*/

		//rho = [Ct * H * C  ]  / [ Ct * C ] 
		//
		//	onde
		//		C = vetor C (candidato à autovetor).
		//		Ct = transposto do vetor C.
		//		H = hamiltoniano.g

		Numerador	= 0;
		CtC			= 0;

	// Numerador [Ct * H * C]

		//Ct * H
		dev_Multiplica_Matrizes(
					devGeracao->individuo[iIndividuo].gene, // Ct
					1,  // número de linhas de Ct
					local_numGenes, // número de colunas de Ct (ordem do hamiltoniano)
					devHamiltoniano, // Hamiltoniano H
					local_numGenes, // número de linhas de H
					local_numGenes, // número de linhas de H
					(dev_Matriz_Auxiliar_Ox1 + iIndividuo * local_numGenes) // resultado da multiplicação
		);
		
		// Ct * H * C
		dev_Multiplica_Matrizes(
					(dev_Matriz_Auxiliar_Ox1 + iIndividuo*local_numGenes), // Ct * H
					1, // número de linhas de Ct * H
					local_numGenes, // número de colunas de Ct * H (ordem do hamiltoniano)
					devGeracao->individuo[iIndividuo].gene, // vetor C
					local_numGenes, // número de linhas de C
					1, // número de colunas de C
					&Numerador
		);

		// Denominador [ Ct * C ]
		dev_Multiplica_Matrizes(
					devGeracao->individuo[iIndividuo].gene,
					1, 
					local_numGenes,
					devGeracao->individuo[iIndividuo].gene,
					local_numGenes,
					1,
					&CtC
		);

		cociente_Rayleigh = (float)(Numerador/CtC);
		
		/*-----------------------------------------------
		** Cálculo do cociente de Rayleigh - Fim
		**-----------------------------------------------
		*/

		/*-----------------------------------------------
		** Cálculo do Gradiente de Rho - Início
		**-----------------------------------------------
		
		2*[H - rho*I]*C / Ct*C

		*/
		
		// rho*I
		dev_multiplica_matriz_por_escalar(
					dev_Matriz_Identidade,
					cociente_Rayleigh,
					local_numGenes*local_numGenes, // número de elementos 
					(dev_Matriz_Auxiliar_OxO + iIndividuo*local_numGenes*local_numGenes) // mAuxOxO);
		);

		// H - rho*I
		dev_Subtrai_Matrizes(
					devHamiltoniano,
					(dev_Matriz_Auxiliar_OxO + iIndividuo*local_numGenes*local_numGenes), // mAuxOxO,
					local_numGenes*local_numGenes, // número de elementos
					(dev_Matriz_Auxiliar_OxO + iIndividuo*local_numGenes*local_numGenes)
		);

		// 2*[H - rho*I]
		dev_multiplica_matriz_por_escalar(
					(dev_Matriz_Auxiliar_OxO + iIndividuo*local_numGenes*local_numGenes),
					2,
					local_numGenes*local_numGenes, // número de elementos
					(dev_Matriz_Auxiliar_OxO + iIndividuo*local_numGenes*local_numGenes)
		);

		// 2*[H - rho*I]*C
		dev_Multiplica_Matrizes(
					(dev_Matriz_Auxiliar_OxO + iIndividuo*local_numGenes*local_numGenes), // mAuxOxO,
					local_numGenes, // ordem do hamiltoniano,
					local_numGenes,
					devGeracao->individuo[iIndividuo].gene, // Autovetor_C,
					local_numGenes, // Ordem_Hamiltoniano,
					1,
					(dev_Matriz_Auxiliar_Ox1 + iIndividuo*local_numGenes)
		);

		// dividindo 2*[H - rho*I]*C por Ct*C
		dev_multiplica_matriz_por_escalar(
					(dev_Matriz_Auxiliar_Ox1 + iIndividuo*local_numGenes), // mAuxOx1, 
					(float)(1/CtC), // o CtC já calculei lá em cima...
					local_numGenes, // Ordem_Hamiltoniano, 
					&devGeracao->individuo[iIndividuo].gradRho[0]
		); //Gradiente_de_rho);


		/*-----------------------------------------------
		** Cálculo do Gradiente de Rho - Fim
		**-----------------------------------------------
		*/

		dev_Multiplica_Matrizes(&devGeracao->individuo[iIndividuo].gradRho[0],
										1,
										local_numGenes,
										&devGeracao->individuo[iIndividuo].gradRho[0],
										local_numGenes,
										1,
										&grad_elevado_ao_quadrado);

		rhos_ao_quadrado = pow( (cociente_Rayleigh - local_rho_minimo), 2);

		fitness = exp( (-1)*local_lambda*(rhos_ao_quadrado + grad_elevado_ao_quadrado));

		// Copia informações para a memória global
		devGeracao->individuo[iIndividuo].cociente_Rayleigh				= cociente_Rayleigh;
		devGeracao->individuo[iIndividuo].grad_elevado_ao_quadrado		= grad_elevado_ao_quadrado;
		devGeracao->individuo[iIndividuo].rho_menos_rho0_ao_quadrado	= rhos_ao_quadrado;
		devGeracao->individuo[iIndividuo].fitness								= fitness;
	}
}

//=======================================================================================
__global__ void dev_CalculaFitnessMedio_kernel(struct generation *dev_geracao,
															  struct parametros *dev_parametrosGA) {
//=======================================================================================

	unsigned short int iThread = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned short int iIndividuo, local_numIndividuos;

	local_numIndividuos = dev_parametrosGA->numIndividuos;

	if (iThread == 0) {

		dev_geracao->sumFitness = 0;
		dev_geracao->Maior_fitness = - 1.0F;
		dev_geracao->FitnessMedio = - 1.0F;
		
		for (iIndividuo = 0; iIndividuo < local_numIndividuos; iIndividuo++) {

			dev_geracao->sumFitness += dev_geracao->individuo[iIndividuo].fitness;

			if (dev_geracao->individuo[iIndividuo].fitness > dev_geracao->Maior_fitness) {

				dev_geracao->Maior_fitness = dev_geracao->individuo[iIndividuo].fitness;
				dev_geracao->Melhor_cociente_Rayleigh = dev_geracao->individuo[iIndividuo].cociente_Rayleigh;
				dev_geracao->idxMelhorIndividuo = iIndividuo;

			} //  do if (geracao->individuo

		} //  do for (iIndividuo = 0

	dev_geracao->FitnessMedio = dev_geracao->sumFitness / local_numIndividuos;

	} //  do if (iThread == 0)
}

//=======================================================================================
void CalculaFitnessMedio(	struct generation *dev_geracao,
									struct parametros *dev_parametrosGA) {
//=======================================================================================
	
	dev_CalculaFitnessMedio_kernel<<<1,32>>>(dev_geracao, dev_parametrosGA);

}

//=======================================================================================
void Fitness(	char tipoFitness,
					float *devHamiltoniano,
					struct generation *devGeracao0,
					struct parametros *hostParametrosGA,
					struct parametros *devParametrosGA,
					struct parametros_Metodo *devParametrosMetodo,
					float *devMatrizIdentidade) {
//=======================================================================================
		
	unsigned long int qtdeElementos;

	struct generation host_Geracao;

	float *devMatrizAuxiliarNxN, *devMatrizAuxiliarNx1;

	if (tipoFitness == '1') {

		//printf("\n"); printf("Fitness tipo 1");
		
		qtdeElementos = hostParametrosGA->numIndividuos * hostParametrosGA->numGenes * hostParametrosGA->numGenes;
		cudaErrorCheck(cudaMalloc((void **)&devMatrizAuxiliarNxN, qtdeElementos*sizeof(float)));

		qtdeElementos = hostParametrosGA->numIndividuos * hostParametrosGA->numGenes;
		cudaErrorCheck(cudaMalloc((void **)&devMatrizAuxiliarNx1, qtdeElementos*sizeof(float)));

		unsigned short int	num_Threads_por_Bloco = 16;
		unsigned short int	num_Blocos = (hostParametrosGA->numIndividuos + num_Threads_por_Bloco)/num_Threads_por_Bloco;

		dim3 BlocosF(num_Blocos);
		dim3 ThreadsF(num_Threads_por_Bloco);

		devFitness<<<BlocosF,ThreadsF>>>(	devHamiltoniano,
														devGeracao0,
														devParametrosGA,
														devParametrosMetodo,
														devMatrizAuxiliarNxN,
														devMatrizAuxiliarNx1,
														devMatrizIdentidade);
	} else if (tipoFitness == '2') {
		
		printf("\n"); printf("Fitness tipo 2");

		qtdeElementos = hostParametrosGA->numGenes * hostParametrosGA->numGenes;
		cudaErrorCheck(cudaMalloc((void **)&devMatrizAuxiliarNxN, qtdeElementos*sizeof(float)));
		
		qtdeElementos = hostParametrosGA->numGenes;
		cudaErrorCheck(cudaMalloc((void **)&devMatrizAuxiliarNx1, qtdeElementos*sizeof(float)));

		unsigned int iIndividuo;
		
		for (iIndividuo = 0; iIndividuo < hostParametrosGA->numIndividuos; iIndividuo++) {

			printf("\n"); printf("PASSO 01\tGeracao inicial antes da dev_CalculaRho: ");
			cudaErrorCheck( cudaMemcpy(&host_Geracao, devGeracao0, sizeof(struct generation), cudaMemcpyDeviceToHost) );
			imprimeGeracao(&host_Geracao, hostParametrosGA);

			dev_CalculaRho(	iIndividuo,
									devHamiltoniano,
									devGeracao0,
									hostParametrosGA,
									devParametrosGA,
									devParametrosMetodo,
									devMatrizAuxiliarNxN,
									devMatrizAuxiliarNx1,
									devMatrizIdentidade
			);

			
			printf("\n"); printf("PASSO N\tGeracao inicial depois da dev_CalculaRho: ");
			cudaErrorCheck( cudaMemcpy(&host_Geracao, devGeracao0, sizeof(struct generation), cudaMemcpyDeviceToHost) );
			imprimeGeracao(&host_Geracao, hostParametrosGA);
			
			GradRho(	iIndividuo,
						devHamiltoniano,
						devGeracao0,
						hostParametrosGA,
						devParametrosGA,
						devParametrosMetodo,
						devMatrizAuxiliarNxN,
						devMatrizAuxiliarNx1,
						devMatrizIdentidade
					);
		}

		/*
		printf("\n"); printf("==================================================");
		printf("\n"); printf("Passo 9: Calculo do inverso de CtC: geracao depois");
		printf("\n"); printf("==================================================");
		

		cudaErrorCheck( cudaMemcpy(&geracao[1], devGeracao0, sizeof(struct generation), cudaMemcpyDeviceToHost) );
		printf("\n"); imprimeGeracao(&geracao[1], &parametrosGA);

		printf("\n"); printf("============================================================");
		printf("\n"); printf("Passo 11: modulo quadrado do gradiente de rho");
		printf("\n"); printf("============================================================");
		
		cudaErrorCheck( cudaMemcpy(&geracao[1], devGeracao0, sizeof(struct generation), cudaMemcpyDeviceToHost) );
		printf("\n"); imprimeGeracao(&geracao[1], &parametrosGA);

		*/

		/*
		printf("\n"); printf("============================================================");
		printf("\n"); printf("Passo 12: calculo do fitness");
		printf("\n"); printf("============================================================");
		
		struct generation hostGeracaoAux;
		*/

		unsigned short int numBlocosRho = (hostParametrosGA->numIndividuos + nThreadsPorBloco)/nThreadsPorBloco;
		dim3 ThreadsRhoInd(nThreadsPorBloco);
		dim3 BlocosRhoInd(numBlocosRho);
		dev_CalculaFitness<<<BlocosRhoInd,ThreadsRhoInd>>>(devGeracao0, devParametrosGA, devParametrosMetodo);

		/*
		cudaErrorCheck( cudaMemcpy(&hostGeracaoAux, devGeracao0, sizeof(struct generation), cudaMemcpyDeviceToHost) );
		imprimeGeracao(&hostGeracaoAux, hostParametrosGA);
		*/
		
	} // final do if principal
	
	cudaErrorCheck(cudaFree(devMatrizAuxiliarNx1));
	cudaErrorCheck(cudaFree(devMatrizAuxiliarNxN));

	// Calcula o Fitness Médio e o Melhor Fitness
	CalculaFitnessMedio(devGeracao0, devParametrosGA);	
}


//=======================================================================================
__global__ void dev_InicializaAtributosGeracao(struct generation *dev_Geracao,
															  struct parametros *dev_ParametrosGA) {
//=======================================================================================
	unsigned int iThread = threadIdx.x + blockIdx.x*blockDim.x;

	if (iThread == 0) {
		dev_Geracao->FitnessMedio = -7.0F;
		dev_Geracao->sumFitness = -7.0F;;
		dev_Geracao->Melhor_cociente_Rayleigh = -7.0F;
		
		dev_Geracao->FitnessMedio = -7.0F;
		dev_Geracao->idxMelhorIndividuo = 0;
		dev_Geracao->Maior_fitness = -7.0F;
		dev_Geracao->Melhor_cociente_Rayleigh = -7.0F;
		dev_Geracao->sumFitness = -7.0F;
	}
}

//=======================================================================================
__global__ void dev_InicializaGenes(struct generation *dev_GeracaoInicial,
												struct parametros *dev_ParametrosGA,
												curandState *dev_estadosCURAND) {
//=======================================================================================

	unsigned int iIndividuo	= blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int iGene		= blockIdx.y * blockDim.y + threadIdx.y;

	unsigned int numIndividuos = dev_ParametrosGA->numIndividuos;
	unsigned int numGenes = dev_ParametrosGA->numGenes;

	if ((iIndividuo < numIndividuos) && (iGene < numGenes)) {

		unsigned int L;
		short int sinal;
		curandState iEstadoCURAND_noGene = *(dev_estadosCURAND + dev_Get_Idx_CURAND(iGene, numGenes, iIndividuo) );
		// ... aqui linearizo todos os genes de todos os indivíduos para
		// acessar globalmente o vetor de estados. Assim, cada gene
		// terá o seu próprio estado cuRAND, necessário para gerar
		// séries de números pesudo-aleatórios independentes.

		L = 2*curand_uniform(&iEstadoCURAND_noGene);
		
		sinal = (short int)pow(-1.0F, (float)L);
		
		dev_GeracaoInicial->individuo[iIndividuo].gene[iGene] = (float)sinal * curand_uniform(&iEstadoCURAND_noGene);		

		__syncthreads();

		if (iGene == 0) {
			float vlrInicializacao = -3.0F;
			dev_GeracaoInicial->individuo[iIndividuo].cociente_Rayleigh = vlrInicializacao;
			//dev_GeracaoInicial->individuo[iIndividuo].gradRho = vlrInicializacao;
			dev_GeracaoInicial->individuo[iIndividuo].grad_elevado_ao_quadrado = vlrInicializacao;
			dev_GeracaoInicial->individuo[iIndividuo].rho_menos_rho0_ao_quadrado = vlrInicializacao;
			dev_GeracaoInicial->individuo[iIndividuo].Numerador = vlrInicializacao;
			dev_GeracaoInicial->individuo[iIndividuo].CtC = vlrInicializacao;
			dev_GeracaoInicial->individuo[iIndividuo].inverso_de_CtC = vlrInicializacao;
			dev_GeracaoInicial->individuo[iIndividuo].fitness = vlrInicializacao;
		}

		*(dev_estadosCURAND + dev_Get_Idx_CURAND(iGene, numGenes, iIndividuo)) = iEstadoCURAND_noGene;
	}
}

//=======================================================================================
void GeraPopulacaoInicial_paralela(	struct generation *dev_GeracaoInicial,
												struct parametros *host_parametrosGA,
												struct parametros *dev_parametrosGA,
												curandState *dev_estadosCURAND) {
//=======================================================================================
	// Como atuarei em apenas uma geração, os atributos específicos da geração devem ser
	// tratados de maneira serial, mas, diretamente na GPU
	dev_InicializaAtributosGeracao<<<1,16>>>(dev_GeracaoInicial, dev_parametrosGA);

	unsigned int numThreadsPorBloco = 16;
	unsigned int numIndividuos = host_parametrosGA->numIndividuos;
	unsigned int numGenes = host_parametrosGA->numGenes;
	unsigned int numBlocosIndividuos = (numIndividuos + numThreadsPorBloco)/numThreadsPorBloco;
	unsigned int numBlocosGenes = (numGenes + numThreadsPorBloco)/numThreadsPorBloco;

	dim3 ThreadsPorBloco(numThreadsPorBloco,numThreadsPorBloco);
	dim3 numBlocos(numBlocosIndividuos, numBlocosGenes);

	dev_InicializaGenes<<<numBlocos,ThreadsPorBloco>>>(
						dev_GeracaoInicial,
						dev_parametrosGA,
						dev_estadosCURAND);
}
//===========================================================================================
__global__ void dev_Gera_Matriz_de_Coope_naGPU(
							float* dev_H,
							struct parametros *dev_ParametrosGA) {
//===========================================================================================
	
	unsigned int iLinha = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int iColuna = threadIdx.y  + blockIdx.y*blockDim.y;
	
	unsigned int numLinhas = dev_ParametrosGA->numGenes;
	unsigned int numColunas = dev_ParametrosGA->numGenes;

	if ( (iLinha < numLinhas) && (iColuna < numColunas)) {	
		unsigned long int Ordem = dev_ParametrosGA->numGenes;
		float valor;
		if (iLinha == iColuna) {
			valor = (float)(2*(iLinha+1) - 1);
			//*(Matriz + Ordem*iLinha + iColuna) = (float) -1 * (20 - iLinha);
				// a definição acima gera uma matriz com todos os autovalores negativos.
		}
		else {
			valor = 1.0F;
		}
		*(dev_H + Ordem*iLinha + iColuna) = valor;
	}
}




//========================================================================================
__global__ void devSelecaoViaTorneio2(
								struct generation *PopulacaoAntes,
								struct generation *PopulacaoDepois,
								struct parametros *parametrosGA,
								curandState *estadoCURAND,
								struct testeGeracao_s *dev_testeSelecao) {
//========================================================================================

	unsigned short int iIndividuo = threadIdx.x + blockIdx.x*blockDim.x;

	if (iIndividuo < parametrosGA->numIndividuos) {

		//unsigned int iIndividuo;
		char iTamanhoDoTorneio;
		unsigned short int iIndividuo_para_Torneio;
		struct individual *Individuo;
		float melhorFitness;
		
		// Faço uma cópia local do stado do CURAND para fins de eficiência
		// Conselho encontrado na documentação do CURAND.
		curandState EstadoCURANDlocal = *(estadoCURAND + dev_Get_Idx_CURAND(0, parametrosGA->numGenes, iIndividuo));

		melhorFitness = -1.0F; // pra sempre garantir que o primeiro indivíduo será o melhor, garantindo
										// que, pelo menos, o primeiro indivíduos será levado pra frente.

		for (iTamanhoDoTorneio = 0; iTamanhoDoTorneio < parametrosGA->tamanho_torneio; iTamanhoDoTorneio++) {
			
			iIndividuo_para_Torneio = devRandomicoUI(
													parametrosGA->numIndividuos,
													&EstadoCURANDlocal);

			// teste - pode ser removido depois
			dev_testeSelecao->individuo[iIndividuo].teste_iIndividuo_para_Torneio[iTamanhoDoTorneio] = iIndividuo_para_Torneio;
			
			Individuo = &PopulacaoAntes->individuo[iIndividuo_para_Torneio];

			if (Individuo->fitness > melhorFitness) {
				melhorFitness = Individuo->fitness;
				PopulacaoDepois->individuo[iIndividuo] = *Individuo;
			};
		};

		//PopulacaoDepois->individuo[iIndividuo] = Melhor;

		// Retorno o estado local para a memória global.
		*(estadoCURAND + dev_Get_Idx_CURAND(0, parametrosGA->numGenes, iIndividuo)) = EstadoCURANDlocal;
	}
}

//========================================================================================
void Selecao(	struct generation *dev_PopulacaoAntes,
					struct generation *dev_PopulacaoDepois,
					struct parametros *host_parametrosGA,
					struct parametros *dev_parametrosGA,
					curandState *dev_estadoCURAND,
					struct testeGeracao_s *dev_testeSelecao) {
//========================================================================================

	unsigned int numBlocosX = (host_parametrosGA->numIndividuos + nThreadsPorBloco)/nThreadsPorBloco;

	dim3 ThreadsI(nThreadsPorBloco);
	dim3 BlocosI(numBlocosX);

	devSelecaoViaTorneio2<<<BlocosI,ThreadsI>>>(
						dev_PopulacaoAntes,
						dev_PopulacaoDepois,
						dev_parametrosGA,
						dev_estadoCURAND,
						dev_testeSelecao);
}

//===========================================================================================
__device__ unsigned int dev_get_iGene_Global(	unsigned int iGene,
																unsigned int iIndividuo,
																unsigned int numGenes) {
//===========================================================================================
	return iGene + numGenes*iIndividuo;	
}


//===========================================================================================
__global__ void devCrossOver1Ponto_paralelo_naGPU_1d(
							struct generation *g0,
							struct generation *g1,
							struct parametros *parametrosGA,
							curandState *dev_estadosCURAND,
							struct testeGeracao_s *dev_TesteGeracao_s) {
//===========================================================================================
	
	unsigned short int iIndividuo = threadIdx.x + blockIdx.x*blockDim.x;

	unsigned short int local_numGenes = parametrosGA->numGenes;
	unsigned short int local_numIndividuos = parametrosGA->numIndividuos;
	float local_probCrossOver = parametrosGA->probCrossOver;

	unsigned short int iGene, ponto_de_corte;
	
	unsigned short int idx_Mae, idx_Pai;
	struct individual *Pai, *Mae;

	float f, p_aux;	// p = probabilidade auxiliar para comparacao com a probabilidade do crossover

	curandState i_EstadoCURAND_noIndividuo = *(dev_estadosCURAND + dev_Get_Idx_CURAND(0, local_numGenes, iIndividuo));

	//for (iIndividuo = 0; iIndividuo < parametrosGA->numIndividuos; iIndividuo = iIndividuo + 1) {	
	if (iIndividuo < parametrosGA->numIndividuos) {
		idx_Pai = iIndividuo;		
		idx_Mae = devRandomicoUI(local_numIndividuos, &i_EstadoCURAND_noIndividuo);
		dev_TesteGeracao_s->individuo[iIndividuo].idx_Mae = idx_Mae;

		Pai = &g0->individuo[idx_Pai];
		Mae = &g0->individuo[idx_Mae];

		p_aux = curand_uniform(&i_EstadoCURAND_noIndividuo);
		dev_TesteGeracao_s->individuo[iIndividuo].p = p_aux;

		f = -1.0F;
		dev_TesteGeracao_s->individuo[iIndividuo].f = f;

		ponto_de_corte = devRandomicoUI(local_numGenes - 1, &i_EstadoCURAND_noIndividuo);
		dev_TesteGeracao_s->individuo[iIndividuo].pontos_de_corte[0] = ponto_de_corte;
		dev_TesteGeracao_s->individuo[iIndividuo].pontos_de_corte[1] = 0;

		f = curand_uniform(&i_EstadoCURAND_noIndividuo);
		dev_TesteGeracao_s->individuo[iIndividuo].f = f;

		if (p_aux <= local_probCrossOver) {
			
			for (iGene = 0; iGene < ponto_de_corte; iGene++) {
				g1->individuo[iIndividuo].gene[iGene] = Pai->gene[iGene];
			}
			
			for (iGene = ponto_de_corte; iGene < parametrosGA->numGenes; iGene++) {				
				g1->individuo[iIndividuo].gene[iGene] = f*Pai->gene[iGene] + (1 - f)*Mae->gene[iGene];
			}
		}
		else {
			// Não houve crossover.
			// Todos os genes dos indivíduos originais são copiados para os indivíduos da próxima geração.
			for (iGene = 0; iGene < parametrosGA->numGenes; iGene++) {
				g1->individuo[iIndividuo].gene[iGene] =	Pai->gene[iGene];
			}
		}
	} // do for no serial e do if no paralelo
	*(dev_estadosCURAND + dev_Get_Idx_CURAND(0, parametrosGA->numGenes, iIndividuo)) = i_EstadoCURAND_noIndividuo;
}

//===========================================================================================
__global__ void devCrossOver1Ponto_paralelo_naGPU_2d(
							struct generation *g0,
							struct generation *g1,
							struct parametros *parametrosGA,
							curandState *dev_estadosCURAND,
							struct testeGeracao_s *dev_TesteGeracao_s) {
//===========================================================================================

	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	// Essa versão 2d não está funcionando e não deve ser utilizada.
	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	unsigned short int iIndividuo = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned short int iGene = threadIdx.y + blockIdx.y*blockDim.y;
	
	unsigned short int local_numGenes = parametrosGA->numGenes;
	unsigned short int local_numIndividuos = parametrosGA->numIndividuos;
	float local_probCrossOver = parametrosGA->probCrossOver;

	unsigned short int idx_Mae, idx_Pai, ponto_de_corte;

	struct individual *Pai, *Mae;

	float f, p_aux;	// p_aux = probabilidade auxiliar para comparacao com a probabilidade do crossover
	
	curandState i_EstadoCURAND_noIndividuo = *(dev_estadosCURAND + dev_Get_Idx_CURAND(0, local_numGenes, iIndividuo));

	if ( (iIndividuo < local_numIndividuos) && (iGene < local_numGenes) ) {
		
		idx_Pai = iIndividuo;

		idx_Mae = devRandomicoUI(local_numIndividuos, &i_EstadoCURAND_noIndividuo);
		dev_TesteGeracao_s->individuo[iIndividuo].idx_Mae = idx_Mae;

		Pai = &g0->individuo[idx_Pai];
		Mae = &g0->individuo[idx_Mae];

		p_aux = curand_uniform(&i_EstadoCURAND_noIndividuo);
		dev_TesteGeracao_s->individuo[iIndividuo].p = p_aux;

		f = -1.0F;
		dev_TesteGeracao_s->individuo[iIndividuo].f = f;

		// GERA PONTO DE CORTE
		ponto_de_corte = devRandomicoUI(local_numGenes - 1, &i_EstadoCURAND_noIndividuo);
		
		//g0->individuo[iIndividuo].pontos_de_corte[0] = ponto_de_corte;
			// Acho que a linha acima não é necessária... A variável ponto_de_corte é definida localmente.
			// Importante: ela só não é necessária aqui pois o crossover é de um ponto!!

		dev_TesteGeracao_s->individuo[iIndividuo].pontos_de_corte[0] = ponto_de_corte;
		dev_TesteGeracao_s->individuo[iIndividuo].pontos_de_corte[1] = 0;

		if (p_aux <= local_probCrossOver) {
			
			f = curand_uniform(&i_EstadoCURAND_noIndividuo);

			if (iGene < ponto_de_corte) {
				g1->individuo[iIndividuo].gene[iGene] = g0->individuo[iIndividuo].gene[iGene];
			}
			else {
				g1->individuo[iIndividuo].gene[iGene] = (float) f*(g0->individuo[iIndividuo].gene[iGene]) + (1 - f)*(g0->individuo[idx_Mae].gene[iGene]);
			}

		}
		else {
				g1->individuo[iIndividuo].gene[iGene] = g0->individuo[iIndividuo].gene[iGene];
		}

	*(dev_estadosCURAND + dev_Get_Idx_CURAND(0, parametrosGA->numGenes, iIndividuo)) = i_EstadoCURAND_noIndividuo;

	}
}

//===========================================================================================
void CrossOver1Ponto(
							struct generation *dev_ger_Antes,
							struct generation *dev_ger_Depois,
							struct parametros *host_parametrosGA,
							struct parametros *dev_parametrosGA,
							curandState *dev_estadosCURAND,
							struct testeGeracao_s *dev_TesteGeracao_s) {
//===========================================================================================

	unsigned int num_Threads_por_Bloco_Cross = 16;
	unsigned int numBlocosX = (host_parametrosGA->numIndividuos + num_Threads_por_Bloco_Cross)/num_Threads_por_Bloco_Cross;
	
	dim3 ThreadsC(num_Threads_por_Bloco_Cross);
	dim3 BlocosC(numBlocosX);

	devCrossOver1Ponto_paralelo_naGPU_1d<<<BlocosC,ThreadsC>>>(
			dev_ger_Antes,
			dev_ger_Depois,
			dev_parametrosGA,
			dev_estadosCURAND,
			dev_TesteGeracao_s);
}

//===========================================================================================
__global__ void devMutacao(
						struct generation *geracao,
						struct parametros *devParametrosGA,
						curandState *estadosCURAND) {
//===========================================================================================

	unsigned int iIndividuo = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int iGene = threadIdx.y + blockIdx.y*blockDim.y;
	
	unsigned int numIndividuos = devParametrosGA->numIndividuos;
	unsigned int numGenes = devParametrosGA->numGenes;
	
	if ((iGene < numGenes) && (iIndividuo < numIndividuos)) {
		/* 
		* A mutação seguiu a equação 10 do artigo "Diagonalization of a real-symmetric
		* Hamiltonian by genetic algorithm".
		*/
		//unsigned int iGene;
		unsigned int L;
		int termo_do_L;
		float r, pAux;
		float probMutacao = devParametrosGA->probMutacao;
		float intensidadeMutacao = devParametrosGA->intensidadeMutacao;
		
		curandState i_estadoCURAND_noGene = *(estadosCURAND + dev_Get_Idx_CURAND(iGene, numGenes, iIndividuo));

		pAux = curand_uniform(&i_estadoCURAND_noGene);
		//geracao->individuo[iIndividuo].gene_pAux[iGene] = pAux; // Teste. Pode ser removido.

		if (pAux <= probMutacao) {
			L = devRandomicoUI(2, &i_estadoCURAND_noGene);
			//geracao->individuo[iIndividuo].gene_L[iGene] = L; // Teste. Pode ser remo = vido.
				// No artigo esse L é definido como um número
				// randômico inteiro. Nada é dito sobre os limites
				// inferior ou superior de L. Decidi arbitrariamente
				// que L = (0, 2].
			r = curand_uniform(&i_estadoCURAND_noGene);
			//geracao->individuo[iIndividuo].gene_r[iGene] = r; // Teste. Pode ser removido.
				// r é um número randômico entre zero e 1. 
			termo_do_L = (int)pow((float)(-1),(int)L);
			//geracao->individuo[iIndividuo].gene_termo_do_L[iGene] = termo_do_L; // Teste. Pode ser removido.
				// é o (-1) elevado ao L. Só pode ser + ou - 1,
				// portanto, um inteiros.
			geracao->individuo[iIndividuo].gene[iGene] =
					geracao->individuo[iIndividuo].gene[iGene] +
					(float)(termo_do_L * r * intensidadeMutacao);
							// parametrosGA->intensidadeMutacao é o 
							// DELTA da equação 10 do artigo.
		}
		*(estadosCURAND + dev_Get_Idx_CURAND(iGene, numGenes, iIndividuo)) = i_estadoCURAND_noGene;
	}	
}

//===========================================================================================
void Mutacao(	struct generation *dev_geracao,
					struct parametros *host_parametrosGA,
					struct parametros *dev_parametrosGA,
					curandState *dev_estadosCURAND) {
//===========================================================================================
	unsigned int numBlocosX = (host_parametrosGA->numIndividuos + nThreadsPorBloco)/nThreadsPorBloco;
	unsigned int numBlocosY = (host_parametrosGA->numGenes + nThreadsPorBloco)/nThreadsPorBloco;

	dim3 ThreadsM(nThreadsPorBloco, nThreadsPorBloco);
	dim3 BlocosM(numBlocosX, numBlocosY);

	devMutacao<<<BlocosM,ThreadsM>>>(
		dev_geracao,
		dev_parametrosGA,
		dev_estadosCURAND);
}

//===========================================================================================
void testePopulacaoInicial (	struct generation *host_Geracao,
										struct generation *dev_Geracao,
										struct parametros *host_ParametrosGA) {
//===========================================================================================
	
	
	printf("\n\n+===================================================+");
	printf("\n+ TESTE POPULACAO INICIAL - INICIO                  +");
	printf("\n+===================================================+");
		
	//printf("\n\nPopulacao inicial na CPU: \n\n");
	//imprimeGeracao(host_Geracao, host_ParametrosGA);

	//printf("\n\nLimpando a da CPU... ");
	
	unsigned int iIndividuo, iGene;

	// Atributos específicos da geração
	host_Geracao->FitnessMedio = -1.0F;
	host_Geracao->idxMelhorIndividuo = 0;
	host_Geracao->Maior_fitness = -1.0F;
	host_Geracao->Melhor_cociente_Rayleigh = -1.0F;
	host_Geracao->sumFitness = -1.0F;
	
	// Indivíduos da geração
	for (iIndividuo = 0; iIndividuo < host_ParametrosGA->numIndividuos; iIndividuo++) {
		
		host_Geracao->individuo[iIndividuo].cociente_Rayleigh = -1.0F;
		host_Geracao->individuo[iIndividuo].grad_elevado_ao_quadrado = -1.0F;
		host_Geracao->individuo[iIndividuo].rho_menos_rho0_ao_quadrado = -1.0F;
		host_Geracao->individuo[iIndividuo].Numerador = -1.0F;
		host_Geracao->individuo[iIndividuo].CtC = -1.0F;
		host_Geracao->individuo[iIndividuo].inverso_de_CtC = -1.0F;
		host_Geracao->individuo[iIndividuo].fitness = -1.0F;
		
		for (iGene = 0; iGene < host_ParametrosGA->numGenes; iGene++) {
			host_Geracao->individuo[iIndividuo].gene[iGene] = - 1.0F;
			host_Geracao->individuo[iIndividuo].gradRho[iGene] = -1.0F;
		}
	}

	printf("Ok. ");

	//printf("\n"); printf("Populacao inicial depois da limpeza: ");
	//imprimeGeracao(host_Geracao, host_ParametrosGA);

	printf("\nCopiando a populacao da GPU pra CPU:");
	cudaMemcpy(host_Geracao, dev_Geracao, sizeof(struct generation), cudaMemcpyDeviceToHost);
	printf("Ok. ");

	printf("\nPopulacao depois da copia:");
	imprimeGeracao(host_Geracao, host_ParametrosGA);
	/*
	printf("\n\n+===================================================+");
	printf("\n+ TESTE POPULACAO INICIAL - FIM                     +");
	printf("\n+===================================================+");
	*/
}

//===========================================================================================
void testeParametrosGA(	struct parametros *host_ParametrosGA,
								struct parametros *dev_ParametrosGA) {
//===========================================================================================

	printf("\n\n+===================================================+");
	printf("\n+ TESTE PARAMETROS GA - INICIO                      +");
	printf("\n+===================================================+");

	printf("\n\nImprimindo os parâmetros do GA... ");
	printf("\nnumGenes = %d", host_ParametrosGA->numGenes);
	printf("\nnumIndividuos = %d", host_ParametrosGA->numIndividuos);
	printf("\nnumGeracoes = %d", host_ParametrosGA->numGeracoes);
	printf("\nprobCrossOver = %f", host_ParametrosGA->probCrossOver);
	printf("\nprobMutacao = %f", host_ParametrosGA->probMutacao);
	printf("\nIntensidade da Mutacao = %d", host_ParametrosGA->intensidadeMutacao);
	printf("\nTamanho do torneio = %d", host_ParametrosGA->tamanho_torneio);
	printf("\nQtde Pontos de corte = %d", host_ParametrosGA->qtde_Pontos_de_Corte);
	printf("\nValor inicial maximo do gene = %d", host_ParametrosGA->vlrMaximoGene);
	printf("\n\nImprimindo os parâmetros do GA... Ok.");

	printf("\n\nLimpando os parametros do GA... ");
	
	host_ParametrosGA->numGenes = 789;
	host_ParametrosGA->numIndividuos = 789;
	host_ParametrosGA->numGeracoes = 789;
	host_ParametrosGA->probCrossOver = 789;
	host_ParametrosGA->probMutacao = 789;
	host_ParametrosGA->intensidadeMutacao = 789;
	host_ParametrosGA->tamanho_torneio = 7;
	host_ParametrosGA->qtde_Pontos_de_Corte = 789;
	host_ParametrosGA->vlrMaximoGene = 789;

	printf("Ok.");

	printf("\n\nImprimindo os parâmetros do GA 'limpos'... ");
	printf("\nnumGenes = %d", host_ParametrosGA->numGenes);
	printf("\nnumIndividuos = %d", host_ParametrosGA->numIndividuos);
	printf("\nnumGeracoes = %d", host_ParametrosGA->numGeracoes);
	printf("\nprobCrossOver = %f", host_ParametrosGA->probCrossOver);
	printf("\nprobMutacao = %f", host_ParametrosGA->probMutacao);
	printf("\nIntensidade da Mutacao = %d", host_ParametrosGA->intensidadeMutacao);
	printf("\nTamanho do torneio = %d", host_ParametrosGA->tamanho_torneio);
	printf("\nQtde Pontos de corte = %d", host_ParametrosGA->qtde_Pontos_de_Corte);
	printf("\nValor inicial maximo do gene = %d", host_ParametrosGA->vlrMaximoGene);
	printf("\n\nImprimindo os parâmetros do GA 'limpos'... Ok.");
	
	printf("\n\nCopiando os parâmetros do GA da GPU pra CPU... ");
	cudaMemcpy(host_ParametrosGA, dev_ParametrosGA, sizeof(struct parametros), cudaMemcpyDeviceToHost);
	printf("Ok.");

	printf("\n\nImprimindo os parâmetros do GA vindos da GPU... ");
	printf("\nnumGenes = %d", host_ParametrosGA->numGenes);
	printf("\nnumIndividuos = %d", host_ParametrosGA->numIndividuos);
	printf("\nnumGeracoes = %d", host_ParametrosGA->numGeracoes);
	printf("\nprobCrossOver = %f", host_ParametrosGA->probCrossOver);
	printf("\nprobMutacao = %f", host_ParametrosGA->probMutacao);
	printf("\nIntensidade da Mutacao = %d", host_ParametrosGA->intensidadeMutacao);
	printf("\nTamanho do torneio = %d", host_ParametrosGA->tamanho_torneio);
	printf("\nQtde Pontos de corte = %d", host_ParametrosGA->qtde_Pontos_de_Corte);
	printf("\nValor inicial maximo do gene = %d", host_ParametrosGA->vlrMaximoGene);
	printf("\n\nImprimindo os parâmetros do GA vindos da GPU... Ok.");
	
	printf("\n\n+===================================================+");
	printf("\n+ TESTE PARAMETROS GA - FIM                         +");
	printf("\n+===================================================+");
}

//===========================================================================================
void testeParametrosMetodo(	struct parametros_Metodo *host_Parametros_Metodo,
										struct parametros_Metodo *dev_Parametros_Metodo) {
//===========================================================================================
	printf("\n\n+===================================================+");
	printf("\n+ TESTE PARAMETROS METODO - INICIO                  +");
	printf("\n+===================================================+");

	printf("\n\nImprimindo os atuais parametros do Metodo...\n");
	printf("\n...lambda = %f", host_Parametros_Metodo->lambda);
	printf("\n...rho_minimo = %f", host_Parametros_Metodo->rho_minimo);
	printf("\nOk.");

	printf("\n\nLimpando os parametros do Metodo na CPU... ");
	host_Parametros_Metodo->lambda = -2;
	host_Parametros_Metodo->rho_minimo = -2;
	printf("Ok.");

	printf("\n\nImprimindo os parametros do Metodo limpos...\n");
	printf("\n...lambda = %f", host_Parametros_Metodo->lambda);
	printf("\n...rho_minimo = %f", host_Parametros_Metodo->rho_minimo);
	printf("\nOk.");

	printf("\n\nCopiando os parametros do metodo da GPU para a CPU... ");
	cudaMemcpy(host_Parametros_Metodo, dev_Parametros_Metodo, sizeof(struct parametros_Metodo), cudaMemcpyDeviceToHost);
	printf("Ok.");

	printf("\n\nImprimindo os novos parametros do Metodo...\n");
	printf("\n...lambda = %f", host_Parametros_Metodo->lambda);
	printf("\n...rho_minimo = %f", host_Parametros_Metodo->rho_minimo);
	printf("\nOk.");

	printf("\n\n+===================================================+");
	printf("\n+ TESTE PARAMETROS METODO - FIM                     +");
	printf("\n+===================================================+");
}

//=======================================================================================
void testeMatrizIdentidade_GPU(	float *dev_MatrizIdentidade,
											struct parametros *host_ParametrosGA) {
//=======================================================================================

	printf("\n\n+===================================================+");
	printf("\n+ TESTE MATRIZ IDENTIDADE - INICIO                  +");
	printf("\n+===================================================+");

	//printf("\n\nMatriz Identidade CPU = \n");

	unsigned int iElemento, qtdeElementos, Ordem;
	float *host_MatrizIdentidade;
	
	Ordem = host_ParametrosGA->numGenes;
	qtdeElementos = Ordem*Ordem;
	host_MatrizIdentidade = (float *)malloc(qtdeElementos*sizeof(float));

	//for (iElemento = 0; iElemento < qtdeElementos; iElemento++) {
	//	if (iElemento % Ordem == 0) printf("\n"); 
	//	printf("%f\t", *(host_MatrizIdentidade + iElemento));
	//}
	
	printf("\n\nLimpando Matriz Identidade CPU... ");
	for (iElemento = 0; iElemento < qtdeElementos; iElemento++) {
		*(host_MatrizIdentidade + iElemento) = -1;
	}
	printf("Ok.");

	//printf("\n\nMatriz Identidade CPU limpa = \n");
	//for (iElemento = 0; iElemento < qtdeElementos; iElemento++) {
	//	if (iElemento % Ordem == 0) printf("\n"); 
	//	printf("%f\t", *(host_MatrizIdentidade + iElemento));
	//}
	
	printf("\n\nCopiando a matriz da GPU para a CPU... ");
	cudaErrorCheck(cudaMemcpy(host_MatrizIdentidade, dev_MatrizIdentidade, qtdeElementos*sizeof(float), cudaMemcpyDeviceToHost) );
	printf("Ok.");

	printf("\n\nMatriz Identidade CPU trazida da GPU = \n");
	for (iElemento = 0; iElemento < qtdeElementos; iElemento++) {
		if (iElemento % Ordem == 0) printf("\n"); 
		printf("%f\t", *(host_MatrizIdentidade + iElemento));
	}

	free(host_MatrizIdentidade);

	printf("\n\n+===================================================+");
	printf("\n+ TESTE MATRIZ IDENTIDADE - FIM                     +");
	printf("\n+===================================================+");
}


//=======================================================================================
__global__ void devPopulaPopGPU(	struct generation *dev_Geracao,
											struct parametros *dev_ParametrosGA) {
//=======================================================================================
	
	unsigned int iThread = threadIdx.x + blockIdx.x*blockDim.x;
	unsigned int t_iIndividuo, t_iGene;

	//unsigned int numGenes = dev_ParametrosGA->numGenes;

	if (iThread == 0) {
		
		dev_Geracao->FitnessMedio = -7.0F;
		dev_Geracao->idxMelhorIndividuo = 0;
		dev_Geracao->Maior_fitness = -7.0F;
		dev_Geracao->Melhor_cociente_Rayleigh = -7.0F;
		dev_Geracao->sumFitness = -7.0F;
		
		for (t_iIndividuo = 0; t_iIndividuo < dev_ParametrosGA->numIndividuos; t_iIndividuo++) {
			dev_Geracao->individuo[t_iIndividuo].cociente_Rayleigh = -7.0F;
			dev_Geracao->individuo[t_iIndividuo].grad_elevado_ao_quadrado = -7.0F;
			dev_Geracao->individuo[t_iIndividuo].rho_menos_rho0_ao_quadrado = -7.0F;
			dev_Geracao->individuo[t_iIndividuo].Numerador = -7.0F;
			dev_Geracao->individuo[t_iIndividuo].CtC = -7.0F;
			dev_Geracao->individuo[t_iIndividuo].inverso_de_CtC = -7.0F;
			dev_Geracao->individuo[t_iIndividuo].fitness = -7.0F;

			for (t_iGene = 0; t_iGene < dev_ParametrosGA->numGenes; t_iGene++) {
				dev_Geracao->individuo[t_iIndividuo].gene[t_iGene] = -7.0F;
				dev_Geracao->individuo[t_iIndividuo].gradRho[t_iGene] = -7.0F;
			}
		}
	}
	//__syncthreads();
		/*
			Imaginei que a sincronização de threads corrigisse um erro
			relacionado à memória. Mas não deu certo. Me parece que,
			nessa função, apesar do if acima, não é necessária a sin-
			cronização.

			Claro que não! como apenas uma thread executa a gravação
			do conteúdo de todos os genes de todos os indivíduos, não
			há risco dos dados serem utilizados antes do término da
			thread zero. Caso isso acontecesse, significaria que a fun-
			ção dev_(...) "terminou" deixando uma thread em execução.
		*/
}

//===========================================================================================
void testeMutacao(int flag_Momento,
						struct generation *dev_Geracao,
						struct parametros *host_parametros) {
//===========================================================================================
	
	struct generation host_Geracao;
	unsigned int iIndividuo, iGene;

	cudaErrorCheck(cudaMemcpy(&host_Geracao, dev_Geracao, sizeof(struct generation), cudaMemcpyDeviceToHost));

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
				printf("\t");	printf("%f", host_Geracao.individuo[iIndividuo].gene[iGene]);
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
				printf("\t");	printf("%f", host_Geracao.individuo[iIndividuo].gene[iGene]);
				//printf("\t");	printf("%f", host_Geracao.individuo[iIndividuo].gene_pAux[iGene]);
				//printf("\t");	printf("%d", host_Geracao.individuo[iIndividuo].gene_L[iGene]);
				//printf("\t");	printf("%f", host_Geracao.individuo[iIndividuo].gene_r[iGene]);
				//printf("\t");	printf("%d", host_Geracao.individuo[iIndividuo].gene_termo_do_L[iGene]);
			}
		}
	}
}


//===========================================================================================
void testeSelecao(unsigned int flag_Momento,
						struct generation *dev_Geracao,
						struct parametros *host_ParametrosGA,
						struct testeGeracao_s *dev_testeSelecao) {
//===========================================================================================

	struct generation host_Geracao;
	struct testeGeracao_s host_testeSelecao;

	cudaErrorCheck(cudaMemcpy(&host_Geracao, dev_Geracao, sizeof(struct generation), cudaMemcpyDeviceToHost));
	cudaErrorCheck(cudaMemcpy(&host_testeSelecao, dev_testeSelecao, sizeof(struct testeGeracao_s), cudaMemcpyDeviceToHost));

	if (flag_Momento == 0) {
		// antes da seleção atuar
		// apenas imprimo a geracao.
		printf("\n");	printf("Teste Selecao - Populacao Antes");
		printf("\n");	printf("-----------------------------------------------------------");
		imprimeGeracao(&host_Geracao, host_ParametrosGA);
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
				printf("\t");	printf("%d", host_testeSelecao.individuo[iIndividuo].teste_iIndividuo_para_Torneio[iTorneio]);
			}
			printf("\t");	printf("%d", host_testeSelecao.individuo[iIndividuo].iIndividuo_Vencedor);
		}
	}
}

//===========================================================================================
void testeCrossOver(	unsigned int flag_Momento,
							struct generation *dev_Geracao,
							struct parametros *host_ParametrosGA,
							struct testeGeracao_s *dev_testeGeracao) {
//===========================================================================================

	struct generation host_Geracao;
	cudaErrorCheck(cudaMemcpy(&host_Geracao, dev_Geracao, sizeof(struct generation), cudaMemcpyDeviceToHost));
	
	if (flag_Momento == 2) { // antes da aplicação da CrossOver
		printf("\n");	printf("Teste Crossover - Proxima Populacao Antes de receber os novos genes");
		printf("\n");	printf("-----------------------------------------------------------");

		imprimeGeracao(&host_Geracao, host_ParametrosGA);
	}

	if (flag_Momento == 0) { // antes da aplicação da CrossOver
		printf("\n");	printf("Teste Crossover - Populacao Antes");
		printf("\n");	printf("-----------------------------------------------------------");

		imprimeGeracao(&host_Geracao, host_ParametrosGA);
	}
	else if (flag_Momento == 1) { // depois da aplicação da CrossOver
		
		printf("\n");	printf("Teste Crossover - Termos auxiliares");
		printf("\n");	printf("-----------------------------------------------------------");

		unsigned int iIndividuo;

		printf("\n");	printf("iIndividuo");
		printf("\t");	printf("idx_Mae");
		printf("\t");	printf("p");
		printf("\t");	printf("f");
		printf("\t");	printf("pontos_de_corte[0]");
		printf("\t");	printf("pontos_de_corte[1]");

		struct testeGeracao_s host_testeGeracao;
		cudaErrorCheck(cudaMemcpy(&host_testeGeracao, dev_testeGeracao, sizeof(struct testeGeracao_s), cudaMemcpyDeviceToHost));

		for (iIndividuo = 0; iIndividuo < host_ParametrosGA->numIndividuos; iIndividuo++) {
			printf("\n");	printf("%d", iIndividuo);
			printf("\t");	printf("%d", host_testeGeracao.individuo[iIndividuo].idx_Mae);
			printf("\t");	printf("%f", host_testeGeracao.individuo[iIndividuo].p);
			printf("\t");	printf("%f", host_testeGeracao.individuo[iIndividuo].f);
			printf("\t");	printf("%d", host_testeGeracao.individuo[iIndividuo].pontos_de_corte[0]);
			printf("\t");	printf("%d", host_testeGeracao.individuo[iIndividuo].pontos_de_corte[1]);
		}

		printf("\n");	printf("Teste Crossover - Populacao Depois");
		printf("\n");	printf("-----------------------------------------------------------");
		imprimeGeracao(&host_Geracao, host_ParametrosGA);
	}
}

//===========================================================================================
void teste_Fitness_GPU(	struct generation *devGeracao,
								struct parametros *host_parametrosGA) {
//===========================================================================================

	struct generation host_Geracao;
	
	cudaErrorCheck(cudaMemcpy(&host_Geracao, devGeracao, sizeof(struct generation), cudaMemcpyDeviceToHost));
	
	printf("\n"); printf("--------------------------------------------------------------------"); 
	printf("\n"); printf("TESTE FITNESS - inicio"); 
	printf("\n"); printf("--------------------------------------------------------------------"); 

	imprimeGeracao(&host_Geracao, host_parametrosGA);

	printf("\n"); printf("--------------------------------------------------------------------"); 
	printf("\n"); printf("TESTE FITNESS - fim"); 
	printf("\n"); printf("--------------------------------------------------------------------"); 

}

//===========================================================================================
void testaHamiltoniano(	float *hostHamiltoniano,
								float *devHamiltoniano,
								struct parametros *hostParametrosGA) {
//===========================================================================================
	

	printf("\n"); printf("//*****************************************************************************************************************");
	printf("\n"); printf("// TESTE HAMILTONIANO - INICIO");
	printf("\n"); printf("//*****************************************************************************************************************");
	
	unsigned long int iElementoH, qtdeElementosH;

	// Primeiro passo: limpar o H da CPU
	printf("\n\nLimpando o H da CPU ... ");
	qtdeElementosH = hostParametrosGA->numGenes*hostParametrosGA->numGenes;
	for (iElementoH = 0; iElementoH < qtdeElementosH; iElementoH++) {
			*(hostHamiltoniano + iElementoH) = -1;
	}
	printf("Ok.");

	// Imprime o H da CPU para verificar se, de fato, está 'limpo'.
	printf("\nVerificando se, de fato, o H da CPU foi limpo: ");

	unsigned long int iLinha, iColuna;

	printf("\nH 'limpo' = \n\n");

	for (iLinha = 0; iLinha < hostParametrosGA->numGenes; iLinha++) {
		for (iColuna = 0; iColuna < hostParametrosGA->numGenes; iColuna++) {
			printf("%f\t", *(hostHamiltoniano + hostParametrosGA->numGenes*iLinha + iColuna));
		}
		printf("\n");
	}
	
	// Copia o Hamiltoniano da GPU para a CPU
	printf("\nCopiando o H da GPU para a CPU... ");
	cudaErrorCheck(cudaMemcpy(hostHamiltoniano, devHamiltoniano, hostParametrosGA->numGenes * hostParametrosGA->numGenes * sizeof(float), cudaMemcpyDeviceToHost));
	printf("Ok.");

	// Imprime o H que foi copiado da GPU
	printf("\nImprimindo o H copiado... ");
	
	printf("\nH copiado = \n\n");

	for (iLinha = 0; iLinha < hostParametrosGA->numGenes; iLinha++) {
		for (iColuna = 0; iColuna < hostParametrosGA->numGenes; iColuna++) {
			printf("%f\t", *(hostHamiltoniano + hostParametrosGA->numGenes*iLinha + iColuna));
		}
		printf("\n");
	}
	
	printf("\n"); printf("//*****************************************************************************************************************");
	printf("\n"); printf("// TESTE HAMILTONIANO - FIM");
	printf("\n"); printf("//*****************************************************************************************************************");

}

//---------------------------------------------------------------------------------------
unsigned short int imprimeComportamentoFitness_GPU(
							unsigned short int flag_Cabecalho,
							unsigned short int cod_Maquina,
							unsigned short int cod_tipoPrograma,
							unsigned short int cod_Tipo_Fitness,
							unsigned short int parm_iGeracao,
							struct generation *dev_Geracao,
							struct parametros_Metodo *dev_parametros_Metodo) {
//---------------------------------------------------------------------------------------

	if (flag_Cabecalho == 0) {
		printf("\n"); printf("Maquina");						// 0 - Adriano, 1 - FT
		printf("\t"); printf("Serial ou Paralelo?");		// 0 - Serial, 1 - Paralelo
		printf("\t"); printf("Quantidade Geracoes");		// constNumGeracoes
		printf("\t"); printf("Quantidade Individuos");	// constNumIndividuos
		printf("\t"); printf("Quantidade Genes");			// constNumGenes
		printf("\t"); printf("Tipo Fitness");
		printf("\t"); printf("Lambda");
		printf("\t"); printf("Rho minimo");
		printf("\t"); printf("iGeracao");				// parm_iGeracao
		printf("\t"); printf("Fitness Medio");
		printf("\t"); printf("Maior Fitness");
		printf("\t"); printf("Melhor rho");
		printf("\t"); printf("Posicao Melhor Individuo");
	}
	else {

		/*
			A cópia da geração e parâmetros do método da GPU para a CPU
			é a única diferença de código para o equivalente serial.
		*/

		struct generation host_Geracao;
		cudaErrorCheck( cudaMemcpy(&host_Geracao, dev_Geracao, sizeof(struct generation), cudaMemcpyDeviceToHost ) );

		struct parametros_Metodo host_parametros_Metodo;
		cudaErrorCheck( cudaMemcpy(&host_parametros_Metodo, dev_parametros_Metodo, sizeof(struct parametros_Metodo), cudaMemcpyDeviceToHost ) );

		char *strMaquina;
		switch (cod_Maquina)	{
			case 0: strMaquina = "00 - GT130M\0"; break;
			case 1: strMaquina = "01 - Beowulf\0"; break;
			default: {
					printf("\n");
					printf("Erro\timprimeTempo\tCodigo da maquina invalido.");
					return 1;
			} break;
		}
		
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

		char *str_Fitness;
		switch (cod_Tipo_Fitness) {
			case  0:	str_Fitness = "00 - Fitness Serial\0"; break;
			case  1:	str_Fitness = "01 - Fitness Paralelo 1d\0";	break;
			case  2:	str_Fitness = "02 - Fitness Paralelo 2d\0"; break;
			default: {
				printf("\n");
				printf("Erro\timprimeTempo\tCodigo invalido para o tipo de Fitness.");
				return 3;
			}; break;
		}


		printf("\n"); printf("%s", strMaquina);			
		printf("\t"); printf("%s", str_TipoPrograma);	
		printf("\t"); printf("%d", constNumGeracoes);	
		printf("\t"); printf("%d", constNumIndividuos);	
		printf("\t"); printf("%d", constNumGenes);		
		printf("\t"); printf("%s", str_Fitness);
		printf("\t"); printf("%f", host_parametros_Metodo.lambda);
		printf("\t"); printf("%f", host_parametros_Metodo.rho_minimo);
		printf("\t"); printf("%d", parm_iGeracao);
		printf("\t"); printf("%f", host_Geracao.FitnessMedio);
		printf("\t"); printf("%f", host_Geracao.Maior_fitness);
		printf("\t"); printf("%f", host_Geracao.Melhor_cociente_Rayleigh);
		printf("\t"); printf("%d", host_Geracao.idxMelhorIndividuo);

	} // final do if flag_Cabealho;

	return 0;
}