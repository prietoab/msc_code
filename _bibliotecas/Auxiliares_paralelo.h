//---------------------------------------------------------------------------------------
__global__ void devInicializaCURAND(struct generation *devGeracaoInicial,
												unsigned int *globalSemente,
												curandState *estadosCURAND,
												struct parametros *dev_ParametrosGA) {
//---------------------------------------------------------------------------------------
	
	unsigned int iThread = threadIdx.x + blockIdx.x*blockDim.x;

	/*
		Como cada indiv�duo tem o seu pr�prio estado cuRAND, devemos
		inicializar todos. A semente � a mesma, conforme especifica��o
		do curand_init. Por�m, o curandState � �nico para cada indiv�duo,
		por meio da vari�vel estado_cuRAND da estrutura individual.
	*/
	unsigned int qtdeMaxThreads = dev_ParametrosGA->numIndividuos*dev_ParametrosGA->numGenes;
	if ( iThread < qtdeMaxThreads ) {
		curand_init(*globalSemente, iThread, 0, (estadosCURAND + iThread) );
	}
			// 0: semente. O mesmo 'seed' do rand() padr�o do C.
			// iIndividuo: identifica o id da sequencia. A mesma semente com diferentes n�meros de sequ�ncia
			//		levam a diferentes sequ�ncias de n�meros pseudorand�micos.
			// 0: offSet
}

//=======================================================================================
__device__ unsigned int dev_Get_Idx_CURAND(unsigned int iGene,
														 unsigned int numGenes,
														 unsigned int iIndividuo) {
//=======================================================================================
	unsigned int idx_CURAND = 0;
															 
	if (numGenes == 0) {
		idx_CURAND = iIndividuo;
	} 
	else {
		idx_CURAND = iGene + numGenes*iIndividuo;
	}

	return idx_CURAND;
}



//=================================================================================
__global__ void dev_copiaGeracoes (	struct generation *dev_GerOrigem,
												struct generation *dev_GerDestino) {
//=================================================================================
	unsigned int iThread = blockIdx.x * blockDim.x + threadIdx.x;

	if (iThread == 0) {
		*dev_GerDestino = *dev_GerDestino;
	}
}

//---------------------------------------------------------------------------------------
#define cudaErrorCheck(ans) { gpuAssert((ans), __FILE__,__LINE__);}
//---------------------------------------------------------------------------------------
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
	if (code != cudaSuccess) {
		fprintf(stderr, "GPUAssert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

//=======================================================================================
__device__ unsigned short int devRandomicoUI(unsigned short int valorMaximo, curandState *estadoCURAND) {								
//=======================================================================================
	return (unsigned short int)floor((float)valorMaximo*curand_uniform(estadoCURAND));	
};