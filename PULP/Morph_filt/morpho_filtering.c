
#include "morpho_filtering.h"
#include "rt/rt_api.h"
#include "../defines.h"

/////////////////////////////////////////////////////////////////////////////////
//-----------------------PRIVATE MEMORY OF CORE 0,1 and 2----------------------//
/////////////////////////////////////////////////////////////////////////////////
//El vector B1 definido en la tecnica se despliega en memoria
#if HIGH_FREQ==1
PULP_L1_DATA  int16_t B1[] = B1_STRUCTURING_ELEMENT;
#endif

//																		SIZE IN BYTES PER LEAD
//										 								(250 Hz)		(500 Hz)
PULP_L1_DATA int16_t bufferEB0[CHANNELS][TAMB0];					//		100				200
PULP_L1_DATA int16_t bufferDB0[CHANNELS][TAMB0];					//		100				200
PULP_L1_DATA int16_t bufferDBc[CHANNELS][TAMBc];					//		150				300
PULP_L1_DATA int16_t bufferEBc[CHANNELS][TAMBc];					//		150				300
PULP_L1_DATA int16_t bufferEnt[CHANNELS][TAMBEnt];					//		300				600
#if HIGH_FREQ==1
PULP_L1_DATA volatile int16_t bufferFBc1[CHANNELS][TAMBFB];			//		10				18
PULP_L1_DATA volatile int16_t bufferFBc1ed[CHANNELS][TAMBFB];		//		10				18
PULP_L1_DATA volatile int16_t bufferFBc1de[CHANNELS][TAMBFB];		//		10				18
#endif

PULP_L1_DATA int16_t bufferEntrada[CHANNELS][BUFFER_SIZE_FILTER];	//		2				2
PULP_L1_DATA int16_t bufferSalida[CHANNELS][BUFFER_SIZE_FILTER];	//		2				2
//volatile int16_t bufferSalida[CHANNELS];

PULP_L1_DATA int16_t samplesRead[CHANNELS];							//		2		 		2
PULP_L1_DATA int16_t actEntrada[CHANNELS];							//		2				2
PULP_L1_DATA int16_t contador[CHANNELS];							//		2				2
PULP_L1_DATA int16_t aBEB0[CHANNELS];								//		2				2
PULP_L1_DATA int16_t pMinEB0[CHANNELS];								//		2				2
PULP_L1_DATA int16_t p2MinEB0[CHANNELS];							//		2				2
PULP_L1_DATA int16_t aBDB0[CHANNELS];								//		2				2
PULP_L1_DATA int16_t pMaxDB0[CHANNELS];								//		2				2
PULP_L1_DATA int16_t p2MaxDB0[CHANNELS];							//		2				2
PULP_L1_DATA int16_t aBDBc[CHANNELS];								//		2				2
PULP_L1_DATA int16_t pMaxDBc[CHANNELS];								//		2				2
PULP_L1_DATA int16_t p2MaxDBc[CHANNELS];							//		2				2
PULP_L1_DATA int16_t aBEBc[CHANNELS];								//		2				2
PULP_L1_DATA int16_t pMinEBc[CHANNELS];								//		2				2
PULP_L1_DATA int16_t p2MinEBc[CHANNELS];							//		2				2
PULP_L1_DATA int16_t aBEnt[CHANNELS];								//		2				2
PULP_L1_DATA int16_t rBEnt[CHANNELS];								//		2				2

#if HIGH_FREQ==1
PULP_L1_DATA volatile int16_t aBFBc1[CHANNELS];						//		2				2
PULP_L1_DATA volatile int16_t aFBc1x[CHANNELS];						//		2				2
PULP_L1_DATA volatile int16_t pMaxFBc1ed[CHANNELS];					//		2				2
PULP_L1_DATA volatile int16_t pMinFBc1de[CHANNELS];					//		2				2
PULP_L1_DATA volatile int16_t p2MinFBc1de[CHANNELS];				//		2				2
PULP_L1_DATA volatile int16_t p2MaxFBc1ed[CHANNELS];				//		2				2
PULP_L1_DATA volatile int16_t indiceB1[CHANNELS];					//		2				2
PULP_L1_DATA volatile int16_t min1E[CHANNELS];						//		2				2
PULP_L1_DATA volatile int16_t max1D[CHANNELS];						//		2				2
#endif

PULP_L1_DATA int16_t vaux[CHANNELS];								//		2				2
PULP_L1_DATA int16_t v2aux[CHANNELS];								//		2				2

																	//TOTAL	890				1714

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


void fill_in_filter(int16_t input, int16_t index) {

	actEntrada[index] = 0;

	bufferEntrada[index][0] = input;

	samplesRead[index] = BUFFER_SIZE_FILTER;
	if (contador[index] < (TAMB0 * 3) + (TAMB0 / 2) - (TAMBcm - TAMB0m) - 3)
		(contador[index])++;
}

void storeOutput(int16_t v, int16_t index) {
	bufferSalida[index][0] = v;
}

//Opening: dilation(erosion(X))
//Closing: erosion(dilation(X))

//Definicion de buffers segun su aparicion en la tecnica:
//bufferEB0:	opening(entrada) -> erosion(entrada, B0)
//bufferDB0:	opening(entrada) -> dilation(erosion(entrada, B0))
//bufferDBc:	closing(fbt)     -> dilation(fbt, Bc)
//bufferEBc:	closing(fbt)     -> erosion(dilation(fbt, Bc))
//bufferEnt:    buffer necesario para superar la diferencia de tamaÃ±o entre
//              los buffers del final del calculo y los iniciales y poder
//				realizar la resta fbc = entrada - fb

// IF  HIGH_FREQ==1:
//BufferFBc1:	buffer con los valores fbc para los siguientes
//BufferFBc1e:  opening(fbc)	 -> erosion(fbc, B1);
//BufferFBc1ed: opening(fbc)	 -> dilation(erosion(fbc, B1));
//BufferFBc1d:  closing(fbc)	 -> dilation(fbc, B1);
//BufferFBc1de: closing(fbc)	 -> erosion(dilation(fbc,B1));
// END IF

//A partir de los valores de BufferFB2de y BufferFB2ed es posible obtener
//el valor resultado del filtrado.


void init_filtering() {

	for(int16_t i = 0; i < CHANNELS; i++) {

	#if HIGH_FREQ==1
		aBFBc1[i]      = 0;
		aFBc1x[i]      = 0;
		pMaxFBc1ed[i]  = 0;
		pMinFBc1de[i]  = 0;
		p2MinFBc1de[i] = -1;
		p2MaxFBc1ed[i] = -1;
		indiceB1[i]    = 0;
	#endif

		contador[i] = 0;
		aBEB0[i]    = TAMB0 - 1;
		pMinEB0[i]  = 0;
		p2MinEB0[i] = -1;
		aBDB0[i]    = TAMB0m;
		pMaxDB0[i]  = 0;
		
		p2MaxDB0[i] = -1;
		aBDBc[i]    = TAMB0m;
		pMaxDBc[i]  = 0;
		p2MaxDBc[i] = -1;
		aBEBc[i]    = TAMB0m;
		pMinEBc[i]  = 0;
		p2MinEBc[i] = -1;
		aBEnt[i]    = 0;
		rBEnt[i]    = -1;
	}
}


void initBuffers(int16_t i, int16_t input_value, int16_t index) {

	#if HIGH_FREQ==1
	int16_t j;
	int16_t component;
	#endif
 
	actEntrada[index] = 0;

	//VICTOR: se incrementa el contador, al ser llamado este procedimiento para cada nuevo sample.
	(contador[index])++;
	samplesRead[index] = 1;
 
	//Llenado inicial de los buffers
	//FRAN:morph --> change bufferEntrada[i] by 0
 	
	bufferEB0[index][i] = input_value;
	if (input_value <= bufferEB0[index][pMinEB0[index]])
  	{
		p2MinEB0[index] = pMinEB0[index];
		pMinEB0[index] = i;
	}
	else if (p2MinEB0[index] != -1 && input_value <= bufferEB0[index][p2MinEB0[index]])
  	{
		p2MinEB0[index] = i;
	}
 
	if (i >= TAMB0m)
		bufferEnt[index][(aBEnt[index])++] = input_value;
	else
  	{
  	
		bufferDB0[index][i] = input_value;
		if (input_value >= bufferDB0[index][pMaxDB0[index]])
    		{
			p2MaxDB0[index] = pMaxDB0[index];
			pMaxDB0[index] = i;
		}
		else if (input_value >= bufferDB0[index][p2MaxDB0[index]])
			p2MaxDB0[index] = i;

		bufferDBc[index][i] = input_value;
		if (input_value >= bufferDBc[index][pMaxDBc[index]])
    		{
			p2MaxDBc[index] = pMaxDBc[index];
			pMaxDBc[index] = i;
		}
		else if (input_value >= bufferDBc[index][p2MaxDBc[index]])
			p2MaxDBc[index] = i;

		bufferEBc[index][i] = input_value;
		if (input_value <= bufferEBc[index][pMinEBc[index]])
    		{
			if (pMinEBc[index] != 0)
				p2MinEBc[index] = pMinEBc[index];
			pMinEBc[index] = i;
		}
		else if (p2MinEBc[index] != -1 && input_value <= bufferEBc[index][p2MinEBc[index]])
			p2MinEBc[index] = i;
 

		#if HIGH_FREQ==1
		bufferFBc1[index][aBFBc1[index]] = 0;
		
		if (aBFBc1[index] != TAMBFB - 1)
			(aBFBc1[index])++;
		else
			aBFBc1[index] = 0;

		if (i >= TAMBFB - 1) //i >= 4
    		{
							   // A partir del bufferFBc1, se calcula erosion y dilation
							   // y se almacenan los valores sobre los buffers FBc1d y FBc1e
							   // cuyos limites se calculan en paralelo.
			j = indiceB1[index];
			if (indiceB1[index] != TAMB1*(TAMB1 - 1))
				indiceB1[index] += TAMB1;
			else
				indiceB1[index] = 0;

			min1E[index] = bufferFBc1[index][0] - B1[j];
			max1D[index] = bufferFBc1[index][0] + B1[j];

			for (component = 1; component < TAMB1; component++)
			{
				min1E[index] = min1E[index]>bufferFBc1[index][component] - B1[j + component] ? bufferFBc1[index][component] - B1[j + component] : min1E[index];
				max1D[index] = max1D[index]<bufferFBc1[index][component] + B1[j + component] ? bufferFBc1[index][component] + B1[j + component] : max1D[index];
			}

			if (i >= TAMBFB - 1 + TAMBFB - 1)	// i >= 8
      			{
												//Paso de valores y mover puntero de bufferFBc1ed
				if (aBFBc1[index] == pMaxFBc1ed[index]) {
					if (min1E[index] < bufferFBc1ed[index][pMaxFBc1ed[index]]) {
						if (p2MaxFBc1ed[index] == -1) {
							bufferFBc1ed[index][aFBc1x[index]] = min1E[index];

							v2aux[index] = INT16_MIN;
							vaux[index] = bufferFBc1ed[index][TAMBFB - 1];
							pMaxFBc1ed[index] = TAMBFB - 1;
							for (j = TAMBFB - 2; j >= 0; j--) {
								if (bufferFBc1ed[index][j] > vaux[index]) {
									v2aux[index] = vaux[index];
									vaux[index] = bufferFBc1ed[index][j];
									p2MaxFBc1ed[index] = pMaxFBc1ed[index];
									pMaxFBc1ed[index] = j;
								}
								else if (bufferFBc1ed[index][j] > v2aux[index]) {
									v2aux[index] = bufferFBc1ed[index][j];
									p2MaxFBc1ed[index] = j;
								}
							}
						}
						else if (min1E[index] < bufferFBc1ed[index][p2MaxFBc1ed[index]]) {
							pMaxFBc1ed[index] = p2MaxFBc1ed[index];
							p2MaxFBc1ed[index] = -1;
						}
					}
				}
				else if (min1E[index] >= bufferFBc1ed[index][pMaxFBc1ed[index]]) {
					p2MaxFBc1ed[index] = pMaxFBc1ed[index];
					pMaxFBc1ed[index] = aFBc1x[index];
				}
				else if (aFBc1x[index] == p2MaxFBc1ed[index]) {
					if (min1E[index] < bufferFBc1ed[index][p2MaxFBc1ed[index]])
						p2MaxFBc1ed[index] = -1;
				}
				else if (p2MaxFBc1ed[index] != -1 && min1E[index] >= bufferFBc1ed[index][p2MaxFBc1ed[index]])
					p2MaxFBc1ed[index] = aFBc1x[index];

				bufferFBc1ed[index][aFBc1x[index]] = min1E[index];

				if (aFBc1x[index] == pMinFBc1de[index]) {
					if (max1D[index] > bufferFBc1de[index][pMinFBc1de[index]]) {
						if (p2MinFBc1de[index] == -1) {
							bufferFBc1de[index][aFBc1x[index]] = max1D[index];

							v2aux[index] = INT16_MAX;
							vaux[index] = bufferFBc1de[index][TAMBFB - 1];
							pMinFBc1de[index] = TAMBFB - 1;
							for (j = TAMBFB - 2; j >= 0; j--) {
								if (bufferFBc1de[index][j] < vaux[index]) {
									v2aux[index] = vaux[index];
									vaux[index] = bufferFBc1de[index][j];
									p2MinFBc1de[index] = pMinFBc1de[index];
									pMinFBc1de[index] = j;
								}
								else if (bufferFBc1de[index][j] < v2aux[index]) {
									v2aux[index] = bufferFBc1de[index][j];
									p2MinFBc1de[index] = j;
								}
							}
						}
						else if (max1D[index] > bufferFBc1de[index][p2MinFBc1de[index]]) {
							pMinFBc1de[index] = p2MinFBc1de[index];
							p2MinFBc1de[index] = -1;
						}
					}
				}
				else if (max1D[index] <= bufferFBc1de[index][pMinFBc1de[index]]) {
					p2MinFBc1de[index] = pMinFBc1de[index];
					pMinFBc1de[index] = aFBc1x[index];
				}
				else if (aFBc1x[index] == p2MinFBc1de[index]) {
					if (max1D[index] > bufferFBc1de[index][p2MinFBc1de[index]])
						p2MinFBc1de[index] = -1;
				}
				else if (p2MinFBc1de[index] != -1 && max1D[index] <= bufferFBc1de[index][p2MinFBc1de[index]])
					p2MinFBc1de[index] = aFBc1x[index];

				bufferFBc1de[index][aFBc1x[index]] = max1D[index];
			}

			if (aFBc1x[index] != TAMBFB - 1)
				(aFBc1x[index])++;
			else
				aFBc1x[index] = 0;
		}
		#endif //HIGH_FREQ==1
	}
}


void proceso(int16_t index) {

	int16_t i; 
	#if HIGH_FREQ==1
	int16_t j, component;
	#endif

	while (actEntrada[index] < samplesRead[index]) {
		//printf("core %d. %d:%d\n", rt_core_id(), actEntrada[index], samplesRead[index]);
 
		if (aBEB0[index] == pMinEB0[index]) {
			if (bufferEntrada[index][actEntrada[index]] > bufferEB0[index][pMinEB0[index]]) {
				if (p2MinEB0[index] == -1) {
					bufferEB0[index][aBEB0[index]] = bufferEntrada[index][actEntrada[index]];
					v2aux[index] = INT16_MAX;
					vaux[index] = bufferEB0[index][TAMB0 - 1];
					pMinEB0[index] = TAMB0 - 1;
					for (i = TAMB0 - 2; i >= 0; i--) {
						if (bufferEB0[index][i] < vaux[index]) {
							v2aux[index] = vaux[index];
							vaux[index] = bufferEB0[index][i];
							p2MinEB0[index] = pMinEB0[index];
							pMinEB0[index] = i;
						}
						else if (bufferEB0[index][i] < v2aux[index]) {
							v2aux[index] = bufferEB0[index][i];
							p2MinEB0[index]  = i;
						}
					}
					// return 0;
				}
				else if (bufferEntrada[index][actEntrada[index]] > bufferEB0[index][p2MinEB0[index]]) {
					pMinEB0[index] = p2MinEB0[index] ;
					p2MinEB0[index]  = -1;
				}
			}
		}
		else if (bufferEntrada[index][actEntrada[index]] <= bufferEB0[index][pMinEB0[index]]) {
			p2MinEB0[index]  = pMinEB0[index];
			pMinEB0[index] = aBEB0[index];
		}
		else if (aBEB0[index] == p2MinEB0[index] ) {
			if (bufferEntrada[index][actEntrada[index]] > bufferEB0[index][p2MinEB0[index]]) {
				p2MinEB0[index]  = -1;
			}
		}
		else if (p2MinEB0[index]  != -1 && bufferEntrada[index][actEntrada[index]] <= bufferEB0[index][p2MinEB0[index]])
		{
			p2MinEB0[index]  = aBEB0[index];
		}

		bufferEnt[index][aBEnt[index]] = bufferEntrada[index][actEntrada[index]];

		if (aBEnt[index] != TAMBEnt - 1) {
			(aBEnt[index])++;
		} else {
			aBEnt[index] = 0;
		}

		bufferEB0[index][aBEB0[index]] = bufferEntrada[index][(actEntrada[index])++];

		if (aBEB0[index] != TAMB0 - 1) {
			(aBEB0[index])++;
		} else {
			(aBEB0[index]) = 0;
		}

		if (aBDB0[index]  == pMaxDB0[index] ) {
			if (bufferEB0[index][pMinEB0[index]] < bufferDB0[index][pMaxDB0[index]]) {
				if (p2MaxDB0[index] == -1) {
					bufferDB0[index][aBDB0[index]] = bufferEB0[index][pMinEB0[index]];
					v2aux[index] = INT16_MIN;
					vaux[index] = bufferDB0[index][TAMB0 - 1];
					pMaxDB0[index]  = TAMB0 - 1;
					for (i = TAMB0 - 2; i >= 0; i--) {
						if (bufferDB0[index][i] > vaux[index]) {
							v2aux[index] = vaux[index];
							vaux[index] = bufferDB0[index][i];
							p2MaxDB0[index]  = pMaxDB0[index] ;
							pMaxDB0[index]  = i;
						}
						else if (bufferDB0[index][i] > v2aux[index]) {
							v2aux[index] = bufferDB0[index][i];
							p2MaxDB0[index]  = i;
						}
					}
				// return 0;
				}
				else if (bufferEB0[index][pMinEB0[index]] < bufferDB0[index][p2MaxDB0[index]]) {
					pMaxDB0[index]  = p2MaxDB0[index] ;
					p2MaxDB0[index]  = -1;
				}
			}
		}
		else if (bufferEB0[index][pMinEB0[index]] >= bufferDB0[index][pMaxDB0[index]]) {
			p2MaxDB0[index]  = pMaxDB0[index] ;
			pMaxDB0[index]  = aBDB0[index] ;
		}
		else if (aBDB0[index]  == p2MaxDB0[index] ) {
			if (bufferEB0[index][pMinEB0[index]] < bufferDB0[index][p2MaxDB0[index]]) {
				p2MaxDB0[index]  = -1;
			}
		}
		else if (p2MaxDB0[index]  != -1 && bufferEB0[index][pMinEB0[index]] >= bufferDB0[index][p2MaxDB0[index]])
		{
			p2MaxDB0[index]  = aBDB0[index];
		}

		bufferDB0[index][aBDB0[index]] = bufferEB0[index][pMinEB0[index]];

		if (aBDB0[index]  != TAMB0 - 1) {
			(aBDB0[index] )++;
		} else {
			aBDB0[index]  = 0;
		}

		//Paso de datos al buffer DBc si se ha llenado el buffer DB0.
		//			if (!firstPass[k] || actEntrada[k] >= TAMB0 + (TAMB0/2) - 1) { // 74
		//VICTOR: he cambiado esto
		if (contador[index] >= TAMB0 + (TAMB0 / 2) - 1) {
		 
			if (aBDBc[index] == pMaxDBc[index]) {
				if (bufferDB0[index][pMaxDB0[index]] < bufferDBc[index][pMaxDBc[index]]) {
					if (p2MaxDBc[index] == -1) {
						bufferDBc[index][aBDBc[index] ] = bufferDB0[index][pMaxDB0[index] ];
						v2aux[index] = INT16_MIN;
						vaux[index] = bufferDBc[index][TAMBc - 1];
						pMaxDBc[index] = TAMBc - 1;
						for (i = TAMBc - 2; i >= 0; i--) {
							if (bufferDBc[index][i] > vaux[index]) {
								v2aux[index] = vaux[index];
								vaux[index] = bufferDBc[index][i];
								p2MaxDBc[index] = pMaxDBc[index];
								pMaxDBc[index] = i;
							}
							else if (bufferDBc[index][i] > v2aux[index]) {
								v2aux[index] = bufferDBc[index][i];
								p2MaxDBc[index] = i;
							}
						}
					// return 0;
					}
					else if (bufferDB0[index][pMaxDB0[index]] < bufferDBc[index][p2MaxDBc[index]]) {
						pMaxDBc[index] = p2MaxDBc[index];
						p2MaxDBc[index] = -1;
					}
				}
			}
			else if (bufferDB0[index][pMaxDB0[index]] >= bufferDBc[index][pMaxDBc[index]]) {
				p2MaxDBc[index] = pMaxDBc[index];
				pMaxDBc[index] = aBDBc[index] ;
			}
			else if (aBDBc[index]  == p2MaxDBc[index]) {
				if (bufferDB0[index][pMaxDB0[index]] < bufferDBc[index][p2MaxDBc[index]]) {
					p2MaxDBc[index] = -1;
				}
			}
			else if (p2MaxDBc[index] != -1 && bufferDB0[index][pMaxDB0[index]] >= bufferDBc[index][p2MaxDBc[index]])
			{
				p2MaxDBc[index] = aBDBc[index];
			}

			bufferDBc[index][aBDBc[index]] = bufferDB0[index][pMaxDB0[index]];

			if (aBDBc[index]  != TAMBc - 1) {
				(aBDBc[index] )++;
			} else {
				aBDBc[index]  = 0;
			}
 
			//Paso de datos al buffer EBc si se ha llenado DBc
			//if (!firstPass[k] || actEntrada[k] >= (TAMB0*2) + (TAMB0/2) - 2) { //123
			//VICTOR: he cambiado esto.
			if (contador[index] >= (TAMB0 * 2) + (TAMB0 / 2) - 2)
			{
				if (aBEBc[index] == pMinEBc[index]) {
					if (bufferDBc[index][pMaxDBc[index]] > bufferEBc[index][pMinEBc[index]]) {
						if (p2MinEBc[index] == -1) {
							bufferEBc[index][aBEBc[index]] = bufferDBc[index][pMaxDBc[index]];
							v2aux[index] = INT16_MAX;
							vaux[index] = bufferEBc[index][TAMBc - 1];
							pMinEBc[index] = TAMBc - 1;
							for (i = TAMBc - 2; i >= 0; i--) {
								if (bufferEBc[index][i] < vaux[index]) {
									v2aux[index] = vaux[index];
									vaux[index] = bufferEBc[index][i];
									p2MinEBc[index] = pMinEBc[index];
									pMinEBc[index] = i;
								}
								else if (bufferEBc[index][i] < v2aux[index]) {
									v2aux[index] = bufferEBc[index][i];
									p2MinEBc[index] = i;
								}
							}
							// return 0;
						}
						else if (bufferDBc[index][pMaxDBc[index]] > bufferEBc[index][p2MinEBc[index]]) {
							pMinEBc[index] = p2MinEBc[index];
							p2MinEBc[index] = -1;
						}
					}
				}
				else if (bufferDBc[index][pMaxDBc[index]] <= bufferEBc[index][pMinEBc[index]]) {
					p2MinEBc[index] = pMinEBc[index];
					pMinEBc[index] = aBEBc[index];
				}
				else if (aBEBc[index] == p2MinEBc[index]) {
					if (bufferDBc[index][pMaxDBc[index]] > bufferEBc[index][p2MinEBc[index]]) {
						p2MinEBc[index] = -1;
					}
				}
				else if (p2MinEBc[index] != -1 && bufferDBc[index][pMaxDBc[index]] <= bufferEBc[index][p2MinEBc[index]])
				{
					p2MinEBc[index] = aBEBc[index];
				}

				bufferEBc[index][aBEBc[index]] = bufferDBc[index][pMaxDBc[index]];

				if (aBEBc[index] != TAMBc - 1) {
					(aBEBc[index])++;
				} else {
					aBEBc[index] = 0;
				}


				//Paso de datos al buffer FBc1 si se ha llenado EBc
				//if (!firstPass[k] || actEntrada[k] >= (TAMB0*3) + (TAMB0/2) - (TAMBcm - TAMB0m) - 3) {
				//VICTOR: he cambiado esto.

				if (contador[index] >= (TAMB0 * 3) + (TAMB0 / 2) - (TAMBcm - TAMB0m) - 3) {
					
					#if HIGH_FREQ==1
					bufferFBc1[index][aBFBc1[index]] = bufferEnt[index][rBEnt[index]] - bufferEBc[index][pMinEBc[index]];

					if (aBFBc1[index] != TAMBFB - 1)
						(aBFBc1[index])++;
					else
						aBFBc1[index] = 0;
					#else
					storeOutput(bufferEnt[index][rBEnt[index]] - bufferEBc[index][pMinEBc[index]], index);
					#endif

					if (rBEnt[index] != TAMBEnt - 1)
						(rBEnt[index])++;
					else
						rBEnt[index] = 0;


					#if HIGH_FREQ==1
					j = indiceB1[index];
					if (indiceB1[index] != TAMB1*(TAMB1 - 1))
						indiceB1[index] += TAMB1;
					else
						indiceB1[index] = 0;

					min1E[index] = bufferFBc1[index][0] - B1[j];
					max1D[index] = bufferFBc1[index][0] + B1[j];

					for (component = 1; component < TAMB1; component++)
					{
						if (min1E[index]>bufferFBc1[index][component] - B1[j + component])
							min1E[index] = bufferFBc1[index][component] - B1[j + component];
						if (max1D[index]<bufferFBc1[index][component] + B1[j + component])
							max1D[index] = bufferFBc1[index][component] + B1[j + component];
					}

					//Paso de valores y mover puntero de bufferFBc1ed
					if (aFBc1x[index] == pMaxFBc1ed[index]) {
						if (min1E[index] < bufferFBc1ed[index][pMaxFBc1ed[index]]) {
							if (p2MaxFBc1ed[index] == -1) {
								bufferFBc1ed[index][aFBc1x[index]] = min1E[index];

								v2aux[index] = INT16_MIN;
								vaux[index] = bufferFBc1ed[index][TAMBFB - 1];
								pMaxFBc1ed[index] = TAMBFB - 1;
								for (j = TAMBFB - 2; j >= 0; j--) {
									if (bufferFBc1ed[index][j] > vaux[index]) {
										v2aux[index] = vaux[index];
										vaux[index] = bufferFBc1ed[index][j];
										p2MaxFBc1ed[index] = pMaxFBc1ed[index];
										pMaxFBc1ed[index] = j;
									}
									else if (bufferFBc1ed[index][j] > v2aux[index]) {
										v2aux[index] = bufferFBc1ed[index][j];
										p2MaxFBc1ed[index] = j;
									}
								}
							}
							else if (min1E[index] < bufferFBc1ed[index][p2MaxFBc1ed[index]]) {
								pMaxFBc1ed[index] = p2MaxFBc1ed[index];
								p2MaxFBc1ed[index] = -1;
							}
						}
					}
					else if (min1E[index] >= bufferFBc1ed[index][pMaxFBc1ed[index]]) {
						p2MaxFBc1ed[index] = pMaxFBc1ed[index];
						pMaxFBc1ed[index] = aFBc1x[index];
					}
					else if (aFBc1x[index] == p2MaxFBc1ed[index]) {
						if (min1E[index] < bufferFBc1ed[index][p2MaxFBc1ed[index]])
							p2MaxFBc1ed[index] = -1;
					}
					else if (p2MaxFBc1ed[index] != -1 && min1E[index] >= bufferFBc1ed[index][p2MaxFBc1ed[index]])
					{
						p2MaxFBc1ed[index] = aFBc1x[index];
					}

					bufferFBc1ed[index][aFBc1x[index]] = min1E[index];

					if (aFBc1x[index] == pMinFBc1de[index]) {
						if (max1D[index] > bufferFBc1de[index][pMinFBc1de[index]]) {
							if (p2MinFBc1de[index] == -1) {
								bufferFBc1de[index][aFBc1x[index]] = max1D[index];

								v2aux[index] = INT16_MAX;
								vaux[index] = bufferFBc1de[index][TAMBFB - 1];
								pMinFBc1de[index] = TAMBFB - 1;
								for (j = TAMBFB - 2; j >= 0; j--) {
									if (bufferFBc1de[index][j] < vaux[index]) {
										v2aux[index] = vaux[index];
										vaux[index] = bufferFBc1de[index][j];
										p2MinFBc1de[index] = pMinFBc1de[index];
										pMinFBc1de[index] = j;
									}
									else if (bufferFBc1de[index][j] < v2aux[index]) {
										v2aux[index] = bufferFBc1de[index][j];
										p2MinFBc1de[index] = j;
									}
								}
							}
							else if (max1D[index] > bufferFBc1de[index][p2MinFBc1de[index]]) {
								pMinFBc1de[index] = p2MinFBc1de[index];
								p2MinFBc1de[index] = -1;
							}
						}
					}
					else if (max1D[index] <= bufferFBc1de[index][pMinFBc1de[index]]) {
						p2MinFBc1de[index] = pMinFBc1de[index];
						pMinFBc1de[index] = aFBc1x[index];
					}
					else if (aFBc1x[index] == p2MinFBc1de[index]) {
						if (max1D[index] > bufferFBc1de[index][p2MinFBc1de[index]])
							p2MinFBc1de[index] = -1;
					}
					else if (p2MinFBc1de[index] != -1 && max1D[index] <= bufferFBc1de[index][p2MinFBc1de[index]])
						p2MinFBc1de[index] = aFBc1x[index];

					bufferFBc1de[index][aFBc1x[index]] = max1D[index];

					if (aFBc1x[index] != TAMBFB - 1)
						(aFBc1x[index])++;
					else
						aFBc1x[index] = 0;

					storeOutput((bufferFBc1ed[index][pMaxFBc1ed[index]] + bufferFBc1de[index][pMinFBc1de[index]]) >> 1, index);
					#endif //HIGH_FREQ==1
				}
	
				//Si no esta lleno el bufferDBc
				//Llenamos el bufferEBc con los valores de DBc en la diferencia
			}
			else if (aBDBc[index]  <= TAMBcm)
			{
				storeOutput(bufferEnt[index][rBEnt[index]] - bufferDBc[index][aBDBc[index]  - 1], index);

				#if HIGH_FREQ==1
				bufferFBc1[index][aBFBc1[index]] = bufferEnt[index][rBEnt[index]] - bufferDBc[index][aBDBc[index]  - 1];
				#else
				storeOutput(bufferEnt[index][rBEnt[index]] - bufferDBc[index][aBDBc[index]  - 1], index);
				#endif

				if (rBEnt[index] != TAMBEnt - 1)
					(rBEnt[index])++;
				else
					rBEnt[index] = 0;

				
				#if HIGH_FREQ==1
				if (aBFBc1[index] != TAMBFB - 1)
					(aBFBc1[index])++;
				else
					aBFBc1[index] = 0;

				if (bufferDBc[index][aBDBc[index]  - 1] <= bufferEBc[index][pMinEBc[index]]) {
					p2MinEBc[index] = pMinEBc[index];
					pMinEBc[index] = aBEBc[index];
				}
				else if (p2MinEBc[index] != -1 && bufferDBc[index][aBDBc[index]  - 1] <= bufferEBc[index][p2MinEBc[index]])
					p2MinEBc[index] = aBEBc[index];

				bufferEBc[index][aBEBc[index]] = bufferDBc[index][aBDBc[index]  - 1];
				(aBEBc[index])++;

				j = indiceB1[index];
				if (indiceB1[index] != TAMB1*(TAMB1 - 1))
					indiceB1[index] += TAMB1;
				else
					indiceB1[index] = 0;

				min1E[index] = bufferFBc1[index][0] - B1[j];
				max1D[index] = bufferFBc1[index][0] + B1[j];

				for (component = 1; component<TAMB1; component++)
				{

					min1E[index] = min1E[index]>bufferFBc1[index][component] - B1[j + component] ? bufferFBc1[index][component] - B1[j + component] : min1E[index];
					max1D[index] = max1D[index]<bufferFBc1[index][component] + B1[j + component] ? bufferFBc1[index][component] + B1[j + component] : max1D[index];
				}


				//Paso de valores y mover puntero de bufferFBc1ed
				if (aFBc1x[index] == pMaxFBc1ed[index]) {
					if (min1E[index] < bufferFBc1ed[index][pMaxFBc1ed[index]]) {
						if (p2MaxFBc1ed[index] == -1) {
							bufferFBc1ed[index][aFBc1x[index]] = min1E[index];

							//vaux = bufferFBc1ed[0][k];
							v2aux[index] = INT16_MIN;
							//pMaxFBc1ed[index][k] = 0;
							//for (j = 1; j < TAMBFB; j++) {
							vaux[index] = bufferFBc1ed[index][TAMBFB - 1];
							pMaxFBc1ed[index] = TAMBFB - 1;
							for (j = TAMBFB - 2; j >= 0; j--) {
								if (bufferFBc1ed[index][j] > vaux[index]) {
									v2aux[index] = vaux[index];
									vaux[index] = bufferFBc1ed[index][j];
									p2MaxFBc1ed[index] = pMaxFBc1ed[index];
									pMaxFBc1ed[index] = j;
								}
								else if (bufferFBc1ed[index][j] > v2aux[index]) {
									v2aux[index] = bufferFBc1ed[index][j];
									p2MaxFBc1ed[index] = j;
								}
							}
						}
						else if (min1E[index] < bufferFBc1ed[index][p2MaxFBc1ed[index]]) {
							pMaxFBc1ed[index] = p2MaxFBc1ed[index];
							p2MaxFBc1ed[index] = -1;
						}
					}
				}
				else if (min1E[index] >= bufferFBc1ed[index][pMaxFBc1ed[index]]) {
					p2MaxFBc1ed[index] = pMaxFBc1ed[index];
					pMaxFBc1ed[index] = aFBc1x[index];
				}
				else if (aFBc1x[index] == p2MaxFBc1ed[index]) {
					if (min1E[index] < bufferFBc1ed[index][p2MaxFBc1ed[index]])
						p2MaxFBc1ed[index] = -1;
				}
				else if (p2MaxFBc1ed[index] != -1 && min1E[index] >= bufferFBc1ed[index][p2MaxFBc1ed[index]])
					p2MaxFBc1ed[index] = aFBc1x[index];

				bufferFBc1ed[index][aFBc1x[index]] = min1E[index];

				if (aFBc1x[index] == pMinFBc1de[index]) {
					if (max1D[index] > bufferFBc1de[index][pMinFBc1de[index]]) {
						if (p2MinFBc1de[index] == -1) {
							bufferFBc1de[index][aFBc1x[index]] = max1D[index];

							v2aux[index] = INT16_MAX;
							vaux[index] = bufferFBc1de[index][TAMBFB - 1];
							pMinFBc1de[index] = TAMBFB - 1;
							for (j = TAMBFB - 2; j >= 0; j--) {
								if (bufferFBc1de[index][j] < vaux[index]) {
									v2aux[index] = vaux[index];
									vaux[index] = bufferFBc1de[index][j];
									p2MinFBc1de[index] = pMinFBc1de[index];
									pMinFBc1de[index] = j;
								}
								else if (bufferFBc1de[index][j] < v2aux[index]) {
									v2aux[index] = bufferFBc1de[index][j];
									p2MinFBc1de[index] = j;
								}
							}
						}
						else if (max1D[index] > bufferFBc1de[index][p2MinFBc1de[index]]) {
							pMinFBc1de[index] = p2MinFBc1de[index];
							p2MinFBc1de[index] = -1;
						}
					}
				}
				else if (max1D[index] <= bufferFBc1de[index][pMinFBc1de[index]]) {
					p2MinFBc1de[index] = pMinFBc1de[index];
					pMinFBc1de[index] = aFBc1x[index];
				}
				else if (aFBc1x[index] == p2MinFBc1de[index]) {
					if (max1D[index] > bufferFBc1de[index][p2MinFBc1de[index]])
						p2MinFBc1de[index] = -1;
				}
				else if (p2MinFBc1de[index] != -1 && max1D[index] <= bufferFBc1de[index][p2MinFBc1de[index]])
					p2MinFBc1de[index] = aFBc1x[index];

				bufferFBc1de[index][aFBc1x[index]] = max1D[index];

				if (aFBc1x[index] != TAMBFB - 1)
					(aFBc1x[index])++;
				else
					aFBc1x[index] = 0;

				storeOutput((bufferFBc1ed[index][pMaxFBc1ed[index]] + bufferFBc1de[index][pMinFBc1de[index]]) >> 1, index);
				#endif //HIGH_FREQ==1
			}
		}
	}//while
}

int16_t filterSample(int16_t offset, int16_t read_data, int16_t index) {

	if (offset < TAMB0) {
		initBuffers(offset, read_data, index);
	} else {
		fill_in_filter(read_data, index);
		proceso(index); 
	}

	return bufferSalida[index][0]; 
}

void filterWindows(int32_t *arg[])
{
	int16_t off;
	int16_t *ecg_buffer = (int16_t*) arg[0];
	int32_t *flag = arg[1];
	int16_t *i_lead = (int16_t*) arg[2];
	int32_t *bufferSize = arg[3];

	if(*flag == 0){
		off = 0;
	} else {
		off = TAMB0;
	}

	for(int indsample = 0; indsample<*bufferSize; indsample++) {
		ecg_buffer[indsample + (dim+LONG_WINDOW) * (*i_lead)] = filterSample(off, ecg_buffer[indsample + (dim+LONG_WINDOW) * (*i_lead)], (*i_lead));
		off++;
	}
}
