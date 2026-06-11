
#ifndef __RELATIVEENERGY_H__
#define __RELATIVEENERGY_H__

#include <stdint.h>
#include "../defines.h"

//long window delay length (in samples): sampling frequency times long window length in seconds
//Long window length = 0.95 seconds
// #define LONG_WINDOW (uint16_t) (0.95*ECG_SAMPLING_FREQUENCY+1)

//RelEn parameters from Sasan's thesis
//Short widow = 140 ms
#define SHORT_WINDOW (uint16_t) (0.14*ECG_SAMPLING_FREQUENCY)
#define SHORT_WINDOW_HALF (uint16_t) (0.14/2*ECG_SAMPLING_FREQUENCY)
#define POWER (uint8_t) 2

//USE THIS ONE!
//Generate the relative energy raw signal used to calculate peaks
//Inputs: ecgBuffer - one long window of the signal. The first index of the buffer is the last value of the previous window
//		  lastShortEnergy - the calculated short energy of the last window
//		  lastLongEnergy - the calculated long energy of the last window
//Output: the relative energy coefficient corresponding to the buffer of values

//float getRelEnCoefficientsMoreOptimized(int16_t ecgBuffer[LONG_WINDOW + 1], int32_t* lastShortEnergy, int32_t* lastLongEnergy);

void relEn_w(int32_t *arg[]);
void clearAndResetRelEn();

#endif // __RELATIVEENERGY_H__
