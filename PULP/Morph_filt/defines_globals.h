#ifndef DEFINES_GLOBALS_H
#define DEFINES_GLOBALS_H

#include "../defines.h"

/************** APPLICATION PARAMETERS *************/
#define SAMPLING_FREQUENCY ECG_SAMPLING_FREQUENCY		//This define is compulsory and must be either 250 or 500	ma
#define HIGH_FREQ 1				//1: do high-freq filtering and baseline removal, 0: do baseline removal only
/***************************************************/

#define CHANNELS NLEADS

#include <stdint.h>

#endif
