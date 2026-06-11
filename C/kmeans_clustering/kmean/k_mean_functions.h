#include <stdint.h>
#include <stdio.h>
#include "../../defines.h"

#define PEAKS_STD 5
#define MAX_PEAKS_IN_WIND (N*15 + PEAKS_STD) 
#define MAX_QRS_DUR 35 //50, 35
#define MIN_RR_DIST 60 //50, 60

type_f ganlogfunc(type_f x, type_f k, type_f b, type_f nu, type_f m);
type_f gaussian(type_f x, type_f mu, type_f sign);
type_f percentile(int16_t *x, int end, int percToFind);
void partial_fit(type_f *x, type_f *hcentr, type_f *lcentr, int *h_el, int *l_el, int *label);
type_f std(type_i *x, int num_peaks, int num_peaks_total);

void rpeaks(int32_t* arg[]);
