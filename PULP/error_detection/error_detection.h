#ifndef ERROR_DETECTION_H_
#define ERROR_DETECTION_H_

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../defines.h"

#define FACTOR_MS 1000
#define FACTOR_RATIO_RR 1000
#define PERCENTILE_LOO_LOW 641 //This is the 0.5th percentile of the RR distribution detected by the clustering using the leave-one-out approach for each subject (it is subject-specific)
#define PERCENTILE_LOO_HIGH 1492 //This is the 99.5th percentile of the RR distribution detected by the clustering using the leave-one-out approach for each subject (it is subject-specific)

int errorDetection(int32_t *arg[]);

#endif // ERROR_DETECTION_H_