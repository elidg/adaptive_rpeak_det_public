#ifndef ERROR_DETECTION_H_
#define ERROR_DETECTION_H_

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../defines.h"

#define FACTOR_MS 1000
#define FACTOR_RATIO_RR 1000
#define PERCENTILE_LOO_LOW 643
#define PERCENTILE_LOO_HIGH 1473

int errorDetection(int32_t *arg[]);

#endif // ERROR_DETECTION_H_
