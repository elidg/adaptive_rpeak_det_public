#ifndef MORPHO_FILTERING_H_
#define MORPHO_FILTERING_H_
#include "defines_globals.h"

#if SAMPLING_FREQUENCY==250

#define BUFFER_BITS_FILTER 0   /* Number of bits in the size of the buffer */
#define BUFFER_SIZE_FILTER (1 << BUFFER_BITS_FILTER)


//R: Size of the raw input buffer
#define TAMBEnt 150

/**************************
 * R: Baseline Removal
 */
//R: Length of the opening structuring element
#define TAMB0 50
// Half length of the opening structuring element
#define TAMB0m 25
//R: Length of the closing structuring element
#define TAMBc 75
//R: Half length of the closing structuring element
#define TAMBcm 37


#if HIGH_FREQ==1

//R: Size of the Buffer after baseline removal
#define TAMBFB 5
//R: Half size of the Buffer after baseline removal
#define TAMBFBm 0

/*************************
 * R: Noise suppression
 */
//R: Length of the B1 structuring element
#define TAMB1 5

#define B1_STRUCTURING_ELEMENT {0,1,5,1,0,\
							   0,0,1,5,1,\
							   1,0,0,1,5,\
							   5,1,0,0,1,\
							   1,5,1,0,0};

#endif  // HIGH_FREQ

#elif SAMPLING_FREQUENCY==500

#define BUFFER_BITS_FILTER 0   /* Number of bits in the size of the buffer */
#define BUFFER_SIZE_FILTER (1 << BUFFER_BITS_FILTER)

//R: Size of the raw input buffer
#define TAMBEnt 300

/**************************
 * R: Baseline Removal
 */
//R: Length of the opening structuring element
#define TAMB0 100
// Half length of the opening structuring element
#define TAMB0m 50
//R: Length of the closing structuring element
#define TAMBc 150
//R: Half length of the closing structuring element
#define TAMBcm 75

#if HIGH_FREQ==1

//R: Size of the Buffer after baseline removal
#define TAMBFB 9
//R: Half size of the Buffer after baseline removal
#define TAMBFBm 4

/*************************
 * R: Noise suppression
 */

//R: Length of the B1 structuring element
#define TAMB1 9

#define B1_STRUCTURING_ELEMENT {0,0,1,2,5,2,1,0,0,\
							    0,0,0,1,2,5,2,1,0,\
							    0,0,0,0,1,2,5,2,1,\
							    1,0,0,0,0,1,2,5,2,\
							    2,1,0,0,0,0,1,2,5,\
							    5,2,1,0,0,0,0,1,2,\
							    2,5,2,1,0,0,0,0,1,\
							    1,2,5,2,1,0,0,0,0,\
							    0,1,2,5,2,1,0,0,0}


#endif  // HIGH_FREQ
#endif  // SAMPLING_FREQUENCY

int16_t filterSample(int16_t offset, int16_t read_data, int16_t index);
void init_filtering();
void filterWindows(int32_t *arg[]);

#endif // MORPHO_FILTERING_H_
