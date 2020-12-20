#include <stdint.h>
#include "rt/rt_api.h"
#include "rt/rt_omp.h"

#define BI
void __attribute__ ((noinline)) shift16();
//#define GPIO
/*GPIO definitions*/
#define PAD_GPIO_NUMBER		9
#define GPIO_NUMBER			9

/*Other definitions*/
#define RESET_FLAG 			0
#define SET_FLAG			1


//FOR PARALLEL EXECUTION
#define STACK_SIZE      2048
#define MOUNT           1
#define UNMOUNT         0
#define CID             0


#if TARGET == 1
#define CL
#else
#define FC 
#endif
