#include "test_double_buffering.h"
#include "rt/rt_api.h"
#include <stdio.h>

void testDoubleBuff(int32_t* arg[]){
	int16_t* ecg_buff = (int16_t*) arg[0];
	int32_t start_index_w = *(arg[1]);
	int32_t end_main_loop = *(arg[2]);

	// printf("start_index_w: %d end_main_loop: %d\n",start_index_w,end_main_loop );
	// for(int i=start_index_w; i<end_main_loop;i++)
	// 	printf("%d\n",ecg_buff[i]);
}