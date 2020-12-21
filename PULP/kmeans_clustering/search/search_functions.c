#include <math.h>
#include "search_functions.h"

type_i argMaxSearch(type_f* x, type_i ind_start, type_i ind_end){
	type_i argmax = 0;
	type_f max_value = x[ind_start];
	for(type_i ix = ind_start+1; ix < ind_end; ix++){
		if(x[ix]>max_value){
			max_value = x[ix];
			argmax = ix;
		}
	}
	return argmax;
}

void argMaxMinSearch(type_f* x, type_i ind_start, type_i ind_end, type_i* argmax_argmin){
	argmax_argmin[0] = ind_start;
	argmax_argmin[1] = ind_start;
	type_f max_value = x[ind_start];
	type_f min_value = x[ind_start];
	for(type_i ix = ind_start+1; ix < ind_end; ix++){
		if(x[ix]>max_value){
			max_value = x[ix];
			argmax_argmin[0] = ix;
		}
		if(x[ix]<min_value){
			min_value = x[ix];
			argmax_argmin[1] = ix;
		}
	}	
}

type_i sign(type_f x){
	if(x>0) return 1;
	if(x<0)	return -1;

	return 0;
}

type_i argMaxProductSignSearch(type_f* x, type_f* x2, type_i ind_start, type_i ind_end){
	type_i argmax = 0;
	type_f max_value = 0;
	type_f product_arrays = 0;

	max_value = x[ind_start]*sign(x2[ind_start]);

	for(type_i ix = ind_start+1; ix < ind_end; ix++){
		product_arrays = x[ix]*sign(x2[ix]);

		if(product_arrays>max_value){
			max_value = product_arrays;
			argmax = ix;
		}
	}
	return argmax;
}

void argMaxMinProductSignSearch(type_f* x, type_f* x2, type_i ind_start, type_i ind_end, type_i* argmax_argmin){
	argmax_argmin[0] = ind_start;
	argmax_argmin[1] = ind_start;
	type_f max_value = 0;
	type_f min_value = 0;
	type_f product_arrays = 0;

	max_value = x[ind_start]*sign(x2[ind_start]);
	min_value = x[ind_start]*sign(x2[ind_start]);

	for(type_i ix = ind_start+1; ix < ind_end; ix++){
		product_arrays = x[ix]*sign(x2[ix]);

		if(product_arrays>max_value){
			max_value = product_arrays;
			argmax_argmin[0] = ix;
		}

		if(product_arrays<min_value){
			min_value = product_arrays;
			argmax_argmin[1] = ix;
		}
	}	
}