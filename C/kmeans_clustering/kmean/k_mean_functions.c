#include "k_mean_functions.h" 
#include <math.h>
#include "../search/search_functions.h"

type_f genlogfunc(type_f x, type_f k, type_f b, type_f nu, type_f m){
	return (k/powf(1+expf(-b*(x-m)),(1/nu)));
}

 
type_f gaussian(type_f x, type_f mu, type_f sign){
	return expf(-powf((x - mu), 2)/(2*powf(sign, 2)));
}

/* Function to sort an array using insertion sort*/
void insertionSort(type_f arr[], int n) 
{ 
	type_f key;
    int i, j; 
    for (i = 1; i < n; i++) { 
        key = arr[i]; 
        j = i - 1; 
  
        /* Move elements of arr[0..i-1], that are 
          greater than key, to one position ahead 
          of their current position */
        while (j >= 0 && arr[j] > key) { 
            arr[j + 1] = arr[j]; 
            j = j - 1; 
        } 
        arr[j + 1] = key; 
    } 
} 

type_f abs_s2_init[SIZE_PERCENTILE_ARR]; //it's the first 1000 samples of abs diff to sort and find the 99 percentile

type_f percentile(int16_t *x, int end, int percToFind){

	int indexPerc = (end-1)*percToFind/100;
	for(int i = 0; i<end-1; i++){
		abs_s2_init[i] = fabs((type_f) (x[i+1]-x[i]));
	}
	
	insertionSort(abs_s2_init, end-1);

	return abs_s2_init[indexPerc];

}


void partial_fit(type_f *x, type_f *hcentr, type_f *lcentr, int *h_el, int *l_el, int *label){

	//Euclidian Distance and update
	type_f dist[2] = {0};
	type_f tmp1, tmp2;

 	tmp1 = (*hcentr-*x); 
 	tmp2 = (*lcentr-*x);
	dist[0] = sqrtf(tmp1*tmp1);
	dist[1] = sqrtf(tmp2*tmp2);

	if(dist[0] < dist[1]){
		 
		*hcentr = ((*hcentr*(*h_el)) + *x)/(*h_el+1);
		*h_el=*h_el+1;
		*label = 1;
	}else{
		*lcentr = ((*lcentr*(*l_el)) + *x)/(*l_el+1);
		*l_el=*l_el+1;
		*label = 0;
	}

}


float std(type_i *x, int num_peaks, int num_peaks_total){

	type_f mean = 0.0;
	type_f std = 0.0;
	int stop_std = 0;
	type_f rr[PEAKS_STD];

	if(num_peaks_total<5){ //Considering 5 peaks it's 4 RR intervals
		stop_std = num_peaks_total-1;
	}else{
		stop_std = 4;
	}

	for(int i = 0; i < stop_std; i++){
		rr[i] = (type_f) (x[num_peaks-i]-x[num_peaks-i-1]);
		mean += rr[i]; 
#ifdef PRINT_DEBUG_CL_STD	
		// printf("rr[i]: %f\n",rr[i] );
#endif
	}
	mean = mean/stop_std;
	 
	for(int i = 0; i<stop_std; i++){
		std = std + (rr[i] - mean)*(rr[i] - mean); 
	}
	 
	if(stop_std==1)
		return 0;
	else
		return sqrtf(std/stop_std); // because the python version divides by N and not N-1

}
  
int flag_initial_centr = 0;
int label = 0;
int zeroctr = 0;
type_f hcentr = 0.0f;
type_f lcentr = 1.0f; 
int h_el=1; 
int l_el=1;
type_f mu = 200.0f;
type_f sd = 25.0f;
type_i last_peak = 0;
type_i new_peak = 0;
type_i end_last_peak = 0;
type_f centroids_state[4];
type_f bf = 0.0f;
type_f bt = 0.0f;
type_i in_qrs = 0;
type_f st[52]; // 52 elements because in_qrs starts after 51 samples max and then the array will always be reset. Once a QRS has been detected, we look for 120ms of 0 output (30 samples) or 200ms to consider that the QRS has finished.
type_f s2[52];
int index_st = 0;
int qrs_init = 0;
int qrs_init_total = 0;
int num_peaks = 0;
int num_peaks_total = 0;
int last_peak_w = 0;
type_i rpks[MAX_PEAKS_IN_WIND];
int end_qrs = 0;
int max_min_slope[2];
int max_min_slope_productSign[2];

void rpeaks(int32_t *arg[]){
	
	int16_t* ecg_buff_full = (int16_t*) arg[0];
	int number_of_window = *(arg[1]); 
	int start_index_buff = *(arg[2]);
	int end_main_loop = *(arg[3]) - start_index_buff;
	int start_index_w = *(arg[4]);
	int flag_prev_error = *(arg[5]);
	int32_t* indicesRpeaks = arg[6];
	int32_t* overlap = arg[7];
	int32_t* overlap_qrs_init = arg[8];
	int ind_peak = 0;

    for(int ix_rr = 0; ix_rr<H_B+1; ix_rr++)
        indicesRpeaks[ix_rr] = 0;

	if(flag_prev_error==1){
		start_index_buff = start_index_buff-*overlap;
		start_index_w = start_index_w-*overlap;
	}else{
		start_index_w = start_index_w - DIM;
		*overlap = 0;
	}

	int16_t *ecg_buff = &ecg_buff_full[start_index_buff];

	type_f x = 0.0f;
	type_f param_log = 0.0f;
	type_f param_log2 = 0.0f;
	type_f new_val = 0.0f;

#ifdef PRINT_DEBUG_CL

	type_i num_wind_to_check = 100;//63;//71; // Subject 3 - 20: Window where peak is below the baseline and there is a problem in detecting it (SOLVED)
								  // Subject 3 - 2: window with another problem at the beginning (check gaussian) (SOLVED)
								  // Subject 3 - 136: SOLVED
#endif	

	if(flag_initial_centr == 0){
		hcentr = percentile(ecg_buff, SIZE_PERCENTILE_ARR+1, 99);
    	flag_initial_centr = 1;
#ifdef PRINT_INITIAL_HCENTR
    	printf("initial hcentr: %f\n",hcentr);
#endif    	
	}

#ifdef PRINT_DEBUG_CL
	printf("number_of_window: %d start_index_buff: %d end_main_loop: %d start_index_w: %d overlap: %d\n", number_of_window,start_index_buff,end_main_loop + *overlap,start_index_w,*overlap);
#endif	

#ifdef PRINT_INPUT_SIG_CL
	if(number_of_window == num_wind_to_check){
		for(int i = 0; i < end_main_loop + *overlap; i++){
			printf("%d,",ecg_buff[i] );
		}
		printf("\n");
	}
#endif

	for(int i = 1; i < end_main_loop + *overlap; i++){

		x = fabs((type_f) (ecg_buff[i] - ecg_buff[i-1]));

#ifdef PRINT_ABS_DIFF		
		// DEBUG ABS DIFF
		printf("%f,", x);
#endif

		//====== Uncomment if you are using end_last_peak for overlap ======//
		if(*overlap==end_main_loop - end_last_peak - 1 && i==1){
			hcentr = centroids_state[0];
			lcentr = centroids_state[1];
			h_el = (type_i) centroids_state[2];
			l_el = (type_i) centroids_state[3];
 		}
 		//====== Uncomment if you are using end_last_peak for overlap ======//

		//====== Uncomment if you are using qrs_init for overlap ======//
// 		if(number_of_window > 0 && i==1 && last_peak > start_index_w){
// 			last_peak = rpks[num_peaks-2];
// 			mu = rpks[num_peaks-2] - rpks[num_peaks-3];
// 			sd = std(rpks, num_peaks-2, num_peaks_total-2);
// #ifdef PRINT_DEBUG_CL_STD		
// 			printf("last_peak: %d mu: %f sd: %f\n",last_peak, mu, sd);
// #endif			
//        		if(sd < 10)
//                 sd = 10;
//             if(sd > 50)
//                 sd = 50;

//             num_peaks = num_peaks-1;
//             last_peak_w = last_peak;
// 		}
		//====== Uncomment if you are using qrs_init for overlap ======//
#ifdef PRINT_DEBUG_CL
		if(number_of_window==num_wind_to_check)  	
			printf("last_peak: %d input_value_gauss: %d mu: %f sd: %f\n",last_peak, start_index_w + (i-1) - last_peak, mu, sd);
#endif		
		bf = gaussian(start_index_w + (i-1) - last_peak, mu, sd); // i - last_peak_w
#ifdef PRINT_GAUSS
		// DEBUG GAUSSIAN
		printf("%f,", bf);
		#endif		       
#ifdef PRINT_GAUSS_MU_SD
		printf("%f %f %f \n", bf,mu,sd);
#endif
		param_log = logf((hcentr - lcentr)/2)/(2*lcentr);
			              
		param_log2 = log2f(hcentr/lcentr);

	 	bt = genlogfunc(x , hcentr, param_log, 1/param_log2, lcentr);

#ifdef PRINT_GENLOGFUN
	 	printf("%f,", bt);
#endif	 	
	            
		new_val = (bt * bf);

		if(x > new_val){
			partial_fit(&x, &hcentr, &lcentr, &h_el, &l_el, &label);
#ifdef PRINT_DEBUG_CL		
			if(number_of_window==num_wind_to_check)	
				printf("index_all: %d i: %d x: %f hcentr: %f lcentr: %f label: %d\n", start_index_w + (i-1), i, x, hcentr, lcentr, label);
#endif		
#ifdef PRINT_BAYFILT
			// DEBUG BAYESIAN FILTERED SIGNAL 
			printf("%f,", x);
#endif			
		}
		else{
			partial_fit(&new_val, &hcentr, &lcentr, &h_el, &l_el, &label);
#ifdef PRINT_DEBUG_CL			
			if(number_of_window==num_wind_to_check)
				printf("index_all: %d i: %d new_val: %f hcentr: %f lcentr: %f label: %d\n", start_index_w + (i-1), i, new_val, hcentr, lcentr, label);
#endif			
#ifdef PRINT_BAYFILT	
			// DEBUG BAYESIAN FILTERED SIGNAL 
			printf("%f,", new_val);
#endif			
		}
#ifdef PRINT_CENTROIDS	
	    // DEBUG CENTROIDS
	    printf("%f %f\n", lcentr,hcentr);
#endif
	          
		if(in_qrs == 1){ 
			if(label == 0){
				zeroctr += 1; 
			}else{
	     		zeroctr = 0;
	    	}

			s2[index_st] = (type_f) (ecg_buff[i] - ecg_buff[i-1]);
    		if(x > new_val){
            	st[index_st] = x;
				index_st++;
            }else{
            	st[index_st] = new_val;
				index_st++;
            }

#ifdef PRINT_DEBUG_CL	
			if(number_of_window==num_wind_to_check)    	
	        	printf("qrs_init: %d zeroctr: %d qrs_init_total: %d\n", qrs_init, zeroctr, qrs_init_total);
#endif	        
	        // Once a QRS has been detected, we look for 120ms of 0 output (30 samples) or 200ms to consider that the QRS has finished.
	   		if(zeroctr==30 || start_index_w + (i-1) - qrs_init_total > MAX_QRS_DUR){
	                  
#ifdef PRINT_DEBUG_CL
	   			if(number_of_window==num_wind_to_check) 
	   				printf("start interval: %d stop interval: %d ecg_buff[start_interval]: %f i: %d zeroctr: %d label: %d  index_st: %d  stop_interval_st: %d \n",qrs_init, (i-1) - zeroctr + 1,ecg_buff[qrs_init],i, zeroctr, label,index_st-1,(index_st-1)-zeroctr+1 );
#endif	   			

				// Search for location of the new peak. Checking maximum and minimum slopes, and move towards the peak
	   			end_qrs = (index_st-1)-zeroctr+1;
	   			argMaxMinProductSignSearch(st, s2, 0, end_qrs,max_min_slope_productSign);
	   			int max_slope = max_min_slope_productSign[0];
	   			int min_slope = max_min_slope_productSign[1];
	   			int max_slope_s2 = 0;
	   			int min_slope_s2 = 0;
	   			if (max_slope < min_slope){
                    new_peak = qrs_init + max_slope + argMaxProductSignSearch(st, s2, max_slope, min_slope+1);
                }else if (max_slope > min_slope){
                	argMaxMinSearch(s2, 0, end_qrs, max_min_slope);
                	max_slope_s2 = max_min_slope[0];
                    min_slope_s2 = max_min_slope[1];
                    if(max_slope_s2 < min_slope_s2){
                       	// Q wave
                        new_peak = qrs_init + max_slope_s2 + argMaxSearch(s2, max_slope_s2, min_slope_s2+1);
                    }else{
                        // S wave
                        new_peak = qrs_init + argMaxSearch(s2, 0, min_slope_s2+1);
                    }
                }else{
                    new_peak = qrs_init + argMaxSearch(st, 0, end_qrs);
                }

#ifdef PRINT_DEBUG_CL
	        	if(number_of_window==num_wind_to_check){  
	        		printf("BEFORE CORRECTION new_peak: %d\n", new_peak);
	        		if (max_slope > min_slope){
	        			printf("end_qrs: %d max_slope: %d min_slope: %d max_slope_s2: %d min_slope_s2: %d\n",end_qrs,max_slope,min_slope,max_slope_s2,min_slope_s2);
	        		}else{
	        			printf("end_qrs: %d max_slope: %d min_slope: %d\n",end_qrs,max_slope,min_slope);
	        		}
	        	}
#endif	        	

	        	if(ecg_buff[new_peak+1] > ecg_buff[new_peak]){
                    while(ecg_buff[new_peak+1] > ecg_buff[new_peak]){
                        new_peak += 1;
                    }
				}
                else{
                    while(ecg_buff[new_peak] < ecg_buff[new_peak-1]){
                        new_peak -= 1;
                    }
                }
	                  	                 	           

#ifdef PRINT_DEBUG_CL
	        	if(number_of_window==num_wind_to_check)  
	        		printf("AFTER CORRECTION new_peak: %d\n", new_peak);
#endif	   	        	

	         	last_peak = new_peak + start_index_w; 	      
	         	if(last_peak > last_peak_w){
					if(last_peak!=rpks[num_peaks]){	         				        		
		        		indicesRpeaks[ind_peak] = last_peak;
		        		ind_peak++;
#ifdef PRINT_RPEAKS_DEBUG	       
		        		printf("ecg_buff_peak: %d last_peak: %d\n", new_peak, last_peak); 
#endif			        		
		        	}
	         		if(num_peaks == MAX_PEAKS_IN_WIND){
						for(int ix=0; ix<PEAKS_STD; ix++){
							rpks[ix] = rpks[MAX_PEAKS_IN_WIND - PEAKS_STD + ix];
						}	
						num_peaks = 5;
					}
		       		rpks[num_peaks] = last_peak;
#ifdef PRINT_DEBUG_CL_STD		
					if(number_of_window==num_wind_to_check)         		
		       			printf("num_peaks: %d rpks[num_peaks]: %d\n", num_peaks,rpks[num_peaks]);
#endif		       		
		       		num_peaks_total++;

		        	if(num_peaks_total > 1){
		          		mu = rpks[num_peaks]-rpks[num_peaks-1]; 
	                	sd = std(rpks, num_peaks, num_peaks_total); 
	#ifdef PRINT_DEBUG_CL_STD
	                	if(number_of_window==num_wind_to_check)  
	                		printf("mu: %f sd: %f\n", mu, sd);
	#endif
	               		if(sd < 10)
	                        sd = 10;
	                    if(sd > 50)
	                        sd = 50;
		            }
		            num_peaks++;
				}
				last_peak_w = last_peak;
				end_last_peak = i-1;
				centroids_state[0] = hcentr;
				centroids_state[1] = lcentr;
				centroids_state[2] = (type_f) h_el; 
				centroids_state[3] = (type_f) l_el;
				zeroctr = 0;
	 			in_qrs = 0;
	 			index_st = 0;
	        }
	               

	    }else{
	        if(label==1 && start_index_w + (i-1) > last_peak + MIN_RR_DIST){
	            in_qrs = 1;
	            qrs_init = i-1; 
	            qrs_init_total = start_index_w + (i-1);
	            s2[index_st] = (type_f) (ecg_buff[i] - ecg_buff[i-1]);
	            if(x > new_val){
	            	st[index_st] = x;
					index_st++;
	            }else{
	            	st[index_st] = new_val;
					index_st++;
	            }
	        }
	    }
	}

#ifdef PRINT_DEBUG_CL
	if(number_of_window==num_wind_to_check-1) 
		printf("END OF WINDOW qrs_init: %d qrs_init_total: %d last_peak: %d mu: %f std: %f index_st: %d new_peak: %d\n",qrs_init,qrs_init_total,last_peak,mu,sd, index_st-1,new_peak);
#endif

	if(in_qrs == 1)
		*overlap_qrs_init = end_main_loop + *overlap - qrs_init - 1;
	else
		*overlap_qrs_init = 0;

#ifdef PRINT_DEBUG_CL
	if(number_of_window==num_wind_to_check-1) 
		printf("end_main_loop: %d overlap: %d qrs_init: %d overlap_qrs_init: %d\n",end_main_loop,*overlap,qrs_init, *overlap_qrs_init);
#endif

	if(last_peak < start_index_w){
		if(in_qrs == 0)
			*overlap = 1;
		else
			*overlap = end_main_loop + *overlap - (qrs_init - MAX_QRS_DUR) - 1;
	}else{
		*overlap = end_main_loop + *overlap - end_last_peak - 1; // - qrs_init; // Using qrs_init gives some problems if the peak is a little bit before than qrs_init. This logically shouln't happen but unfortunately it happens because the python code is not perfect (probably needing some tweeking in the variable representing the duration of the qrs)
	}

	in_qrs = 0;
	zeroctr = 0; //if overlap = end_main_loop - end_last_peak - 1; 
	qrs_init = 0;	
	index_st = 0;	

}
 
