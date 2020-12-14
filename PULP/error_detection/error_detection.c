#include "./error_detection.h"

int errorDetection(int32_t *arg[]){
	int32_t *r_counter = arg[0];
    int32_t *rWindow = arg[1];
    int32_t *lastRRp = arg[2];
    int32_t *indRpeaks = arg[3];
    int32_t *offset_ind_r = arg[4];    
    int32_t *lastPeak = arg[5];
    int32_t RR_intervals[H_B+2];
    int32_t ratioConsecutiveRR = 0;
    int32_t offset_ind_rr = -1;
    int32_t rr_counter = 0;

    for(int32_t ix_rr = 0; ix_rr < H_B ; ix_rr++) {
       RR_intervals[ix_rr] = 0;
    }

    if(*r_counter>0){
#ifdef PRINT_DEBUG_ERRDET
        printf("rWindow: %d r_counter: %d\n", *rWindow,*r_counter );
#endif        
        if(!(*rWindow == 0)){
            if(*lastRRp == 0){
                RR_intervals[0] = (FACTOR_MS*(indRpeaks[0] + *offset_ind_r - *lastPeak))/ECG_SAMPLING_FREQUENCY;
                offset_ind_rr = 0;
#ifdef PRINT_DEBUG_ERRDET
                printf("R[0]: %d RR[0]: %d\n", (FACTOR_MS*(indRpeaks[0] + *offset_ind_r))/ECG_SAMPLING_FREQUENCY, RR_intervals[0]);
#endif                
            }
            else{
                RR_intervals[0] = *lastRRp;
                RR_intervals[1] = (FACTOR_MS*(indRpeaks[0] + *offset_ind_r - *lastPeak))/ECG_SAMPLING_FREQUENCY;
                offset_ind_rr = 1;
#ifdef PRINT_DEBUG_ERRDET
                printf("R[0]: %d RR[0]: %d RR[1]: %d\n", (FACTOR_MS*(indRpeaks[0] + *offset_ind_r))/ECG_SAMPLING_FREQUENCY, RR_intervals[0], RR_intervals[1]);
#endif               
            }
        }

        for(int32_t ix_rp = 1; ix_rp < *r_counter; ix_rp++) {
            RR_intervals[ix_rp+offset_ind_rr] = (FACTOR_MS*(indRpeaks[ix_rp]-indRpeaks[ix_rp-1]))/ECG_SAMPLING_FREQUENCY;
#ifdef PRINT_DEBUG_ERRDET
            printf("R[%d]: %d RR[%d]: %d\n", ix_rp,(FACTOR_MS*(indRpeaks[ix_rp] + *offset_ind_r))/ECG_SAMPLING_FREQUENCY,ix_rp+offset_ind_rr,RR_intervals[ix_rp+offset_ind_rr]);           
#endif             
        }

        *lastPeak = indRpeaks[*r_counter-1] + *offset_ind_r;
        *lastRRp = RR_intervals[*r_counter-1];
#ifdef PRINT_DEBUG_ERRDET
        printf("lastPeak: %d lastRR: %d\n",(FACTOR_MS*(*lastPeak))/ECG_SAMPLING_FREQUENCY,*lastRRp );
#endif        

        if(*rWindow==0)
            rr_counter = *r_counter-1;
        else
            rr_counter = *r_counter;
        
        for(int32_t ix_rr = 1; ix_rr < rr_counter; ix_rr++){
            ratioConsecutiveRR = (FACTOR_RATIO_RR*RR_intervals[ix_rr])/RR_intervals[ix_rr-1];
#ifdef PRINT_DEBUG_ERRDET
            printf("ratioConsecutiveRR: %d\n", ratioConsecutiveRR);
#endif            
            if(ratioConsecutiveRR<PERCENTILE_LOO_LOW || ratioConsecutiveRR>PERCENTILE_LOO_HIGH)
                return 1;
        }
    }else{
        return 1;
    }

    return 0;
}