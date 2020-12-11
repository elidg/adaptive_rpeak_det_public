#include "adaptiveRpeakDetection.h"
#include "../defines.h"
#include "../profiling/profile.h"
#include "../data/signal.h"
#include "../Morph_filt/morpho_filtering.h"
#include "../Morph_filt/defines_globals.h"

#define N_WINDOWS (int) (2*((ECG_VECTOR_SIZE-LONG_WINDOW)/dim)+1)// Counting the worst case scenario when overlap is dim
 
RT_L2_DATA int16_t ecg_buff[(LONG_WINDOW+dim)*(NLEADS+1)];

RT_L2_DATA rt_perf_t perf[NUM_CORES];

#ifdef MODULE_MF
RT_L2_DATA int32_t *argMF[4];
RT_L2_DATA int32_t buffSize_windowMF;
#endif

#ifdef MODULE_RELEN
RT_L2_DATA int32_t *argRelEn[4];
RT_L2_DATA int32_t start_RelEn = 1; 
RT_L2_DATA int32_t buffSize_windowRelEn;
#endif

#ifdef MODULE_RPEAK_REWARD
RT_L2_DATA int32_t *argRW_Rpeak[2];
RT_L2_DATA int32_t indicesRpeaks[H_B+1];
#endif

#ifdef MODULE_ERROR_DETECTION
RT_L2_DATA int32_t *argErrDet[4];
RT_L2_DATA int32_t lastRpeak = 0;
RT_L2_DATA int32_t lastRR = 0;
RT_L2_DATA int32_t error_RWindow = 0;

#endif

RT_L2_DATA int32_t overlap;
RT_L2_DATA int32_t rWindow;


void clearRelEn() {
    clearAndResetRelEn();
    resetPeakDetection();
}

int errorDetection(int32_t *arg[]){
    int32_t *indRpeaks = arg[0];
    int32_t *r_counter = arg[1];
    int32_t *lastPeak = arg[2];
    int32_t *lastRRp = arg[3];
    int32_t RR_intervals[H_B+2];
    int32_t ratioConsecutiveRR = 0;
    int32_t offset_ind_rr = 1;

    for(int32_t ix_rr = 0; ix_rr < H_B ; ix_rr++) {
       RR_intervals[ix_rr] = 0;
    }

    if(*r_counter>0){

        if(!(rWindow == 1 && *r_counter == 1)){
            if(*lastRRp = 0){
                RR_intervals[0] = (FACTOR_MS*(indRpeaks[0] - *lastPeak))/ECG_SAMPLING_FREQUENCY;
                offset_ind_rr = 0;
            }
            else{
                RR_intervals[0] = *lastRRp;
                RR_intervals[1] = (FACTOR_MS*(indRpeaks[0] - *lastPeak))/ECG_SAMPLING_FREQUENCY;
                offset_ind_rr = 1;
            }
            printf("RR[0]: %d\n",RR_intervals[0] );
            printf("RR[1]: %d\n",RR_intervals[1] );
        }

        for(int32_t ix_rp = 1; ix_rp < *r_counter; ix_rp++) {
            RR_intervals[ix_rp+offset_ind_rr] = (FACTOR_MS*(indRpeaks[ix_rp]-indRpeaks[ix_rp-1]))/ECG_SAMPLING_FREQUENCY;
        }

        for(int32_t ix_rr = 1; ix_rr < *r_counter; ix_rr++){
            ratioConsecutiveRR = (FACTOR_RATIO_RR*RR_intervals[ix_rr])/RR_intervals[ix_rr-1];
            if(ratioConsecutiveRR<PERCENTILE_LOO_LOW || ratioConsecutiveRR>PERCENTILE_LOO_HIGH)
                return 1;
        }

        *lastPeak = indRpeaks[*r_counter-1];
        *lastRRp = RR_intervals[*r_counter-1];
    }else{
        return 1;
    }

    return 0;
}

void adaptiveRpeakDetection(){
   
    int32_t count_sample = 0;
    int32_t offset_window = 0;
    int32_t offset_ind = LONG_WINDOW/2+1;
    int32_t tot_overlap = 0;

    overlap = 0;

#ifdef MODULE_MF
    int32_t flagMF = 0;
    int32_t i_lead = 0;

    buffSize_windowMF = LONG_WINDOW+dim;

    argMF[0] = (int32_t*) ecg_buff;
    argMF[1] = &flagMF;
    argMF[2] = &i_lead;
    argMF[3] = &buffSize_windowMF;    

    init_filtering();
#endif

#ifdef MODULE_RELEN    
    buffSize_windowRelEn = LONG_WINDOW+dim;

    argRelEn[0] = (int32_t*) ecg_buff;
    argRelEn[1] = (int32_t*) &ecg_buff[(LONG_WINDOW+dim)*NLEADS];
    argRelEn[2] = &start_RelEn;
    argRelEn[3] = &buffSize_windowRelEn;

    clearRelEn();
#endif

#ifdef MODULE_RPEAK_REWARD 
    int32_t rpeaks_counter = 0;

    argRW_Rpeak[0] = (int32_t*) &ecg_buff[LONG_WINDOW+(LONG_WINDOW + dim)*NLEADS];
    argRW_Rpeak[1] = indicesRpeaks;
#endif


    for(rWindow=0; rWindow<N_WINDOWS; rWindow++)
    {
        if((rWindow+1)*dim + LONG_WINDOW -1 - tot_overlap >= ECG_VECTOR_SIZE){
            return;
        }

        if(rWindow > 0){
            offset_window = LONG_WINDOW;
        }
        else{
            offset_window = 0;
        }        

        if (rWindow > 0) {
            for(int32_t i=0; i<overlap; i++) {
#ifdef OVERLAP_MF
                ecg_buff[i+offset_window] = ecg_buff[(LONG_WINDOW + dim - overlap + i + offset_window)];
#endif
#ifdef OVERLAP_RELEN
                ecg_buff[i + offset_window + (LONG_WINDOW + dim)*NLEADS] = ecg_buff[(2*(LONG_WINDOW + dim)*NLEADS - overlap + i + offset_window)];
#endif                
            }
        }

        for(int32_t lead=0; lead<NLEADS; lead++) {
            for(int32_t i=overlap + offset_window; i< LONG_WINDOW + dim; i++) {
                ecg_buff[i + (dim+LONG_WINDOW)*lead] = ecg_1l[rWindow*dim + i - tot_overlap]; 
            }
        }

#ifdef PRINT_DEBUG
        printf("start_window: %d end_window: %d overlap: %d\n", offset_window + rWindow*dim - tot_overlap,LONG_WINDOW -1 + (rWindow+1)*dim - tot_overlap,overlap);
#endif

#ifdef MODULE_MF

        if (rWindow == 0) {
            // Needed to initialize the MF filter properly
            for(int32_t lead=0; lead<NLEADS; lead++) {
                for(int32_t i=0; i<=OFFSET_MF; i++) {
                    ecg_buff[i + (dim+LONG_WINDOW)*lead] = 0;
                }
            }
        }

        // MF on FC because not enough memory on Cluster
        argMF[0] = (int32_t*) &ecg_buff[offset_window + overlap];
        buffSize_windowMF = LONG_WINDOW - offset_window + dim-overlap;

        for(i_lead = 0; i_lead < NLEADS; i_lead++){
            filterWindows(argMF);
        }

        flagMF=1;

    #ifdef PRINT_SIG_MF
        for(int32_t sample = offset_window; sample<LONG_WINDOW + dim; sample++) {
            printf("%d\n", ecg_buff[sample]);
        }
    #endif

#endif

#ifdef HWPERF_FULL
        profile_start(perf);
#endif

#ifdef MODULE_RELEN

    #ifdef HWPERF_MODULES
        profile_start(perf);
    #endif

#ifdef OVERLAP_MF        
        clearAndResetRelEn();
        start_RelEn = 1;
#endif       

#ifdef OVERLAP_RELEN
        if(rWindow > 0)
            start_RelEn = 0;

        argRelEn[0] = (int32_t*) &ecg_buff[offset_window + overlap];
        argRelEn[1] = (int32_t*) &ecg_buff[offset_window + (LONG_WINDOW + dim)*NLEADS+overlap];
        buffSize_windowRelEn = LONG_WINDOW - offset_window + dim-overlap;
#endif        

        relEn_w(argRelEn);       

    #ifdef HWPERF_MODULES
        profile_stop(perf);
    #endif

    #ifdef PRINT_RELEN
        for(int32_t sample = offset_window + (LONG_WINDOW + dim)*NLEADS; sample< (LONG_WINDOW + dim)*(NLEADS+1); sample++) {
            printf("%d\n", ecg_buff[sample]);
        }
    #endif

#endif

#ifdef MODULE_RPEAK_REWARD

    #ifdef HWPERF_MODULE_RPEAK_REWARD || HWPERF_MODULES
        profile_start(perf);
    #endif

        getPeaks_w(argRW_Rpeak);

        rpeaks_counter = 0;

        while(indicesRpeaks[rpeaks_counter]!=0) {
            rpeaks_counter++;
        }

    #ifdef HWPERF_MODULE_RPEAK_REWARD || HWPERF_MODULES
        profile_stop(perf);
    #endif

    #ifdef PRINT_RPEAKS
        for(int32_t indR=0; indR<rpeaks_counter; indR++) {
            printf("%d\n", (indicesRpeaks[indR]) + offset_ind);
        }
    #endif

#endif

#ifdef MODULE_ERROR_DETECTION

        argErrDet[0] = indicesRpeaks;
        argErrDet[1] = &rpeaks_counter;
        argErrDet[2] = &lastRpeak;
        argErrDet[3] = &lastRR;
        error_RWindow = errorDetection(argErrDet);

    #ifdef PRINT_ERROR_RPEAKS
        printf("%d\n", error_RWindow);
    #endif
#endif        

#ifdef ONLY_FIRST_WINDOW //Only for debug
    return;
#endif

#ifdef OVERLAP_MF    
        overlap = LONG_WINDOW + LONG_WINDOW/2 + 1;            
#endif
#ifdef OVERLAP_RELEN
        overlap = 0;
#endif

        tot_overlap += overlap;
        offset_ind = offset_ind + dim - tot_overlap;        

#ifdef MODULE_RPEAK_REWARD        
        rpeaks_counter = 0;

        for(int32_t ix_rp = 0; ix_rp < H_B+1 ; ix_rp++) {
           indicesRpeaks[ix_rp] = 0;
        }
#endif        
    }

}
