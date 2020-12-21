#include "rt/rt_api.h"
#include "adaptiveRpeakDetection.h"
#include "../defines.h"
#include "../profiling/profile.h"
#include "../profiling/profile_cl.h"
#include "../profiling/defines.h"
#include "../data/signal.h"
#include "../Morph_filt/morpho_filtering.h"
#include "../Morph_filt/defines_globals.h"
#include "../error_detection/error_detection.h"
#include "../kmeans_clustering/kmean/k_mean_functions.h"

#define N_WINDOWS (int) (2*((ECG_VECTOR_SIZE-LONG_WINDOW)/DIM)+1)// Counting the worst case scenario when overlap is DIM
 
RT_L2_DATA int16_t ecg_L2buff[(LONG_WINDOW+DIM)*(NLEADS+1)];

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
RT_L2_DATA int32_t *argRW_Rpeak[3];
RT_L2_DATA int32_t indicesRpeaks[H_B+1];
#endif

#ifdef MODULE_ERROR_DETECTION
RT_L2_DATA int32_t *argErrDet[5];
RT_L2_DATA int32_t lastRpeak = 0;
RT_L2_DATA int32_t lastRR = 0;
RT_L2_DATA int32_t error_RWindow = 0;
#endif

#ifdef MODULE_CLUSTERING
RT_L2_DATA int32_t rL2BufferIndex;
RT_L2_DATA int32_t rL1BufferIndex = 0;
RT_L1_DATA int16_t ecg_L1buff[DIM*(NLEADS+1)]; 
RT_L2_DATA int32_t start_index_buff = 0;
RT_L2_DATA int32_t end_main_loop;
RT_L2_DATA int32_t flag_prev_error = 0;
RT_L2_DATA int32_t* argCL[7];
RT_L1_DATA int32_t* indicesRpeaksCL[H_B+1];
RT_L2_DATA int32_t rpeaks_counter_cl;
RT_L2_DATA rt_event_sched_t * psched = 0;
RT_L2_DATA int32_t done = 0;
#endif

RT_L2_DATA int32_t overlap;
RT_L2_DATA int32_t rWindow;


void clearRelEn() {
    clearAndResetRelEn();
    resetPeakDetection();
}

#ifdef MODULE_CLUSTERING

static void cluster_Rpeaks(int32_t *arg[])
{
  rt_team_fork(NUM_CORES, rpeaks, arg);
}

extern void end_of_call(void *arg)
{
  done = 1;
}

static void fCore0_DmaTransfer_Window(void *arg)
{
  rt_dma_copy_t dmaCp;
      
    // Copy data block window from L2 to shared L1 memory using the cluster DMA
    rt_dma_memcpy((unsigned int)&ecg_L2buff[rL2BufferIndex], (unsigned int)&ecg_L1buff[rL1BufferIndex], 2*DIM, RT_DMA_DIR_EXT2LOC, 0, &dmaCp);

    // Wait for dma to finish
    rt_dma_wait(&dmaCp); 
}

#endif //#ifdef MODULE_CLUSTERING

void adaptiveRpeakDetection(){

#ifdef MODULE_CLUSTERING
    rt_cluster_mount(MOUNT, 0, 0, NULL);
#endif   

    int32_t count_sample = 0;
    int32_t offset_window = 0;
    int32_t offset_ind = LONG_WINDOW/2+1;
    int32_t tot_overlap = 0;

    overlap = 0;

#ifdef MODULE_MF
    int32_t flagMF = 0;
    int32_t i_lead = 0;

    buffSize_windowMF = LONG_WINDOW+DIM;

    argMF[0] = (int32_t*) ecg_L2buff;
    argMF[1] = &flagMF;
    argMF[2] = &i_lead;
    argMF[3] = &buffSize_windowMF;    

    init_filtering();
#endif

#ifdef MODULE_RELEN    
    buffSize_windowRelEn = LONG_WINDOW+DIM;

    argRelEn[0] = (int32_t*) ecg_L2buff;
    argRelEn[1] = (int32_t*) &ecg_L2buff[(LONG_WINDOW+DIM)*NLEADS];
    argRelEn[2] = &start_RelEn;
    argRelEn[3] = &buffSize_windowRelEn;

    clearRelEn();
#endif

#ifdef MODULE_RPEAK_REWARD 
    int32_t rpeaks_counter = 0;

    argRW_Rpeak[0] = (int32_t*) &ecg_L2buff[LONG_WINDOW+(LONG_WINDOW + DIM)*NLEADS];
    argRW_Rpeak[1] = indicesRpeaks;
    argRW_Rpeak[2] = &offset_ind;
#endif

#ifdef MODULE_CLUSTERING
    rL2BufferIndex = 0;

    // Allocate event on the default scheduler
    if (rt_event_alloc(NULL, 1)) return -1;
    rt_event_t *event;

#endif    

    for(rWindow=0; rWindow<N_WINDOWS; rWindow++)
    {
        if((rWindow+1)*DIM + LONG_WINDOW -1 - tot_overlap >= ECG_VECTOR_SIZE){
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
                ecg_L2buff[i+offset_window] = ecg_L2buff[(LONG_WINDOW + DIM - overlap + i + offset_window)];
#endif
#ifdef OVERLAP_RELEN
                ecg_L2buff[i + offset_window + (LONG_WINDOW + DIM)*NLEADS] = ecg_L2buff[(2*(LONG_WINDOW + DIM)*NLEADS - overlap + i + offset_window)];
#endif                
            }
        }

        for(int32_t lead=0; lead<NLEADS; lead++) {
            for(int32_t i=overlap + offset_window; i< LONG_WINDOW + DIM; i++) {
                ecg_L2buff[i + (DIM+LONG_WINDOW)*lead] = ecg_1l[rWindow*DIM + i - tot_overlap]; 
            }
        }

#ifdef PRINT_DEBUG_WINDOW
        printf("start_window: %d end_window: %d overlap: %d\n", offset_window + rWindow*DIM - tot_overlap,LONG_WINDOW -1 + (rWindow+1)*DIM - tot_overlap,overlap);
#endif

#ifdef MODULE_MF

        if (rWindow == 0) {
            // Needed to initialize the MF filter properly
            for(int32_t lead=0; lead<NLEADS; lead++) {
                for(int32_t i=0; i<=OFFSET_MF; i++) {
                    ecg_L2buff[i + (DIM+LONG_WINDOW)*lead] = 0;
                }
            }
        }

        // MF on FC because not enough memory on Cluster
        argMF[0] = (int32_t*) &ecg_L2buff[offset_window + overlap];
        buffSize_windowMF = LONG_WINDOW - offset_window + DIM-overlap;

        for(i_lead = 0; i_lead < NLEADS; i_lead++){
            filterWindows(argMF);
        }

        flagMF=1;

    #ifdef PRINT_SIG_MF
        for(int32_t sample = offset_window; sample<LONG_WINDOW + DIM; sample++) {
            printf("%d\n", ecg_L2buff[sample]);
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

        argRelEn[0] = (int32_t*) &ecg_L2buff[offset_window + overlap];
        argRelEn[1] = (int32_t*) &ecg_L2buff[offset_window + (LONG_WINDOW + DIM)*NLEADS+overlap];
        buffSize_windowRelEn = LONG_WINDOW - offset_window + DIM-overlap;
    #endif              

        relEn_w(argRelEn);       

    #ifdef HWPERF_MODULES
        profile_stop(perf);
    #endif

    #ifdef PRINT_RELEN
        for(int32_t sample = offset_window + (LONG_WINDOW + DIM)*NLEADS; sample< (LONG_WINDOW + DIM)*(NLEADS+1); sample++) {
            printf("%d\n", ecg_L2buff[sample]);
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
            printf("%d\n", indicesRpeaks[indR]);
        }
    #endif

#endif        

#ifdef MODULE_ERROR_DETECTION

    #ifdef MODULE_CLUSTERING
        if(error_RWindow == 0){
            start_index_buff = 0;
        }else{
            start_index_buff = DIM;
        }
    #endif
        argErrDet[0] = &rpeaks_counter;
        argErrDet[1] = &rWindow; 
        argErrDet[2] = &lastRR;        
        argErrDet[3] = indicesRpeaks;
        argErrDet[4] = &lastRpeak;
        error_RWindow = errorDetection(argErrDet);

    #ifdef PRINT_ERROR_RPEAKS
        printf("%d\n", error_RWindow);
    #endif
#endif        

#ifdef MODULE_CLUSTERING   

        rL2BufferIndex = LONG_WINDOW+(LONG_WINDOW + DIM)*NLEADS;
        end_main_loop = DIM*(NLEADS+1);

    #ifdef MODULE_ERROR_DETECTION    
        if(error_RWindow == 0 || rWindow == 0){ //Previous window or first window
            rL1BufferIndex = 0;
            // ----------------------------Copy previous window ecg buffer from L2 to L1 memory if error was 0 ------------------------------ //
            // Initialize event
            event = rt_event_get_blocking(NULL);

            // Run function on Core 0 of the cluster
            rt_cluster_call(NULL, 0, fCore0_DmaTransfer_Window, NULL, NULL, STACK_SIZE, STACK_SIZE, 1, event);

            // Wait for event
            rt_event_wait(event);
            // ------------------------------------------------------------------------------------------------------------------------------//
            flag_prev_error = 0;
        }else{
            rL1BufferIndex = DIM;
        }
        if(rWindow > 0 && error_RWindow == 1){
    #else
        if(rWindow == 0){
            rL1BufferIndex = 0;
            start_index_buff = 0;
            // -------------------------------- Copy first window ecg buffer from L2 to L1 memory ------------------------------------------ //
            // Initialize event
            event = rt_event_get_blocking(NULL);

            // Run function on Core 0 of the cluster
            rt_cluster_call(NULL, 0, fCore0_DmaTransfer_Window, NULL, NULL, STACK_SIZE, STACK_SIZE, 1, event);

            // Wait for event
            rt_event_wait(event);
            // ------------------------------------------------------------------------------------------------------------------------------//
        }else{
            rL1BufferIndex = DIM;
            start_index_buff = DIM;
    #endif       
            // ----------------Copy current window ecg buffer from L2 to L1 memory if error was 1 or after first window--------------------- //
            // Initialize event
            event = rt_event_get_blocking(NULL);

            // Run function on Core 0 of the cluster
            rt_cluster_call(NULL, 0, fCore0_DmaTransfer_Window, NULL, NULL, STACK_SIZE, STACK_SIZE, 1, event);

            // Wait for event
            rt_event_wait(event);
            // ------------------------------------------------------------------------------------------------------------------------------//
            argCL[0] = (int32_t*) ecg_L1buff; //The start index is the one in argCL[1]
            argCL[1] = &rWindow;
            argCL[2] = &start_index_buff;
            argCL[3] = &end_main_loop;
            argCL[4] = &offset_ind;
            argCL[5] = &flag_prev_error;
            argCL[6] = indicesRpeaksCL;
            rt_cluster_call(NULL, CID, cluster_Rpeaks, argCL, NULL, STACK_SIZE, STACK_SIZE, NUM_CORES, rt_event_get(psched, end_of_call, (void *) CID));
            while(!done)
                rt_event_execute(psched, 1);
            done = 0;

            // Move last sample to index DIM-1 of L1 buffer for the next window (this is necessary if the clustering module runs without error detection 
            // or if the previous error was 1 and the clustering must keep running)
            // The clustering uses the approximated signal derivative (a diff function), so it needs the last sample of the previous window to not miss anything
            ecg_L1buff[DIM-1] = ecg_L1buff[end_main_loop];

            flag_prev_error = 1;

    #ifdef PRINT_RPEAKS_CL
            rpeaks_counter_cl = 0;

            while(indicesRpeaksCL[rpeaks_counter_cl]!=0) {
                rpeaks_counter_cl++;
            }

            for(int ix_rr=0; ix_rr<rpeaks_counter_cl; ix_rr++){
                printf("%d\n",indicesRpeaksCL[ix_rr]);
            }
    #endif            
        }
#endif

#ifdef ONLY_FIRST_WINDOW //Only for debug
    return;
#endif

#ifndef MODULE_CLUSTERING
    #ifdef OVERLAP_MF    
        overlap = LONG_WINDOW + LONG_WINDOW/2 + 1;            
    #endif

    #ifdef OVERLAP_RELEN
        overlap = 0;
    #endif

        tot_overlap += overlap;      
#endif
        offset_ind = offset_ind + DIM - tot_overlap;  

#ifdef MODULE_RPEAK_REWARD        
        rpeaks_counter = 0;

        for(int32_t ix_rp = 0; ix_rp < H_B+1 ; ix_rp++) {
           indicesRpeaks[ix_rp] = 0;
        }
#endif        
    }

#ifdef MODULE_CLUSTERING
    rt_cluster_mount(UNMOUNT, 0, 0, NULL);
#endif    

}
