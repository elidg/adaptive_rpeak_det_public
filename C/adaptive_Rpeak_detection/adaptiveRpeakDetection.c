#include "adaptiveRpeakDetection.h"
#include "../defines.h"
#include "../data/signal.h"
#include "../Morph_filt/morpho_filtering.h"
#include "../Morph_filt/defines_globals.h"
#include "../error_detection/error_detection.h"
#include "../kmeans_clustering/kmean/k_mean_functions.h"

#define N_WINDOWS (int) (2*((ECG_VECTOR_SIZE-LONG_WINDOW)/DIM)+1)// Counting the worst case scenario when overlap is DIM
 
int16_t ecg_L2buff[(LONG_WINDOW+DIM)*(NLEADS+1)];

#ifdef MODULE_MF
int32_t *argMF[4];
int32_t buffSize_windowMF;
#endif

#ifdef MODULE_RELEN
int32_t *argRelEn[4];
int32_t start_RelEn = 1; 
int32_t buffSize_windowRelEn;
#endif

#ifdef MODULE_RPEAK_REWARD
int32_t *argRW_Rpeak[3];
int32_t indicesRpeaks[H_B+1];
#endif

#ifdef MODULE_ERROR_DETECTION
int32_t *argErrDet[5];
int32_t lastRpeak = 0;
int32_t lastRR = 0;
int32_t error_RWindow = 0;
#endif

#ifdef MODULE_CLUSTERING
int32_t rL2BufferIndex = 0;
int32_t rL1BufferIndex = 0;
int32_t sizeBufferTransfer = DIM; // In PULP version this is in bytes, here is in elements
int16_t ecg_L1buff[DIM*(NLEADS+1)]; 
int32_t start_index_buff = 0;
int32_t end_main_loop;
int32_t flag_prev_error = 0;
int32_t* argCL[9];
int32_t overlapCL = 1;
int32_t indicesRpeaksCL[H_B+1];
int32_t overlap_qrs_init = 0;
int32_t rpeaks_counter_cl;
int32_t done = 0;
#endif

int32_t overlap;
int32_t rWindow;


void clearRelEn() {
    clearAndResetRelEn();
    resetPeakDetection();
}

#ifdef MODULE_CLUSTERING

void copy_fromL2toL1_memory()
{
    // This function is only a transfer simulation of a buffer segment from a virtual L2 (core) to L1 (cluster) memory
    for(int32_t i=0; i<sizeBufferTransfer; i++) {
        ecg_L1buff[rL1BufferIndex+i] = ecg_L2buff[rL2BufferIndex+i];
    }
}

#endif

void adaptiveRpeakDetection(){

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
        printf("start_window: %d end_window: %d overlapCL: %d offset_ind: %d\n", offset_window + rWindow*DIM - tot_overlap,LONG_WINDOW -1 + (rWindow+1)*DIM - tot_overlap,overlapCL, offset_ind);
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


#ifdef MODULE_RELEN

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

    #ifdef PRINT_RELEN
        for(int32_t sample = offset_window + (LONG_WINDOW + DIM)*NLEADS; sample< (LONG_WINDOW + DIM)*(NLEADS+1); sample++) {
            printf("%d\n", ecg_L2buff[sample]);
        }
    #endif 

#endif

#ifdef MODULE_RPEAK_REWARD

        getPeaks_w(argRW_Rpeak);

        rpeaks_counter = 0;

        while(indicesRpeaks[rpeaks_counter]!=0) {
            rpeaks_counter++;
        }

    #ifdef PRINT_RPEAKS
        for(int32_t indR=0; indR<rpeaks_counter; indR++) {
            printf("RW: %d\n", indicesRpeaks[indR]);
        }
    #endif

#endif        

#ifdef MODULE_ERROR_DETECTION

    #ifdef MODULE_CLUSTERING
        if(error_RWindow == 0 || rWindow == 0){
            start_index_buff = 0;
        }else if(error_RWindow == 1 && rWindow>1){
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
        printf("Err: %d\n", error_RWindow);
    #endif
#endif        

#ifdef MODULE_CLUSTERING   

        rL2BufferIndex = LONG_WINDOW+(LONG_WINDOW + DIM)*NLEADS;
        end_main_loop = DIM*(NLEADS+1);

    #ifdef MODULE_ERROR_DETECTION    
        if(error_RWindow == 0 || rWindow == 0){ //Previous window or first window

            // ======= Find the r peaks at the border of CL if it was run in the previous window (err=1) and the peak was not finished ======= //
            if(flag_prev_error==1 && overlap_qrs_init!=0){
                rL1BufferIndex = DIM;
                sizeBufferTransfer = (MAX_QRS_DUR + 2 - overlap_qrs_init);
                // --------------- Copy the few samples to finish the peak from the previous window from L2 to L1 memory if error was 1 ----------------- //
                copy_fromL2toL1_memory();
                // ------------------------------------------------------------------------------------------------------------------------------------- //   

                start_index_buff = DIM;
                end_main_loop = DIM + (MAX_QRS_DUR + 2 - overlap_qrs_init);

                argCL[0] = (int32_t*) ecg_L1buff; //The start index is the one in argCL[2]
                argCL[1] = &rWindow;
                argCL[2] = &start_index_buff;
                argCL[3] = &end_main_loop;
                argCL[4] = &offset_ind;
                argCL[5] = &flag_prev_error;
                argCL[6] = indicesRpeaksCL;
                argCL[7] = &overlapCL;
                argCL[8] = &overlap_qrs_init;

                rpeaks(argCL);

    #ifdef PRINT_RPEAKS_CL
                rpeaks_counter_cl = 0;

                while(indicesRpeaksCL[rpeaks_counter_cl]!=0) {
                    rpeaks_counter_cl++;
                }

                for(int ix_rr=0; ix_rr<rpeaks_counter_cl; ix_rr++){
                    printf("CL: %d\n",indicesRpeaksCL[ix_rr]);
                }
    #endif           
                //Reset start and end variables in case the CL has to run in this window
                start_index_buff = 0;
                end_main_loop = DIM*(NLEADS+1);
                sizeBufferTransfer = DIM;
                overlap_qrs_init = 0;
            }
            // ======================================================================================================================= //

            rL1BufferIndex = 0;
            // ----------------------------Copy previous window ecg buffer from L2 to L1 memory if error was 0 ------------------------------ //
            copy_fromL2toL1_memory();
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
            copy_fromL2toL1_memory();
            // ------------------------------------------------------------------------------------------------------------------------------//
        }else{
            rL1BufferIndex = DIM;
            if(rWindow == 1)
                start_index_buff = 0;
            else
                start_index_buff = DIM;
    #endif       
            // ----------------Copy current window ecg buffer from L2 to L1 memory if error was 1 or after first window--------------------- //
            copy_fromL2toL1_memory();
            // ------------------------------------------------------------------------------------------------------------------------------//
            argCL[0] = (int32_t*) ecg_L1buff; //The start index is the one in argCL[2]
            argCL[1] = &rWindow;
            argCL[2] = &start_index_buff;
            argCL[3] = &end_main_loop;
            argCL[4] = &offset_ind;
            argCL[5] = &flag_prev_error;
            argCL[6] = indicesRpeaksCL;
            argCL[7] = &overlapCL;
            argCL[8] = &overlap_qrs_init;

            rpeaks(argCL);           
            // Move last sample to index DIM-1 of L1 buffer for the next window (this is necessary if the clustering module runs without error detection 
            // or if the previous error was 1 and the clustering must keep running)
            // The clustering uses the approximated signal derivative (a diff function), so it needs the last sample of the previous window to not miss anything
            for(int ix_ov = 0; ix_ov < overlapCL; ix_ov++){
                ecg_L1buff[DIM-overlapCL+ix_ov] = ecg_L1buff[end_main_loop-overlapCL+ix_ov];
            }

            flag_prev_error = 1;

    #ifdef PRINT_RPEAKS_CL
            rpeaks_counter_cl = 0;

            while(indicesRpeaksCL[rpeaks_counter_cl]!=0) {
                rpeaks_counter_cl++;
            }

            for(int ix_rr=0; ix_rr<rpeaks_counter_cl; ix_rr++){
                #ifdef MODULE_RPEAK_REWARD
                    printf("BS: %d\n",indicesRpeaksCL[ix_rr]);
                #else
                    printf("%d\n",indicesRpeaksCL[ix_rr]);
                #endif
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

}
