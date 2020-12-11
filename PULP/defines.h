#ifndef DEFINES_H_
#define DEFINES_H_ 

//========= DEFINE PRINT OUTPUT ========//
// #define PRINT_SIG_MF
// #define PRINT_RELEN
// #define PRINT_RPEAKS
#define PRINT_ERROR_RPEAKS

//========= DEFINE PRINT DEBUG =========//
// #define PRINT_DEBUG
// #define PRINT_SIG_INPUT_PEAKS
// #define PRINT_THR
// #define PRINT_RPEAKS_BEFORE_T_CHECK

//=========== Sampling frequency ========== //
#define ECG_SAMPLING_FREQUENCY 250

//=========== Number of leads/channels ============//
#define NLEADS 1

//=========== Buffer size for R peak detection ========== //
//The minumum distance between two ECG peaks based on physiological limits: 200 miliseconds * sampling frequency
#define BUFFER_LENGTH_SECONDS (float) 1.75 //2.048 //
#define BUFFER_SIZE (int16_t) (BUFFER_LENGTH_SECONDS*ECG_SAMPLING_FREQUENCY)

//=========== Overlap for Rel-En ========== //
//long window delay length (in samples): sampling frequency times long window length in seconds
//Long window length = 0.95 seconds
#define LONG_WINDOW (uint16_t) (0.95*ECG_SAMPLING_FREQUENCY+1)

//=========== R PEAK DETECTION PARAMETERS ===========//
// #define NEGATIVE_PEAK //For ISSUL-EPFL dataset, comment this line. For QTDB used for publication in EMBC 2019 and ESWEEK 2020, uncomment.

//======== DEFINE MODULES ==========//
#define MODULE_MF
#define MODULE_RELEN
#define MODULE_RPEAK_REWARD
#ifdef MODULE_RPEAK_REWARD
#define MODULE_ERROR_DETECTION
#endif
//======== DEFINE WINDOW OVERLAP ==========//
// Copy the MF or the RelEn signal in the overlapping window. The two defines are mutually exclusive (only one can be uncommented)
// #define OVERLAP_MF 
#define OVERLAP_RELEN

//======== DEFINE WINDOWS (ONLY FOR DEBUG) ==========//
// #define ONLY_FIRST_WINDOW

#define N 1
#define H_B 30
#define dim  (int16_t) (BUFFER_SIZE * N) //((BUFFER_SIZE * N) + LONG_WINDOW)
#if ECG_SAMPLING_FREQUENCY == 250
#define OFFSET_MF 150
#else
#define OFFSET_MF 300
#endif
#define FACTOR_MS 1000
#define FACTOR_RATIO_RR 1000
#define PERCENTILE_LOO_LOW 641 //This is the 0.5th percentile of the RR distribution detected by the clustering using the leave-one-out approach for each subject (it is subject-specific)
#define PERCENTILE_LOO_HIGH 1492 //This is the 99.5th percentile of the RR distribution detected by the clustering using the leave-one-out approach for each subject (it is subject-specific)


//========= DEFINE PROFILING ================//
// #define HWPERF_MODULE_RPEAK_REWARD	//start profiling module RPEAK
// #define HWPERF_MODULES	//start profiling separate modules
// #define HWPERF_FULL	//start profiling full app (N.B. it will profile also some buffering)
// #define ACTIVE
//#define EXTACC		//# of loads and stores in EXT memory (L2)
//#define INTACC		//# of loads and stores in INT memory (L1)
//#define STALL			//# number of core stalls
//#define INSTRUCTION	//# number of instructions
//#define TCDM			//# of conflicts in TCDM (L1 memory) between cores 

#define PULP_L1_DATA RT_L2_DATA
#define PULP_L2_DATA RT_L2_DATA

#define MUL 
#define SCALE 100//6//64

#define CGRA_OFF

#endif // DEFINES_H_
