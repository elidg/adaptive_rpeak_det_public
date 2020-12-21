PULP_APP = adaptive_rpeak_detection

PULP_APP_FC_SRCS  = main.c \
				       adaptive_Rpeak_detection/adaptiveRpeakDetection.c \
				       Morph_filt/morpho_filtering.c \
				       REWARD_R_peak_detection/relativeEnergy.c \
				       REWARD_R_peak_detection/peakDetection.c \
				       error_detection/error_detection.c \
				       profiling/profile.c \
				       profiling/profile_cl.c  

PULP_APP_CL_OMP_SRCS  = kmeans_clustering/kmean/k_mean_functions.c \
						kmeans_clustering/search/search_functions.c

CORES  ?= 1
CTARGET ?= 1

stackSize ?= 2048

PULP_CFLAGS += -DTARGET=$(CTARGET) -DNUM_CORES=$(CORES) -O3 -g3 -w
PULP_LDFLAGS = -lm #-lg #uncomment -lg for pulpissimo

include $(PULP_SDK_HOME)/install/rules/pulp_rt.mk
