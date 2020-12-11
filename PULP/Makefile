PULP_APP = delineation_app

PULP_APP_FC_SRCS  = main.c \
				       adaptive_Rpeak_detection/adaptiveRpeakDetection.c \
				       Morph_filt/morpho_filtering.c \
				       REWARD_R_peak_detection/relativeEnergy.c \
				       REWARD_R_peak_detection/peakDetection.c \
				       profiling/profile.c  	
CORES  ?= 1
CTARGET ?= 0

PULP_CFLAGS += -DTARGET=$(CTARGET) -DNUM_CORES=$(CORES) -O3 -g3 -w
PULP_LDFLAGS = -lm #-lg #uncomment -lg for pulpissimo

include $(PULP_SDK_HOME)/install/rules/pulp_rt.mk
