CC         = gcc
#SRC        = 

TARGET = ./bin/adaptive_rpeak_detection

CFLAGS  = -O3
LDFLAGS = -lm -lg 
#VPATH      = $(SRC)

SOURCES    = main.c \
				adaptive_Rpeak_detection/adaptiveRpeakDetection.c \
				Morph_filt/morpho_filtering.c \
				REWARD_R_peak_detection/relativeEnergy.c \
				REWARD_R_peak_detection/peakDetection.c \
				error_detection/error_detection.c \
				kmeans_clustering/kmean/k_mean_functions.c \
				kmeans_clustering/search/search_functions.c				

				 
OBJECTS    = $(SOURCES:%/%.c=%/%.o)
 

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $+ $(LDFLAGS)
	#cp -f $(TARGET) /bin/

%.o: %.c
	$(CC) $(CFLAGS) -c $<

# run: $(TARGET)
# 	$(TARGET) ecg.txt > output.txt

# --- Uncomment the following two lines if you are using a static array as input signal
run: $(TARGET)
	./$(TARGET) > output.txt
# ---

clean:
	rm -rf ./bin
	mkdir ./bin

