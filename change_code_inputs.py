# Library to change inputs, parameters and outputs of C code

folderCode = './C/'

# Change input signal in C code
def change_input_signal_C(sig):
	# read a list of lines into data
	with open(folderCode + 'data/signal.h', 'r') as file:
		data = file.readlines()

	# Change the input signal and related parameters in signal.h
	data[6] = "#define ECG_VECTOR_SIZE " + str(len(sig)) + "\n"
	data[8] = "int16_t ecg_1l[ECG_VECTOR_SIZE] = {"
	ecgstr = ''.join([str(int(i))+"," for i in sig])
	data[8] = ''.join([data[8],ecgstr,"};\n"])

	# and write everything back
	with open(folderCode + 'data/signal.h', 'w') as file:
		file.writelines( data )

# Change error detection input (percentile of Leave-One-Subject-Out distribution)
def change_input_error_detection_C(subject, perc_low, perc_high):
	# read a list of lines into data
	with open(folderCode + 'error_detection/error_detection.h', 'r') as file:
		data = file.readlines()

	# Change the subject-specific percentiles in error_detection.h
	data[10] = "#define PERCENTILE_LOO_LOW " + str(perc_low) + "\n"
	data[11] = "#define PERCENTILE_LOO_HIGH " + str(perc_high) + "\n"

	# and write everything back
	with open(folderCode + 'error_detection/error_detection.h', 'w') as file:
		file.writelines( data )

# Change defines.h to run only BayeSlope or full Adaptive R Peak Detection
def change_algorithm(flagBS):
	# read a list of lines into data
	with open(folderCode + 'defines.h', 'r') as file:
		data = file.readlines()

	# Comment or uncomment module of REWARD algorithm to run only BayeSlope or full Adaptive R Peak Detection
	if flagBS == 1:
		data[57] = "// #define MODULE_RPEAK_REWARD" + "\n"
	else:
		data[57] = "#define MODULE_RPEAK_REWARD" + "\n"

	# and write everything back
	with open(folderCode + 'defines.h', 'w') as file:
		file.writelines( data )

# Change defines.h to run only BayeSlope or full Adaptive R Peak Detection
def change_output(flagOutput):
	# read a list of lines into data
	with open(folderCode + 'defines.h', 'r') as file:
		data = file.readlines()

	# Comment or uncomment module of REWARD algorithm to run only BayeSlope or full Adaptive R Peak Detection
	if flagOutput == 'MF':
		data[64] = "// #define ACCURACY" + "\n"
		data[67] = "#define PRINT_SIG_MF" + "\n"
		data[68] = "// #define PRINT_RELEN" + "\n"
	elif flagOutput == 'RELEN':
		data[64] = "// #define ACCURACY" + "\n"
		data[67] = "// #define PRINT_SIG_MF" + "\n"
		data[68] = "#define PRINT_RELEN" + "\n"
	else:
		data[64] = "#define ACCURACY" + "\n"
		data[67] = "// #define PRINT_SIG_MF" + "\n"
		data[68] = "// #define PRINT_RELEN" + "\n"

	# and write everything back
	with open(folderCode + 'defines.h', 'w') as file:
		file.writelines( data ) 