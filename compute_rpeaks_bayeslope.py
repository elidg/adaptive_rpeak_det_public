# Script to run BayeSlope on provided data (longer segment extracted from "ECG in High Intensity Exercise Dataset": https://doi.org/10.5281/zenodo.5727800)

import os.path
import subprocess
import shutil
import re
from natsort import natsorted 
import pandas as pd
import numpy as np
import change_code_inputs as cc
import matplotlib.pyplot as plt

folderCode = './C/'
folderData = './C/data/input_signals/'
folderResults = './Results/'

os.makedirs(folderResults, exist_ok=True) # create folder if it does not exist

df_errdet = pd.read_csv(folderCode + 'error_detection/input_errdet.csv')
input_files = natsorted([f for f in os.listdir(folderData)]) # ['sub3_seg4.csv']: one of the best cases; ['sub7_seg4.csv']: one of the worst cases
fs = 250 # Hz, sampling frequency of the signals in the analyzed dataset and validated in the R peak detection code
flagBS = 1 # if 1, run only BayeSlope; if 0, run full Adaptive R Peak Detection
if flagBS == 1:
	str_alg = 'bs_'
else:
	str_alg = 'adaptive_'
flagOutput = 'Rpeaks' # default 'Rpeaks'; 'MF': to print morphofological filtered signal; 'RELEN': to print relative energy signal 
offset_relen = 118 # samples considering 250 Hz, LONG_WINDOW/2-1 of Relative Energy algorithm

for file in input_files:

	filename = file.split('.')[0]
	print("Analyzing input " + filename + "...")

	match_sub = re.search(r'sub(\d+)', filename)
	subject = int(match_sub.group(1))
	match_seg = re.search(r'seg(\d+)', filename)
	segment = int(match_seg.group(1))

	# Load subject-specific percentiles for error detection
	perc_loo_low = df_errdet['LowPercentile'].values[df_errdet['Subject'] == subject][0]
	perc_loo_high = df_errdet['HighPercentile'].values[df_errdet['Subject'] == subject][0]

	# Load ecg signal from data folder
	df_sig = pd.read_csv(folderData + file, header=None)
	ecg_raw = df_sig.values[:,0]

	# Change input signal and error detection parameters in C code
	cc.change_input_signal_C(ecg_raw)
	cc.change_input_error_detection_C(subject,perc_loo_low,perc_loo_high)

	# Output morphological filtered signal
	cc.change_output('MF')
	subprocess.call(["make", "clean", "all","run"], cwd=folderCode) 
	# Copy output file to results folder
	shutil.copy(folderCode + "output.txt",folderResults + "mf_" + filename + ".txt") 
	print("Copied file in " + folderResults + "mf_" + filename + ".txt")
	df_mf = pd.read_csv(folderResults + "mf_" + filename + ".txt", header=None)
	ecg_mf = df_mf.values[:,0]
	t = np.linspace(0,len(ecg_mf)/fs,len(ecg_mf))

	# Output relative energy signal
	cc.change_output('RELEN')
	subprocess.call(["make", "clean", "all","run"], cwd=folderCode) 
	# Copy output file to results folder
	shutil.copy(folderCode + "output.txt",folderResults + "relen_" + filename + ".txt") 
	print("Copied file in " + folderResults + "relen_" + filename + ".txt")
	df_relen = pd.read_csv(folderResults + "relen_" + filename + ".txt", header=None)
	ecg_relen = df_relen.values[:,0]	

	# Run the C code for the R peak detection
	cc.change_algorithm(flagBS)
	cc.change_output('Rpeaks')
	subprocess.call(["make", "clean", "all","run"], cwd=folderCode) 
	# Copy output file to results folder
	shutil.copy(folderCode + "output.txt",folderResults + "rpeaks_" + str_alg + filename + ".txt") 
	print("Copied file in " + folderResults + "rpeaks_" + str_alg + filename + ".txt")
	if flagBS == 1:
		df_rpks = pd.read_csv(folderResults + "rpeaks_" + str_alg + filename + ".txt", header=None)
		rpks = df_rpks.values[:,0]
	else:
		df_rpks = pd.read_csv(folderResults + "rpeaks_" + str_alg + filename + ".txt", sep = ':', header=None)

	if len(input_files) == 1:
		plt.figure(figsize=(10, 6))
		plt.plot(t,ecg_mf, zorder = -1,label="ECG processed with Morphological Filtering")
		plt.scatter(t[rpks],ecg_mf[rpks], marker = "o",color = "red", zorder = 1,label = "R peaks detected via BayeSlope")
		plt.legend(loc="upper right")
		plt.tight_layout()
		plt.title("R peaks detection via BayeSlope - Subject " + str(subject) + " Segment" + str(segment))
		plt.xlabel("Time (s)")
		plt.ylabel("Amplitude (mV)")

		plt.figure(figsize=(10, 6))
		plt.plot(t,ecg_relen, zorder = -1,label="ECG processed with Relative Energy")
		plt.scatter(t[rpks+offset_relen],ecg_relen[rpks+offset_relen], marker = "o",color = "red", zorder = 1,label = "R peaks detected via BayeSlope")
		plt.legend(loc="upper right")
		plt.tight_layout()
		plt.title("R peaks detection via BayeSlope - Subject " + str(subject) + " Segment" + str(segment))
		plt.xlabel("Time (s)")
		plt.ylabel("Amplitude (mV)")