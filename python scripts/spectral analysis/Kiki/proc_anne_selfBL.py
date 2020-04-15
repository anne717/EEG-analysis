# A script to analyse spectral data from exported SleepSign files

# Make sure Python is looking in the correct directory for the file
import os
os.chdir('E:/spectral analysis/jan2020_ott')

# Open SleepSign exported file and extract data in a sensible fashion
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mb
import numpy.ma as ma

filename = 'j1_sham_fft.txt'
file_output = filename.replace('.txt','')
export_folder = 'E:/spectral analysis/jan2020_ott/results/'
indx_freq_start = 1
indx_freq_end = -1
hz60_start = 58
hz60_end = 62
epoch_avg_number =500

file_extras = np.genfromtxt(filename,skip_header = 19,dtype = [('epoch','float'),('stage','|S1'),('date','|S11'),('time','|S11')])
whole_file = np.genfromtxt(filename,skip_header = 19);
data = whole_file[:,4:];
stage = file_extras['stage'];
file_freq=np.genfromtxt(filename,skip_header = 18,invalid_raise =False, dtype = 'string')[3:];
frequency = np.char.strip(file_freq,  'Hz').astype(np.float)

stage_num = np.zeros(len(stage))

# Manipulate data so that each value is expressed as a percentage of the whole spectrum
data_trunc = data[:,indx_freq_start:indx_freq_end]    #remove first and last frequency band as it may contain artifacts
num_freq_values = np.size(data_trunc[0,:])
                    
# Get rid of 60 Hz noise  
use_freq = frequency[indx_freq_start:indx_freq_end]
noisy = np.logical_and(use_freq > hz60_start, use_freq < hz60_end)
noise = mb.find(noisy == True)

data_neg = data_trunc
data_neg[:,noise[0]:noise[-1]] = -1

masked_data = ma.masked_less(data_neg,0)
data_sum = np.sum(masked_data, axis = 1)
normalized_data = 100*(masked_data/data_sum[:,np.newaxis])
 
#stage_numlist = np.array(stage_num).tolist()
#num_wake = stage_numlist.count('W')    # number of wake episodes
#num_NREM = stage_numlist.count('N')    # number of NREM episodes
#num_REM = stage_numlist.count('R')    # number of REM episodes
#num_art = stage_numlist.count('M')    # number of artefacts

num_wake = np.sum(stage == 'W')
num_NREM = np.sum(stage == 'N')
num_REM = np.sum(stage == 'R')
num_art = np.sum(stage == 'M')

art_percent = float(num_art) / float(np.size(stage)) * 100

if num_wake < 1:
    num_wake = 1
    
if num_NREM < 1:
    num_NREM = 1

if num_REM < 1:
    num_REM = 1

wake = []
wake_indx = []
sws = []
sws_indx = []
rem = []
rem_indx = []

wake_indx = np.where(stage == 'W');
wake = np.zeros((num_wake, np.size(normalized_data[0,:])));
wake_mat = normalized_data[wake_indx[0],:]

sws_indx = np.where(stage == 'N')
sws = np.zeros((num_NREM, np.size(normalized_data[0,:])));
sws_mat = normalized_data[sws_indx[0],:] 

rem_indx = np.where(stage == 'R')
rem = np.zeros((num_REM, np.size(normalized_data[0,:])));
rem_mat = normalized_data[rem_indx[0],:]          
        
if np.size(rem_mat) < np.size(normalized_data[0,:]):
    rem = np.zeros(num_freq_values)

if np.size(sws_mat) < np.size(normalized_data[0,:]):
    sws = np.zeros(num_freq_values)  

if np.size(wake_mat) < np.size(normalized_data[0,:]):
    wake = np.zeros(num_freq_values)               

wake_mask = ma.masked_less(wake_mat, 0)
sws_mask = ma.masked_less(sws_mat, 0)
rem_mask = ma.masked_less(rem_mat, 0)

wake_mean = ma.mean(wake_mask, axis = 0)      
sws_mean = ma.mean(sws_mask, axis = 0) 
rem_mean = ma.mean(rem_mask, axis = 0) 

np.savetxt((str(export_folder) + str(file_output) +'_wake'), wake_mean)
np.savetxt((str(export_folder) + str(file_output) +'_sws'), sws_mean)
np.savetxt((str(export_folder) + str(file_output) +'_rem'), rem_mean)
np.savetxt((str(export_folder) + str(file_output) +'_frequency'), use_freq)

figure1 = plt.figure()
plt.hold(True)
line1, = plt.plot(use_freq, wake_mean, color = 'r')
line2, = plt.plot(use_freq, sws_mean, color = 'b')
line3, = plt.plot(use_freq, rem_mean, color = 'g')
label1 = 'wake'
label2 = 'sws'
label3 = 'rem'
plt.figlegend((line1, line2, line3),(label1, label2, label3), loc = 'upper right', fontsize = 10, frameon = False)

title = plt.title('raw power spectra: ' + file_output)
xlab = plt.xlabel('frequency (Hz)')
ylab = plt.ylabel('mean spectral power ($V^{2}$)')
xlim = plt.xlim(0,80)

plt.figtext(0.1, 0.02, "num wake = " + repr(num_wake) + ", num sws = " + repr(num_NREM) + ", num_rem = " + repr(num_REM) + ", num artefacts = " + repr(num_art) + ", % artefacts = ~" + repr(int(art_percent)), fontsize = 10)

plt.hold(False)