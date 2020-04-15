# -*- coding: utf-8 -*-
"""
Created on Mon Mar 04 14:06:34 2013

@author: avenner
"""
#timestamp = np.genfromtxt('timestamp_' + str(filename), dtype = [('date','|S10'),('time','|S9')])
#timestamp = np.loadtxt('timestamp_hdc104_cno10.txt', dtype = '|S19', delimiter = ',') 

# CALCULATES MEAN SPECTRA AGAINST AS A PERCENTAGE OF BASELINE AND THEN CALCULATES 
# BANDWIDTH BY AVERAGING OVER NUMBER OF FREQUENCY POINTS

import os
import numpy as np
import matplotlib
import matplotlib.mlab as mb
import datetime as dt
import matplotlib.dates as md
import matplotlib.pyplot as plt
import scipy.stats as st
import numpy.ma as ma

os.chdir('E:/spectral analysis/jan2020_ott/')

filename = 'j1_sham_fft.txt'
fig_title = 'Average: 10am CNO injection in SERT-M3 mice'
path = 'E:/spectral analysis/jan2020_ott/results/'  

stage1 = 'wake'
stage2 = 'sws'
stage3 = 'rem'

condition1 = 'sham'
condition2 = 'tbi1'

indx_freq_start = 1
indx_freq_end = -1
hz60_start = 58
hz60_end = 62

value_threshold = 0

delta = [0,4]
lotheta = [4,8]
hitheta = [8,12]
alph = [12,20]
gamma = [20,50] 
higamma = [70,200]

file_freq=np.genfromtxt(filename,skip_header = 18,invalid_raise =False, dtype = 'string')[3:];
frequency = np.char.strip(file_freq,  'Hz').astype(np.float)
use_freq = frequency[indx_freq_start:indx_freq_end] 

noisy = np.logical_and(use_freq > hz60_start, use_freq < hz60_end)
noise = mb.find(noisy == True)

dirList= os.listdir(path)    # insert the path to the directory of interest
filename = []
file_dict = {}

for fname in dirList:
    if fname.startswith('j'):
        filename = np.append(filename, fname)

fdata = list(filename)

for i in np.arange(filename.size):        
    fdata[i] = np.loadtxt(str(path)+str(filename[i]))
    
for i in np.arange(filename.size):
        file_dict[filename[i]] = fdata[i]
        
index_keys = np.arange(len(file_dict.keys()))
num_wake_sal = []
num_wake_cno = []
num_sws_sal = []
num_sws_cno = []
num_rem_sal = []
num_rem_cno = []
    
for i in np.arange(index_keys.size):
    if stage1 in file_dict.keys()[i] and condition1 in file_dict.keys()[i]:
        num_wake_sal = np.append(num_wake_sal, 1)
    else:
        if stage1 in file_dict.keys()[i] and condition2 in file_dict.keys()[i]:
            num_wake_cno = np.append(num_wake_cno, 1)
        else:
            if stage2 in file_dict.keys()[i] and condition1 in file_dict.keys()[i]:
               num_sws_sal = np.append(num_sws_sal, 1)
            else:
                if stage2 in file_dict.keys()[i] and condition2 in file_dict.keys()[i]:
                    num_sws_cno = np.append(num_sws_cno, 1)
                else:
                    if stage3 in file_dict.keys()[i] and condition1 in file_dict.keys()[i]:
                        num_rem_sal = np.append(num_rem_sal, 1)
                    else:
                        if stage3 in file_dict.keys()[i] and condition2 in file_dict.keys()[i]:
                            num_rem_cno = np.append(num_rem_cno, 1)
                            
arraylen = 0 
x = 0      

for i in np.arange(index_keys.size):
        x = np.size(file_dict.values()[i])
        if x > arraylen:
            arraylen = np.size(file_dict.values()[i])

null = -1
extra = []
app_values = [] 
       
for i in np.arange(index_keys.size):
    if arraylen > np.size(file_dict.values()[i]):
        extra = arraylen - np.size(file_dict.values()[i])
        app_values = np.tile(null, extra)
        file_dict[file_dict.keys()[i]] = np.append(file_dict.values()[i], app_values)
        
wake_sal_values = np.zeros((len(num_wake_sal), arraylen))
wake_sal_keys = range(0, len(num_wake_sal))
wake_cno_values = np.zeros((len(num_wake_cno), arraylen))
wake_cno_keys = range(0, len(num_wake_cno))
sws_sal_values = np.zeros((len(num_sws_sal), arraylen))
sws_sal_keys = range(0, len(num_sws_sal))
sws_cno_values = np.zeros((len(num_sws_cno), arraylen))
sws_cno_keys = range(0, len(num_sws_cno))
rem_sal_values = np.zeros((len(num_rem_sal), arraylen))
rem_sal_keys = range(0, len(num_rem_sal))
rem_cno_values = np.zeros((len(num_rem_cno), arraylen))
rem_cno_keys = range(0, len(num_rem_cno))

q = -1
p = -1
r = -1
s = -1
t = -1
u = -1

for i in np.arange(index_keys.size):
    if stage1 in file_dict.keys()[i] and condition1 in file_dict.keys()[i]:
        q = q + 1
        wake_sal_keys[q] = file_dict.keys()[i] 
        wake_sal_values[q,:] = file_dict.values()[i]
    else:
        if stage1 in file_dict.keys()[i] and condition2 in file_dict.keys()[i]:
            p = p + 1
            wake_cno_keys[p] = file_dict.keys()[i] 
            wake_cno_values[p,:] = file_dict.values()[i]
        else:
            if stage2 in file_dict.keys()[i] and condition1 in file_dict.keys()[i]:
                r = r + 1
                sws_sal_keys[r] = file_dict.keys()[i] 
                sws_sal_values[r,:] = file_dict.values()[i]
            else:
                if stage2 in file_dict.keys()[i] and condition2 in file_dict.keys()[i]:
                    s = s + 1
                    sws_cno_keys[s] = file_dict.keys()[i] 
                    sws_cno_values[s,:] = file_dict.values()[i]
                else:
                    if stage3 in file_dict.keys()[i] and condition1 in file_dict.keys()[i]:
                        t = t + 1
                        rem_sal_keys[t] = file_dict.keys()[i] 
                        rem_sal_values[t,:] = file_dict.values()[i]
                    else:
                        if stage3 in file_dict.keys()[i] and condition2 in file_dict.keys()[i]:
                            u = u + 1
                            rem_cno_keys[u] = file_dict.keys()[i] 
                            rem_cno_values[u,:] = file_dict.values()[i]

sorted_wake_sal_keys = np.sort(wake_sal_keys)
order_index_wake_sal = np.arange(len(num_wake_sal))
sorted_wake_cno_keys = np.sort(wake_cno_keys)
order_index_wake_cno = np.arange(len(num_wake_cno))

sorted_sws_sal_keys = np.sort(sws_sal_keys)
order_index_sws_sal = np.arange(len(num_sws_sal))
sorted_sws_cno_keys = np.sort(sws_cno_keys)
order_index_sws_cno = np.arange(len(num_sws_cno))

sorted_rem_sal_keys = np.sort(rem_sal_keys)
order_index_rem_sal = np.arange(len(num_rem_sal))
sorted_rem_cno_keys = np.sort(rem_cno_keys)
order_index_rem_cno = np.arange(len(num_rem_cno))
    
for i in np.arange(num_wake_sal.size):
    order_index_wake_sal[i] = mb.find(sorted_wake_sal_keys == wake_sal_keys[i])
    
for i in np.arange(num_wake_cno.size):
    order_index_wake_cno[i] = mb.find(sorted_wake_cno_keys == wake_cno_keys[i])

for i in np.arange(num_sws_sal.size):
    order_index_sws_sal[i] = mb.find(sorted_sws_sal_keys == sws_sal_keys[i])
    
for i in np.arange(num_sws_cno.size):
    order_index_sws_cno[i] = mb.find(sorted_sws_cno_keys == sws_cno_keys[i])

for i in np.arange(num_rem_sal.size):
    order_index_rem_sal[i] = mb.find(sorted_rem_sal_keys == rem_sal_keys[i])
    
for i in np.arange(num_rem_cno.size):
    order_index_rem_cno[i] = mb.find(sorted_rem_cno_keys == rem_cno_keys[i])          
        
sorted_wake_sal_values = np.zeros((len(wake_sal_keys), arraylen))
sorted_wake_cno_values = np.zeros((len(wake_cno_keys), arraylen))
sorted_sws_sal_values = np.zeros((len(sws_sal_keys), arraylen))
sorted_sws_cno_values = np.zeros((len(sws_cno_keys), arraylen))
sorted_rem_sal_values = np.zeros((len(rem_sal_keys), arraylen))
sorted_rem_cno_values = np.zeros((len(rem_cno_keys), arraylen))
    
for i in np.arange(num_wake_sal.size):
    sorted_wake_sal_values[order_index_wake_sal[i],:] = wake_sal_values[i,:]    
    
for i in np.arange(num_wake_cno.size):
    sorted_wake_cno_values[order_index_wake_cno[i],:] = wake_cno_values[i,:]
    
for i in np.arange(num_sws_sal.size):
    sorted_sws_sal_values[order_index_sws_sal[i],:] = sws_sal_values[i,:]    
    
for i in np.arange(num_sws_cno.size):
    sorted_sws_cno_values[order_index_sws_cno[i],:] = sws_cno_values[i,:]
    
for i in np.arange(num_rem_sal.size):
    sorted_rem_sal_values[order_index_rem_sal[i],:] = rem_sal_values[i,:]    
    
for i in np.arange(num_rem_cno.size):
    sorted_rem_cno_values[order_index_rem_cno[i],:] = rem_cno_values[i,:]
 
# Mask out 60 Hz noise and excess data points
wake_sal_neg = sorted_wake_sal_values
wake_cno_neg = sorted_wake_cno_values
sws_sal_neg = sorted_sws_sal_values
sws_cno_neg = sorted_sws_cno_values
rem_sal_neg = sorted_rem_sal_values
rem_cno_neg = sorted_rem_cno_values


wake_sal_neg[noise[0]:noise[-1]] = -1
wake_cno_neg[noise[0]:noise[-1]] = -1
sws_sal_neg[noise[0]:noise[-1]] = -1
sws_cno_neg[noise[0]:noise[-1]] = -1
rem_sal_neg[noise[0]:noise[-1]] = -1
rem_cno_neg[noise[0]:noise[-1]] = -1
   
masked_wake_sal_values = ma.masked_less(wake_sal_neg, value_threshold)   
masked_wake_cno_values = ma.masked_less(wake_cno_neg, value_threshold)

masked_sws_sal_values = ma.masked_less(sws_sal_neg, value_threshold)   
masked_sws_cno_values = ma.masked_less(sws_cno_neg, value_threshold)

masked_rem_sal_values = ma.masked_less(rem_sal_neg, value_threshold)   
masked_rem_cno_values = ma.masked_less(rem_cno_neg, value_threshold)

mean_wake_sal = ma.mean(masked_wake_sal_values, axis = 0)
mean_wake_cno = ma.mean(masked_wake_cno_values, axis = 0)
mean_sws_sal = ma.mean(masked_sws_sal_values, axis = 0)
mean_sws_cno = ma.mean(masked_sws_cno_values, axis = 0)
mean_rem_sal = ma.mean(masked_rem_sal_values, axis = 0)
mean_rem_cno = ma.mean(masked_rem_cno_values, axis = 0)

wakesal_sem_line = st.sem(masked_wake_sal_values, axis = 0)
wakecno_sem_line = st.sem(masked_wake_cno_values, axis = 0)
swssal_sem_line = st.sem(masked_sws_sal_values, axis = 0)
swscno_sem_line = st.sem(masked_sws_cno_values, axis = 0)
remsal_sem_line = st.sem(masked_rem_sal_values, axis = 0)
remcno_sem_line = st.sem(masked_rem_cno_values, axis = 0)
                
# Plot BL, saline and CNO spectra on the same figure
fig = plt.figure(facecolor = 'w')
ax = fig.add_subplot(111)

plt.hold(True)
ax1 = fig.add_subplot(311)
wakesal_fig, = plt.plot(use_freq, mean_wake_sal, color = 'b')
wakecno_fig, = plt.plot(use_freq, mean_wake_cno, color = 'r')
plt.fill_between(use_freq, mean_wake_sal-wakesal_sem_line, mean_wake_sal+wakesal_sem_line,
    alpha=0.2, edgecolor='b', facecolor='b')
plt.fill_between(use_freq, mean_wake_cno-wakecno_sem_line, mean_wake_cno+wakecno_sem_line,
    alpha=0.2, edgecolor='r', facecolor='r')

ax1.spines['top'].set_color('none')
ax1.spines['right'].set_color('none')
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
plt.title((stage1), fontsize = 12, x = 0.5, fontweight = 'demi')
ax1.set_xlim([0,50])
#ax1.set_ylim([0,300])

ax2 = fig.add_subplot(312)
swssal_fig, = plt.plot(use_freq, mean_sws_sal, color = 'b')
swscno_fig, = plt.plot(use_freq, mean_sws_cno, color = 'r')
plt.fill_between(use_freq, mean_sws_sal-swssal_sem_line, mean_sws_sal+swssal_sem_line,
    alpha=0.2, edgecolor='b', facecolor='b')
plt.fill_between(use_freq, mean_sws_cno-swscno_sem_line, mean_sws_cno+swscno_sem_line,
    alpha=0.2, edgecolor='r', facecolor='r')


ax2.spines['top'].set_color('none')
ax2.spines['right'].set_color('none')
ax2.xaxis.set_ticks_position('bottom')
ax2.yaxis.set_ticks_position('left')
plt.title(stage2, fontsize = 12, x = 0.5, y = 0.8, fontweight = 'demi')
ax2.set_xlim([0, 50])
#ax2.set_ylim([0,300])


ax3 = fig.add_subplot(313)
remsal_fig, = plt.plot(use_freq, mean_rem_sal, color = 'b')
remcno_fig = plt.plot(use_freq, mean_rem_cno, color = 'r')
plt.fill_between(use_freq, mean_rem_sal-remsal_sem_line, mean_rem_sal+remsal_sem_line,
    alpha=0.2, edgecolor='b', facecolor='b')
plt.fill_between(use_freq, mean_rem_cno-remcno_sem_line, mean_rem_cno+remcno_sem_line,
    alpha=0.2, edgecolor='r', facecolor='r')

ax3.spines['top'].set_color('none')
ax3.spines['right'].set_color('none')
ax3.xaxis.set_ticks_position('bottom')
ax3.yaxis.set_ticks_position('left')
plt.title(stage3, fontsize = 12, x = 0.5, y = 0.8, fontweight = 'demi')
ax3.set_xlim([0, 50])
#ax3.set_ylim([0,200])


# Turn off axis lines and ticks of the big subplot
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

plt.suptitle(fig_title, fontsize = 15, color = 'b')
plt.figlegend((wakesal_fig, wakecno_fig),(condition1, condition2), loc = 'upper right', fontsize = 10, frameon = False) 
ax.set_xlabel('frequency (Hz)', fontsize = 14)
ax.set_ylabel('spectral power (% epoch total power)', fontsize = 14)

plt.hold(False)

delta_lower = delta[0]
delta_upper = max(mb.find(use_freq < delta[1]))
lotheta_lower = min(mb.find(use_freq > lotheta[0]))
lotheta_upper = max(mb.find(use_freq < lotheta[1]))
hitheta_lower = min(mb.find(use_freq > hitheta[0]))
hitheta_upper = max(mb.find(use_freq < hitheta[1]))
alph_lower = min(mb.find(use_freq > alph[0]))
alph_upper = max(mb.find(use_freq < alph[1]))
gamma_lower = min(mb.find(use_freq > gamma[0]))    
gamma_upper = max(mb.find(use_freq < gamma[1]))
higamma_lower = max(mb.find(use_freq < higamma[0]))    
higamma_upper = max(mb.find(use_freq < higamma[1]))
    
wakesal_delta = masked_wake_sal_values[:,delta_lower:delta_upper]
wakesal_lotheta = masked_wake_sal_values[:,lotheta_lower:lotheta_upper]
wakesal_hitheta = masked_wake_sal_values[:,hitheta_lower:hitheta_upper]
wakesal_alph = masked_wake_sal_values[:,alph_lower:alph_upper]
wakesal_gamma = masked_wake_sal_values[:,gamma_lower:gamma_upper]
wakesal_higamma = masked_wake_sal_values[:,higamma_lower:higamma_upper]

wakesal_mean_delta = np.mean(wakesal_delta, axis = 1)
wakesal_mean_lotheta = np.mean(wakesal_lotheta, axis = 1)
wakesal_mean_hitheta = np.mean(wakesal_hitheta, axis = 1)
wakesal_mean_alph = np.mean(wakesal_alph, axis = 1)
wakesal_mean_gamma = np.mean(wakesal_gamma, axis = 1)
wakesal_mean_higamma = np.mean(wakesal_higamma, axis = 1)

wakesal_mean_delta = list(wakesal_mean_delta)
wakesal_mean_lotheta = list(wakesal_mean_lotheta)
wakesal_mean_hitheta = list(wakesal_mean_hitheta)
wakesal_mean_alph = list(wakesal_mean_alph)
wakesal_mean_gamma = list(wakesal_mean_gamma)
wakesal_mean_higamma = list(wakesal_mean_higamma)

wakecno_delta = masked_wake_cno_values[:,delta_lower:delta_upper]
wakecno_lotheta = masked_wake_cno_values[:,lotheta_lower:lotheta_upper]
wakecno_hitheta = masked_wake_cno_values[:,hitheta_lower:hitheta_upper]
wakecno_alph = masked_wake_cno_values[:,alph_lower:alph_upper]
wakecno_gamma = masked_wake_cno_values[:,gamma_lower:gamma_upper]
wakecno_higamma = masked_wake_cno_values[:,higamma_lower:higamma_upper]

wakecno_mean_delta = np.mean(wakecno_delta, axis = 1)
wakecno_mean_lotheta = np.mean(wakecno_lotheta, axis = 1)
wakecno_mean_hitheta = np.mean(wakecno_hitheta, axis = 1)
wakecno_mean_alph = np.mean(wakecno_alph, axis = 1)
wakecno_mean_gamma = np.mean(wakecno_gamma, axis = 1)
wakecno_mean_higamma = np.mean(wakecno_higamma, axis = 1)

swssal_delta = masked_sws_sal_values[:,delta_lower:delta_upper]
swssal_lotheta = masked_sws_sal_values[:,lotheta_lower:lotheta_upper]
swssal_hitheta = masked_sws_sal_values[:,hitheta_lower:hitheta_upper]
swssal_alph = masked_sws_sal_values[:,alph_lower:alph_upper]
swssal_gamma = masked_sws_sal_values[:,gamma_lower:gamma_upper]
swssal_higamma = masked_sws_sal_values[:,higamma_lower:higamma_upper]


swssal_mean_delta = np.mean(swssal_delta, axis = 1)
swssal_mean_lotheta = np.mean(swssal_lotheta, axis = 1)
swssal_mean_hitheta = np.mean(swssal_hitheta, axis = 1)
swssal_mean_alph = np.mean(swssal_alph, axis = 1)
swssal_mean_gamma = np.mean(swssal_gamma, axis = 1)
swssal_mean_higamma = np.mean(swssal_higamma, axis = 1)

swscno_delta = masked_sws_cno_values[:,delta_lower:delta_upper]
swscno_lotheta = masked_sws_cno_values[:,lotheta_lower:lotheta_upper]
swscno_hitheta = masked_sws_cno_values[:,hitheta_lower:hitheta_upper]
swscno_alph = masked_sws_cno_values[:,alph_lower:alph_upper]
swscno_gamma = masked_sws_cno_values[:,gamma_lower:gamma_upper]
swscno_higamma = masked_sws_cno_values[:,higamma_lower:higamma_upper]

swscno_mean_delta = np.mean(swscno_delta, axis = 1)
swscno_mean_lotheta = np.mean(swscno_lotheta, axis = 1)
swscno_mean_hitheta = np.mean(swscno_hitheta, axis = 1)
swscno_mean_alph = np.mean(swscno_alph, axis = 1)
swscno_mean_gamma = np.mean(swscno_gamma, axis = 1)
swscno_mean_higamma = np.mean(swscno_higamma, axis = 1)

remsal_delta = masked_rem_sal_values[:,delta_lower:delta_upper]
remsal_lotheta = masked_rem_sal_values[:,lotheta_lower:lotheta_upper]
remsal_hitheta = masked_rem_sal_values[:,hitheta_lower:hitheta_upper]
remsal_alph = masked_rem_sal_values[:,alph_lower:alph_upper]
remsal_gamma = masked_rem_sal_values[:,gamma_lower:gamma_upper]
remsal_higamma = masked_rem_sal_values[:,higamma_lower:higamma_upper]

remsal_mean_delta = np.mean(remsal_delta, axis = 1)
remsal_mean_lotheta = np.mean(remsal_lotheta, axis = 1)
remsal_mean_hitheta = np.mean(remsal_hitheta, axis = 1)
remsal_mean_alph = np.mean(remsal_alph, axis = 1)
remsal_mean_gamma = np.mean(remsal_gamma, axis = 1)
remsal_mean_higamma = np.mean(remsal_higamma, axis = 1)

remcno_delta = masked_rem_cno_values[:,delta_lower:delta_upper]
remcno_lotheta = masked_rem_cno_values[:,lotheta_lower:lotheta_upper]
remcno_hitheta = masked_rem_cno_values[:,hitheta_lower:hitheta_upper]
remcno_alph = masked_rem_cno_values[:,alph_lower:alph_upper]
remcno_gamma = masked_rem_cno_values[:,gamma_lower:gamma_upper]
remcno_higamma = masked_rem_cno_values[:,higamma_lower:higamma_upper]

remcno_mean_delta = np.mean(remcno_delta, axis = 1)
remcno_mean_lotheta = np.mean(remcno_lotheta, axis = 1)
remcno_mean_hitheta = np.mean(remcno_hitheta, axis = 1)
remcno_mean_alph = np.mean(remcno_alph, axis = 1)
remcno_mean_gamma = np.mean(remcno_gamma, axis = 1)
remcno_mean_higamma = np.mean(remcno_higamma, axis = 1)

wakesal_bands = np.array([[wakesal_mean_delta], [wakesal_mean_lotheta], [wakesal_mean_hitheta], [wakesal_mean_alph],[wakesal_mean_gamma],[wakesal_mean_higamma]])
wakesal_bands = np.reshape(wakesal_bands,(np.size(wakesal_bands, axis = 0), np.size(wakesal_bands, axis = 2)))
wakecno_bands = np.array([[wakecno_mean_delta], [wakecno_mean_lotheta], [wakecno_mean_hitheta], [wakecno_mean_alph],[wakecno_mean_gamma],[wakecno_mean_higamma]])
wakecno_bands = np.reshape(wakecno_bands,(np.size(wakecno_bands, axis = 0), np.size(wakecno_bands, axis = 2)))

swssal_bands = np.array([[swssal_mean_delta], [swssal_mean_lotheta], [swssal_mean_hitheta], [swssal_mean_alph],[swssal_mean_gamma],[swssal_mean_higamma]])
swssal_bands = np.reshape(swssal_bands,(np.size(swssal_bands, axis = 0), np.size(swssal_bands, axis = 2)))
swscno_bands = np.array([[swscno_mean_delta], [swscno_mean_lotheta], [swscno_mean_hitheta], [swscno_mean_alph],[swscno_mean_gamma],[swscno_mean_higamma]])
swscno_bands = np.reshape(swscno_bands,(np.size(swscno_bands, axis = 0), np.size(swscno_bands, axis = 2)))

remsal_bands = np.array([[remsal_mean_delta], [remsal_mean_lotheta], [remsal_mean_hitheta], [remsal_mean_alph],[remsal_mean_gamma],[remsal_mean_higamma]])
remsal_bands = np.reshape(remsal_bands,(np.size(remsal_bands, axis = 0), np.size(remsal_bands, axis = 2)))
remcno_bands = np.array([[remcno_mean_delta], [remcno_mean_lotheta], [remcno_mean_hitheta], [remcno_mean_alph],[remcno_mean_gamma],[remcno_mean_higamma]])
remcno_bands = np.reshape(remcno_bands,(np.size(remcno_bands, axis = 0), np.size(remcno_bands, axis = 2)))

wakesal_means = np.mean(wakesal_bands, axis = 1)
wakecno_means = np.mean(wakecno_bands, axis = 1)
swssal_means = np.mean(swssal_bands, axis = 1)
swscno_means = np.mean(swscno_bands, axis = 1)
remsal_means = np.mean(remsal_bands, axis = 1)
remcno_means = np.mean(remcno_bands, axis = 1)

wakesal_sem = st.sem(wakesal_bands, axis = 1)
wakecno_sem = st.sem(wakecno_bands, axis = 1)
swssal_sem = st.sem(swssal_bands, axis = 1)
swscno_sem = st.sem(swscno_bands, axis = 1)
remsal_sem = st.sem(remsal_bands, axis = 1)
remcno_sem = st.sem(remcno_bands, axis = 1)

index = np.arange(np.size(wakesal_means))
bar_width = 0.35


fig2 = plt.figure(facecolor = 'w')
bax = fig2.add_subplot(111)

plt.hold(True)

bax1 = fig2.add_subplot(311)
wakesal_rects = plt.bar(index, wakesal_means, bar_width, color ='b', yerr = [np.zeros(np.size(wakesal_sem)),wakesal_sem], ecolor = 'b', label = condition1)
wakecno_rects = plt.bar(index + bar_width, wakecno_means, bar_width, color ='r', yerr = [np.zeros(np.size(wakecno_sem)),wakecno_sem], ecolor = 'r', label = condition2)

bax1.spines['top'].set_color('none')
bax1.spines['right'].set_color('none')
bax1.xaxis.set_ticks_position('bottom')
bax1.yaxis.set_ticks_position('none')
bax1.set_xticklabels([])
bax1.set_title((stage1), fontsize = 12, x = 0.5, fontweight = 'demi')


bax2 = fig2.add_subplot(312)
swssal_rects = plt.bar(index, swssal_means, bar_width, color ='b', yerr = [np.zeros(np.size(swssal_sem)),swssal_sem], ecolor = 'b', label = condition1)
swscno_rects = plt.bar(index + bar_width, swscno_means, bar_width, color ='r', yerr = [np.zeros(np.size(swscno_sem)),swscno_sem], ecolor = 'r', label = condition2)

bax2.spines['top'].set_color('none')
bax2.spines['right'].set_color('none')
bax2.xaxis.set_ticks_position('bottom')
bax2.yaxis.set_ticks_position('none')
bax2.set_xticklabels([])
plt.title((stage2), fontsize = 12, x = 0.5, fontweight = 'demi')


bax3 = fig2.add_subplot(313)
remsal_rects = plt.bar(index, remsal_means, bar_width, color ='b', yerr = [np.zeros(np.size(remsal_sem)),remsal_sem], ecolor = 'b', label = condition1)
remcno_rects = plt.bar(index + bar_width, remcno_means, bar_width, color ='r', yerr = [np.zeros(np.size(remcno_sem)),remcno_sem], ecolor = 'r', label = condition2)

bax3.spines['top'].set_color('none')
bax3.spines['right'].set_color('none')
bax3.xaxis.set_ticks_position('bottom')
bax3.yaxis.set_ticks_position('none')
plt.title((stage3), fontsize = 12, x = 0.5, fontweight = 'demi')
plt.xticks(index + bar_width, (str(delta[0]) + '-' + str(delta[1]), str(lotheta[0]) + '-' + str(lotheta[1]), str(hitheta[0]) + '-' + str(hitheta[1]), str(alph[0]) + '-' + str(alph[1]), str(gamma[0]) + '-' + str(gamma[1])))

bax.spines['top'].set_color('none')
bax.spines['bottom'].set_color('none')
bax.spines['left'].set_color('none')
bax.spines['right'].set_color('none')
bax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

plt.suptitle('Spectral power comparison between band widths \n' + fig_title, fontsize = 15, color = 'b')
plt.subplots_adjust(top=0.85)
bax.set_ylabel('mean spectral power (% epoch)', fontsize = 14)

plt.figlegend((wakesal_rects, wakecno_rects),(condition1, condition2),loc = 'upper right', fontsize = 10, frameon = False) 

plt.hold(False)