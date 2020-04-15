# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 15:16:51 2020

@author: annevenner
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mb
import numpy.ma as ma
import scipy.stats as st

path = 'E:/ls project/analysis_mar2020/spectral/Tapedown FFT'
fig_folder  = 'E:/ls project/analysis_mar2020/spectral/Tapedown FFT/individual mouse plots/'
fig_title = 'Baseline vs Tapedown Spectra: Crfr2-cre mice'

### USER-DEFINED VARIABLES ####
stage1 = '_wake'
stage2 = '_nrem'
stage3 = '_rem'
condition1 = '_bl'   
condition2 = '_td'
time1 = '_pre'
time2 = '_post'

indx_freq_start = 1
indx_freq_end = -1
hz60_start = 58
hz60_end = 62

value_threshold = 0

delta = [0.5,4]
lotheta = [4,6]
hitheta = [6,10]
alph = [10,20]
gamma = [20,50] 
higamma = [70,120]

########## 

dirList= os.listdir(path)    # insert the path to the directory of interest
filename = []

dirList = [i.lower() for i in dirList];

for fname in dirList:
    if fname.endswith('.txt'):        
        filename = np.append(filename, fname)
        
freq = np.loadtxt(filename[0],skiprows = 1)[:,0]
use_freq = freq[indx_freq_start:indx_freq_end] 

noisy = np.logical_and(freq > hz60_start, freq < hz60_end)
noise = mb.find(noisy == True)


xName = list(filename) #converts  filename into a list (rather than strings)
masterName = []
listName = []
numMice = []
phenotypic = []

for i in np.arange(filename.size):
    masterName = np.append(masterName,xName[i].split('_')[0])
    
masterName = list(set(masterName))
numMice = np.arange(0,len(masterName),1)

all_wake_bl_pre = np.empty((numMice.size, 512))
all_wake_bl_post = np.empty((numMice.size, 512))
all_wake_td_pre = np.empty((numMice.size, 512))
all_wake_td_post = np.empty((numMice.size, 512))

all_nrem_bl_pre = np.empty((numMice.size, 512))
all_nrem_bl_post = np.empty((numMice.size, 512))
all_nrem_td_pre = np.empty((numMice.size, 512))
all_nrem_td_post = np.empty((numMice.size, 512))

all_rem_bl_pre = np.empty((numMice.size, 512))
all_rem_bl_post = np.empty((numMice.size, 512))
all_rem_td_pre = np.empty((numMice.size, 512))
all_rem_td_post = np.empty((numMice.size, 512))
    
for i in np.arange(numMice.size):
    mouse = masterName[i]
    num_wake_bl_pre = 0
    num_wake_bl_post = 0
    num_wake_td_pre = 0
    num_wake_td_post = 0
    num_nrem_bl_pre = 0
    num_nrem_bl_post = 0
    num_nrem_td_pre = 0
    num_nrem_td_post = 0
    num_rem_bl_pre = 0
    num_rem_bl_post = 0
    num_rem_td_pre = 0
    num_rem_td_post = 0          
    for fname in dirList:
        if fname.endswith('.txt'):
            if (mouse in fname) and (stage1 in fname) and (condition1 in fname) and (time1 in fname):
                wake_bl_pre = np.loadtxt(fname, skiprows = 1)[:,1]
                num_wake_bl_pre += 1
            elif (mouse in fname) and (stage1 in fname) and (condition1 in fname) and (time2 in fname):
                wake_bl_post = np.loadtxt(fname, skiprows = 1)[:,1]
                num_wake_bl_post += 1
            elif (mouse in fname) and (stage1 in fname) and (condition2 in fname) and (time1 in fname):
                wake_td_pre = np.loadtxt(fname, skiprows = 1)[:,1]
                num_wake_td_pre += 1
            elif (mouse in fname) and (stage1 in fname) and (condition2 in fname) and (time2 in fname):
                wake_td_post = np.loadtxt(fname, skiprows = 1)[:,1]
                num_wake_td_post += 1
            elif (mouse in fname) and (stage2 in fname) and (condition1 in fname) and (time1 in fname):
                nrem_bl_pre = np.loadtxt(fname, skiprows = 1)[:,1]
                num_nrem_bl_pre += 1
            elif (mouse in fname) and (stage2 in fname) and (condition1 in fname) and (time2 in fname):
                nrem_bl_post = np.loadtxt(fname, skiprows = 1)[:,1]
                num_nrem_bl_post += 1
            elif (mouse in fname) and (stage2 in fname) and (condition2 in fname) and (time1 in fname):
                nrem_td_pre = np.loadtxt(fname, skiprows = 1)[:,1]
                num_nrem_td_pre += 1
            elif (mouse in fname) and (stage2 in fname) and (condition2 in fname) and (time2 in fname):
                nrem_td_post = np.loadtxt(fname, skiprows = 1)[:,1]
                num_nrem_td_post += 1
            elif (mouse in fname) and (stage3 in fname) and (condition1 in fname) and (time1 in fname):
                rem_bl_pre = np.loadtxt(fname, skiprows = 1)[:,1]
                num_rem_bl_pre += 1
            elif (mouse in fname) and (stage3 in fname) and (condition1 in fname) and (time2 in fname):
                rem_bl_post = np.loadtxt(fname, skiprows = 1)[:,1]
                num_rem_bl_post += 1
            elif (mouse in fname) and (stage3 in fname) and (condition2 in fname) and (time1 in fname):
                rem_td_pre = np.loadtxt(fname, skiprows = 1)[:,1]
                num_rem_td_pre += 1
            elif (mouse in fname) and (stage3 in fname) and (condition2 in fname) and (time2 in fname):
                rem_td_post = np.loadtxt(fname, skiprows = 1)[:,1]
                num_rem_td_post+= 1
        
    if  num_wake_bl_pre == 0:
        all_wake_bl_pre[i,:] = -1
    else:
        all_wake_bl_pre[i,:] = wake_bl_pre
        
    if  num_wake_bl_post == 0:
        all_wake_bl_post[i,:] = -1
    else:
        all_wake_bl_post[i,:] = wake_bl_post
        
    if  num_wake_td_pre == 0:    
        all_wake_td_pre[i,:] = -1
    else:
        all_wake_td_pre[i,:] = wake_td_pre
        
    if  num_wake_td_post == 0:  
        all_wake_td_post[i,:] = -1
    else:
        all_wake_td_post[i,:] = wake_td_post
     
    if  num_nrem_bl_pre == 0:
        all_nrem_bl_pre[i,:] = -1
    else: 
        all_nrem_bl_pre[i,:] = nrem_bl_pre
        
    if  num_nrem_bl_post == 0:  
        all_nrem_bl_post[i,:] = -1
    else:
        all_nrem_bl_post[i,:] = nrem_bl_post
        
    if  num_nrem_td_pre == 0: 
        all_nrem_td_pre[i,:] = -1
    else:
        all_nrem_td_pre[i,:] = nrem_td_pre
    
    if  num_nrem_td_post == 0:   
        all_nrem_td_post[i,:] = -1
    else:
        all_nrem_td_post[i,:] = nrem_td_post
   
    if  num_rem_bl_pre == 0:
        all_rem_bl_pre[i,:] = -1
    else: 
        all_rem_bl_pre[i,:] = rem_bl_pre
    
    if  num_rem_bl_post == 0:    
        all_rem_bl_post[i,:] = -1
    else:
        all_rem_bl_post[i,:] = rem_bl_post
        
    if num_rem_td_pre == 0:   
        all_rem_td_pre[i,:] = -1
    else:
        all_rem_td_pre[i,:] = rem_td_pre
    
    if  num_rem_td_post == 0: 
        all_rem_td_post[i,:] = -1
    else:
        all_rem_td_post[i,:] = rem_td_post
        
baseline_wake_bl = all_wake_bl_post/all_wake_bl_pre*100      
baseline_wake_td = all_wake_td_post/all_wake_td_pre*100

baseline_nrem_bl = all_nrem_bl_post/all_nrem_bl_pre*100      
baseline_nrem_td = all_nrem_td_post/all_nrem_td_pre*100

baseline_rem_bl = all_rem_bl_post/all_rem_bl_pre*100      
baseline_rem_td = all_rem_td_post/all_rem_td_pre*100

#for i in np.arange(numMice.size):
#    mouse = masterName[i]
#    plt.figure()
#    plt.suptitle(mouse, fontsize = 15)
#    
#    plt.subplot(3,1,1)
#    plt.plot(freq, baseline_wake_bl[i,:],'b')  
#    plt.plot(freq, baseline_wake_td[i,:],'r')
#    plt.legend(('baseline', 'tapedown'))
#    plt.xlim(0,120)
#    plt.title('wake')
#    
#    plt.subplot(3,1,2)
#    plt.plot(freq, baseline_nrem_bl[i,:],'b')  
#    plt.plot(freq, baseline_nrem_td[i,:],'r')
#    plt.ylabel('spectral power (% baseline)')
#    plt.xlim(0,120)
#    plt.title('nrem')
#    
#    plt.subplot(3,1,3)
#    plt.plot(freq, baseline_rem_bl[i,:],'b')  
#    plt.plot(freq, baseline_rem_td[i,:],'r')
#    plt.xlim(0,120)
#    plt.xlabel('frequency (Hz)')
#    plt.title('rem')
#    
#    plt.savefig(fig_folder+mouse + '.png')

baseline_wake_bl[:,noise[0]:noise[-1]] = -1    
baseline_wake_td[:,noise[0]:noise[-1]] = -1
baseline_nrem_bl[:,noise[0]:noise[-1]] = -1      
baseline_nrem_td[:,noise[0]:noise[-1]] = -1
baseline_rem_bl[:,noise[0]:noise[-1]] = -1      
baseline_rem_td[:,noise[0]:noise[-1]] = -1
    
masked_wake_bl = ma.masked_less(baseline_wake_bl, 0)[:,indx_freq_start:indx_freq_end]   
masked_wake_td = ma.masked_less(baseline_wake_td, 0)[:,indx_freq_start:indx_freq_end]
masked_nrem_bl = ma.masked_less(baseline_nrem_bl, 0)[:,indx_freq_start:indx_freq_end]   
masked_nrem_td = ma.masked_less(baseline_nrem_td, 0)[:,indx_freq_start:indx_freq_end]
masked_rem_bl = ma.masked_less(baseline_rem_bl, 0)[:,indx_freq_start:indx_freq_end]   
masked_rem_td = ma.masked_less(baseline_rem_td, 0)[:,indx_freq_start:indx_freq_end]

mean_wake_bl = ma.mean(masked_wake_bl,axis = 0)
mean_wake_td = ma.mean(masked_wake_td,axis = 0)
mean_nrem_bl = ma.mean(masked_nrem_bl,axis = 0)
mean_nrem_td = ma.mean(masked_nrem_td,axis = 0)
mean_rem_bl = ma.mean(masked_rem_bl,axis = 0)
mean_rem_td = ma.mean(masked_rem_td,axis = 0)

sem_wake_bl = ma.std(masked_wake_bl, axis = 0)/np.sqrt(ma.count(masked_wake_bl,axis = 0))
sem_wake_td = ma.std(masked_wake_td, axis = 0)/np.sqrt(ma.count(masked_wake_bl,axis = 0))
sem_nrem_bl = ma.std(masked_nrem_bl, axis = 0)/np.sqrt(ma.count(masked_wake_bl,axis = 0))
sem_nrem_td = ma.std(masked_nrem_td, axis = 0)/np.sqrt(ma.count(masked_wake_bl,axis = 0))
sem_rem_bl = ma.std(masked_rem_bl, axis = 0)/np.sqrt(ma.count(masked_wake_bl,axis = 0))
sem_rem_td = ma.std(masked_rem_td, axis = 0)/np.sqrt(ma.count(masked_rem_td,axis = 0))


# Plot baseline and tapedown on the same figure
fig = plt.figure(facecolor = 'w')
ax = fig.add_subplot(111)

plt.hold(True)
ax1 = fig.add_subplot(311)
wake_bl_fig, = plt.plot(use_freq, mean_wake_bl, color = 'k')
wake_td_fig, = plt.plot(use_freq, mean_wake_td, color = 'r')
plt.fill_between(use_freq, mean_wake_bl-sem_wake_bl, mean_wake_bl+sem_wake_bl,
    alpha=0.2, edgecolor='k', facecolor='k')
plt.fill_between(use_freq, mean_wake_td-sem_wake_td, mean_wake_td+sem_wake_td,
    alpha=0.2, edgecolor='r', facecolor='r')

ax1.spines['top'].set_color('none')
ax1.spines['right'].set_color('none')
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
plt.title((stage1), fontsize = 12, x = 0.5, fontweight = 'demi')
ax1.set_xlim([0,20])
ax1.set_ylim([50,175])

ax2 = fig.add_subplot(312)
nrem_bl_fig, = plt.plot(use_freq, mean_nrem_bl, color = 'k')
nrem_td_fig, = plt.plot(use_freq, mean_nrem_td, color = 'r')
plt.fill_between(use_freq, mean_nrem_bl-sem_nrem_bl, mean_nrem_bl+sem_nrem_bl,
    alpha=0.2, edgecolor='k', facecolor='k')
plt.fill_between(use_freq, mean_nrem_td-sem_nrem_td, mean_nrem_td+sem_nrem_td,
    alpha=0.2, edgecolor='r', facecolor='r')

ax2.spines['top'].set_color('none')
ax2.spines['right'].set_color('none')
ax2.xaxis.set_ticks_position('bottom')
ax2.yaxis.set_ticks_position('left')
plt.title(stage2, fontsize = 12, x = 0.5, y = 0.8, fontweight = 'demi')
ax2.set_xlim([0, 20])
ax2.set_ylim([50,150])

ax3 = fig.add_subplot(313)
rem_bl_fig, = plt.plot(use_freq, mean_rem_bl, color = 'k')
rem_td_fig = plt.plot(use_freq, mean_rem_td, color = 'r')
plt.fill_between(use_freq, mean_rem_bl-sem_rem_bl, mean_rem_bl+sem_rem_bl,
    alpha=0.2, edgecolor='k', facecolor='k')
plt.fill_between(use_freq, mean_rem_td-sem_rem_td, mean_rem_td+sem_rem_td,
    alpha=0.2, edgecolor='r', facecolor='r')
    
# Turn off axis lines and ticks of the big subplot
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['right'].set_color('none')
ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

plt.suptitle(fig_title, fontsize = 15, color = 'b')
plt.figlegend((wake_bl_fig, wake_td_fig),(condition1, condition2), loc = 'upper right', fontsize = 10, frameon = False) 
ax.set_xlabel('frequency (Hz)', fontsize = 14)
ax.set_ylabel('spectral power as a percentage of baseline (%)', fontsize = 14)

plt.hold(False)

ax3.spines['top'].set_color('none')
ax3.spines['right'].set_color('none')
ax3.xaxis.set_ticks_position('bottom')
ax3.yaxis.set_ticks_position('left')
plt.title(stage3, fontsize = 12, x = 0.5, y = 0.8, fontweight = 'demi')
ax3.set_xlim([0, 120])
ax3.set_ylim([50,150])

delta_lower = 0
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

mean_wake_bl_delta = ma.mean(masked_wake_bl[:,delta_lower:delta_upper],axis = 1)
mean_wake_bl_lotheta = ma.mean(masked_wake_bl[:,lotheta_lower:lotheta_upper],axis = 1)
mean_wake_bl_hitheta = ma.mean(masked_wake_bl[:,hitheta_lower:hitheta_upper],axis = 1)
mean_wake_bl_alph = ma.mean(masked_wake_bl[:,alph_lower:alph_upper],axis = 1)
mean_wake_bl_gamma = ma.mean(masked_wake_bl[:,gamma_lower:gamma_upper],axis = 1)
mean_wake_bl_higamma = ma.mean(masked_wake_bl[:,higamma_lower:higamma_upper],axis = 1)

mean_wake_td_delta = ma.mean(masked_wake_td[:,delta_lower:delta_upper],axis = 1)
mean_wake_td_lotheta = ma.mean(masked_wake_td[:,lotheta_lower:lotheta_upper],axis = 1)
mean_wake_td_hitheta = ma.mean(masked_wake_td[:,hitheta_lower:hitheta_upper],axis = 1)
mean_wake_td_alph = ma.mean(masked_wake_td[:,alph_lower:alph_upper],axis = 1)
mean_wake_td_gamma = ma.mean(masked_wake_td[:,gamma_lower:gamma_upper],axis = 1)
mean_wake_td_higamma = ma.mean(masked_wake_td[:,higamma_lower:higamma_upper],axis = 1)

mean_nrem_bl_delta = ma.mean(masked_nrem_bl[:,delta_lower:delta_upper],axis = 1)
mean_nrem_bl_lotheta = ma.mean(masked_nrem_bl[:,lotheta_lower:lotheta_upper],axis = 1)
mean_nrem_bl_hitheta = ma.mean(masked_nrem_bl[:,hitheta_lower:hitheta_upper],axis = 1)
mean_nrem_bl_alph = ma.mean(masked_nrem_bl[:,alph_lower:alph_upper],axis = 1)
mean_nrem_bl_gamma = ma.mean(masked_nrem_bl[:,gamma_lower:gamma_upper],axis = 1)
mean_nrem_bl_higamma = ma.mean(masked_nrem_bl[:,higamma_lower:higamma_upper],axis = 1)

mean_nrem_td_delta = ma.mean(masked_nrem_td[:,delta_lower:delta_upper],axis = 1)
mean_nrem_td_lotheta = ma.mean(masked_nrem_td[:,lotheta_lower:lotheta_upper],axis = 1)
mean_nrem_td_hitheta = ma.mean(masked_nrem_td[:,hitheta_lower:hitheta_upper],axis = 1)
mean_nrem_td_alph = ma.mean(masked_nrem_td[:,alph_lower:alph_upper],axis = 1)
mean_nrem_td_gamma = ma.mean(masked_nrem_td[:,gamma_lower:gamma_upper],axis = 1)
mean_nrem_td_higamma = ma.mean(masked_nrem_td[:,higamma_lower:higamma_upper],axis = 1)

mean_rem_bl_delta = ma.mean(masked_rem_bl[:,delta_lower:delta_upper],axis = 1)
mean_rem_bl_lotheta = ma.mean(masked_rem_bl[:,lotheta_lower:lotheta_upper],axis = 1)
mean_rem_bl_hitheta = ma.mean(masked_rem_bl[:,hitheta_lower:hitheta_upper],axis = 1)
mean_rem_bl_alph = ma.mean(masked_rem_bl[:,alph_lower:alph_upper],axis = 1)
mean_rem_bl_gamma = ma.mean(masked_rem_bl[:,gamma_lower:gamma_upper],axis = 1)
mean_rem_bl_higamma = ma.mean(masked_rem_bl[:,higamma_lower:higamma_upper],axis = 1)

mean_rem_td_delta = ma.mean(masked_rem_td[:,delta_lower:delta_upper],axis = 1)
mean_rem_td_lotheta = ma.mean(masked_rem_td[:,lotheta_lower:lotheta_upper],axis = 1)
mean_rem_td_hitheta = ma.mean(masked_rem_td[:,hitheta_lower:hitheta_upper],axis = 1)
mean_rem_td_alph = ma.mean(masked_rem_td[:,alph_lower:alph_upper],axis = 1)
mean_rem_td_gamma = ma.mean(masked_rem_td[:,gamma_lower:gamma_upper],axis = 1)
mean_rem_td_higamma = ma.mean(masked_rem_td[:,higamma_lower:higamma_upper],axis = 1)


wake_bl_bands = np.vstack((mean_wake_bl_delta, mean_wake_bl_lotheta, mean_wake_bl_hitheta, mean_wake_bl_alph,mean_wake_bl_gamma,mean_wake_bl_higamma))
wake_td_bands = np.vstack((mean_wake_td_delta, mean_wake_td_lotheta, mean_wake_td_hitheta, mean_wake_td_alph,mean_wake_td_gamma,mean_wake_td_higamma))
nrem_bl_bands = np.vstack((mean_nrem_bl_delta, mean_nrem_bl_lotheta, mean_nrem_bl_hitheta, mean_nrem_bl_alph,mean_nrem_bl_gamma,mean_nrem_bl_higamma))
nrem_td_bands = np.vstack((mean_nrem_td_delta, mean_nrem_td_lotheta, mean_nrem_td_hitheta, mean_nrem_td_alph,mean_nrem_td_gamma,mean_nrem_td_higamma))
rem_bl_bands = np.vstack((mean_rem_bl_delta, mean_rem_bl_lotheta, mean_rem_bl_hitheta, mean_rem_bl_alph,mean_rem_bl_gamma,mean_rem_bl_higamma))
rem_td_bands = np.vstack((mean_rem_td_delta, mean_rem_td_lotheta, mean_rem_td_hitheta, mean_rem_td_alph,mean_rem_td_gamma,mean_rem_td_higamma))

wake_bl_means = np.mean(wake_bl_bands, axis = 1)
wake_td_means = np.mean(wake_td_bands, axis = 1)
nrem_bl_means = np.mean(nrem_bl_bands, axis = 1)
nrem_td_means = np.mean(nrem_td_bands, axis = 1)
rem_bl_means = np.mean(rem_bl_bands, axis = 1)
rem_td_means = np.mean(rem_td_bands, axis = 1)

wake_bl_sem = st.sem(wake_bl_bands, axis = 1)
wake_td_sem = st.sem(wake_td_bands, axis = 1)
nrem_bl_sem = st.sem(nrem_bl_bands, axis = 1)
nrem_td_sem = st.sem(nrem_td_bands, axis = 1)
rem_bl_sem = st.sem(rem_bl_bands, axis = 1)
rem_td_sem = st.sem(rem_td_bands, axis = 1)

wake_bl_bands100 = wake_bl_bands - 100
wake_td_bands100 = wake_td_bands - 100
nrem_bl_bands100 = nrem_bl_bands - 100
nrem_td_bands100 = nrem_td_bands - 100 
rem_bl_bands100 = rem_bl_bands - 100
rem_td_bands100 = rem_td_bands - 100

index = np.arange(np.size(wake_bl_means))
bar_width = 0.35


fig2 = plt.figure(facecolor = 'w')
bax = fig2.add_subplot(111)

plt.hold(True)

bax1 = fig2.add_subplot(311)
wake_bl_rects = plt.bar(index, wake_bl_means, bar_width, color ='b', yerr = [np.zeros(np.size(wake_bl_sem)),wake_bl_sem], ecolor = 'b', label = condition1)
wake_td_rects = plt.bar(index + bar_width, wake_td_means, bar_width, color ='r', yerr = [np.zeros(np.size(wake_td_sem)),wake_td_sem], ecolor = 'r', label = condition2)

bax1.spines['top'].set_color('none')
bax1.spines['right'].set_color('none')
bax1.xaxis.set_ticks_position('bottom')
bax1.yaxis.set_ticks_position('none')
bax1.set_xticklabels([])
bax1.set_title((stage1), fontsize = 12, x = 0.5, fontweight = 'demi')


bax2 = fig2.add_subplot(312)
nrem_bl_rects = plt.bar(index, nrem_bl_means, bar_width, color ='b', yerr = [np.zeros(np.size(nrem_bl_sem)),nrem_bl_sem], ecolor = 'b', label = condition1)
nrem_td_rects = plt.bar(index + bar_width, nrem_td_means, bar_width, color ='r', yerr = [np.zeros(np.size(nrem_td_sem)),nrem_td_sem], ecolor = 'r', label = condition2)

bax2.spines['top'].set_color('none')
bax2.spines['right'].set_color('none')
bax2.xaxis.set_ticks_position('bottom')
bax2.yaxis.set_ticks_position('none')
bax2.set_xticklabels([])
plt.title((stage2), fontsize = 12, x = 0.5, fontweight = 'demi')


bax3 = fig2.add_subplot(313)
rem_bl_rects = plt.bar(index, rem_bl_means, bar_width, color ='b', yerr = [np.zeros(np.size(rem_bl_sem)),rem_bl_sem], ecolor = 'b', label = condition1)
rem_td_rects = plt.bar(index + bar_width, rem_td_means, bar_width, color ='r', yerr = [np.zeros(np.size(rem_td_sem)),rem_td_sem], ecolor = 'r', label = condition2)

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
bax.set_ylabel('mean spectral power as a percentage of baseline', fontsize = 14)

plt.figlegend((wake_bl_rects, wake_td_rects),(condition1, condition2),loc = 'upper right', fontsize = 10, frameon = False) 

