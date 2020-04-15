# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 15:16:51 2020

@author: annevenner
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import spectralFuncs as sf

path = 'E:/ls project/analysis_mar2020/spectral/Tapedown FFT'

stage1 = '_wake'
stage2 = '_nrem'
stage3 = '_rem'
condition1 = '_bl'   
condition2 = '_td'
time1 = '_pre'
time2 = '_post'

dirList= os.listdir(path)    # insert the path to the directory of interest
filename = []

dirList = [i.lower() for i in dirList];

for fname in dirList:
    if fname.endswith('.txt'):        
        filename = np.append(filename, fname)

xName = list(filename) #converts  filename into a list (rather than strings)
masterName = []
listName = []
numMice = []
phenotypic = []

for i in np.arange(filename.size):
    masterName = np.append(masterName,xName[i].split('_')[0])
    
masterName = list(set(masterName))
numMice = np.arange(0,len(masterName),1)
   
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
                num_wake_bl_pre += 1
                sf.print_error_log(fname, mouse, stage1)
            elif (mouse in fname) and (stage1 in fname) and (condition1 in fname) and (time2 in fname):
                num_wake_bl_post += 1
                sf.print_error_log(fname, mouse, stage1)
            elif (mouse in fname) and (stage1 in fname) and (condition2 in fname) and (time1 in fname):
                num_wake_td_pre += 1
                sf.print_error_log(fname, mouse, stage1)
            elif (mouse in fname) and (stage1 in fname) and (condition2 in fname) and (time2 in fname):
                num_wake_td_post += 1
                sf.print_error_log(fname, mouse, stage1)
            elif (mouse in fname) and (stage2 in fname) and (condition1 in fname) and (time1 in fname):
                num_nrem_bl_pre += 1
                sf.print_error_log(fname, mouse, stage2)
            elif (mouse in fname) and (stage2 in fname) and (condition1 in fname) and (time2 in fname):
                num_nrem_bl_post += 1
                sf.print_error_log(fname, mouse, stage2)
            elif (mouse in fname) and (stage2 in fname) and (condition2 in fname) and (time1 in fname):
                num_nrem_td_pre += 1
                sf.print_error_log(fname, mouse, stage2)
            elif (mouse in fname) and (stage2 in fname) and (condition2 in fname) and (time2 in fname):
                num_nrem_td_post += 1
                sf.print_error_log(fname, mouse, stage2)
            elif (mouse in fname) and (stage3 in fname) and (condition1 in fname) and (time1 in fname):
                num_rem_bl_pre += 1
                sf.print_error_log(fname, mouse, stage3)
            elif (mouse in fname) and (stage3 in fname) and (condition1 in fname) and (time2 in fname):
                num_rem_bl_post += 1
                sf.print_error_log(fname, mouse, stage3)
            elif (mouse in fname) and (stage3 in fname) and (condition2 in fname) and (time1 in fname):
                num_rem_td_pre += 1
                sf.print_error_log(fname, mouse, stage3)
            elif (mouse in fname) and (stage3 in fname) and (condition2 in fname) and (time2 in fname):
                num_rem_td_post += 1
                sf.print_error_log(fname, mouse, stage3)
                
    if  num_wake_bl_pre == 0:
        print(mouse+stage1+condition1+time1 + ' does not exist')
    if  num_wake_bl_post == 0:
        print(mouse+stage1+condition1+time2 + ' does not exist')
    if  num_wake_td_pre == 0:
        print(mouse+stage1+condition2+time1 + ' does not exist')
    if  num_wake_td_post == 0:
        print(mouse+stage1+condition2+time2 + ' does not exist')
    if  num_nrem_bl_pre == 0:
        print(mouse+stage2+condition1+time1 + ' does not exist')
    if  num_nrem_bl_post == 0:
        print(mouse+stage2+condition1+time2 + ' does not exist')
    if  num_nrem_td_pre == 0:
        print(mouse+stage2+condition2+time1 + ' does not exist')
    if  num_nrem_td_post == 0:
        print(mouse+stage2+condition2+time2 + ' does not exist')
    if  num_rem_bl_pre == 0:
        print(mouse+stage3+condition1+time1 + ' does not exist')
    if  num_rem_bl_post == 0:
        print(mouse+stage3+condition1+time2 + ' does not exist')
    if  num_rem_td_pre == 0:
        print(mouse+stage3+condition2+time1 + ' does not exist')
    if  num_rem_td_post == 0:
        print(mouse+stage3+condition2+time2 + ' does not exist')
      

#pre_data = np.loadtxt('ls28_td_pre_nrem.txt',skiprows = 1)
#frequency = pre_data[:,0];
#pre = pre_data[:,1];
#
#with open('ls28_td_pre_nrem.txt') as f:
#    first_line = f.readline()
#time_processed_pre = float(first_line.split('Time processed: ')[1].split(' s')[0])
#
#post_data = np.loadtxt('ls28_td_post_nrem.txt',skiprows = 1);
#post = post_data[:,1];
#
#with open('ls28_td_post_nrem.txt') as f:
#    first_line = f.readline()
#time_processed_post = float(first_line.split('Time processed: ')[1].split(' s')[0])  
#
#normalized = post/pre*100;
#
#plt.plot(frequency, normalized)
