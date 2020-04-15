# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 13:02:06 2020

@author: annevenner
"""

import numpy as np

def print_error_log(fname, mouse, stage): 
    freq_start = np.loadtxt(fname,skiprows = 1)[0,0]
    freq_end = np.loadtxt(fname,skiprows = 1)[-1,0]   
    
    with open(fname) as f:
        first_line = f.readline()                           
    time_processed = float(first_line.split('Time processed: ')[1].split(' s')[0])
    fft_type = first_line.split('FFT: ')[1].split('"')[0]
    mouseName = first_line.split(' ')[1].split('_')[0].lower()
    fileStage = first_line.split(': ')[4].split('. ')[1].lower()
    startTime = float(first_line.split(': ')[1].split(' ')[0])
    endTime = float(first_line.split(': ')[1].split(' ')[2])
    timeWindow = endTime-startTime
    
    if time_processed < 300:
        print(fname + ': only ' + str(time_processed) + ' secs of data')
    if fft_type != '1024 Hanning':
        print(fname + ': wrong FFT type (not 1024 Hanning')
    if mouseName != mouse:
        print(fname + ': mouse name in file is ' + mouseName)
    if fileStage != stage[1:]:
        print(fname + ' does not match stage in file: ' + fileStage)
    if (timeWindow > 7201.0) or (timeWindow < 2700.0):
        print(fname + ': timeWindow out of range, ' + str(timeWindow))
    if (freq_start != 0.0) or (freq_end < 249):
        print(fname + ': spectra out of range, ' + str(freq_start) + ' to ' + str(freq_end) + ' Hz' )