AV 13th February 2020
Requires Python 2.7 to run
Scripts developed to analyse FFT data exported from SleepSign2, for different sleep-wake states for different conditions. FFT records could not be baselined to previous recordings from the same individual so the normalization method is to express each frequency bin as a percent of the summed power for all frequency bin for each epoch.

proc_anne_selfBL.py   Loads a txt file exported from Sleepsign2 (containing the FFT for each epoch evaluated during the recording) using 'loadtxt', applies the normalization function, plots average spectra for each sleep-wake state for that individual and exports the average spectra for each sleep-wake state and the frequency binning into separate text files 

group_spectral_selfBL   Loads average spectra for each sleep-wake state and for each individual and compares and plots 1) average (with SEM) sleep-wake spectra for individuals from the same group (e.g. 'sham' vs 'tbi1') and 2) averaged spectral power for a range of FFT frequency bands
