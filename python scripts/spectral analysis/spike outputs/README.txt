AV 30th March 2020. 
Requires Python 2.7 to run.
These scripts were developed to analyse FFT generated from Spike2 recordings, during different sleep-wake states and in different conditions. These scripts were specifically developed to analysed how spectral activity changes before and after a particular conditions (e.g. a sustained optogenetic stimulation or behavioral intervention)

check_spike_spectra.py    Checks the files exported from the Spike2 Online Sleep Detection Script Gated Power function.

spike_spectra.py    Compares (and plots) FFT spectra from each sleep-wake state for 2 different conditions for each individual. The data is normalized by expressing each frequency bin as a percentage of the same frequency bin from data taken during the pre-experimental time window.
The script then averages the information from all individuals and comapres and plots 1) an average spectra with SEM for each condition and 2) a bar chart detailing averaged spectral power for a range of FFT bands.

spectralFuncs.py:   contains functions used by the other scripts 
