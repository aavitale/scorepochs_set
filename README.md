# scorepochs_set

This repository is based on the "scorEpochs" toolbox created by Matteo Fraschini and Matteo Demuru
(paper ref)

The purpose of this tool is to provide an automatic scoring of (resting state) MEEG epochs 
based on the power spectrum density 

The main improvements in this package are:
1) INPUT:  eeg_structure in the eeglab format (.set) with also the channel information
2) OUTPUT: score value not only averaged across channels but also for each channel (averaged across epochs)
          ![GitHub Logo](/images/logo.png)
3)         multichannel scroll and topoplot (for a single epoch) of the outlier channels 
           in the low or high frequencies

