# scorepochs_set

This repository is based on the "scorEpochs" toolbox created by Matteo Fraschini and Matteo Demuru
(paper: https://www.biorxiv.org/content/10.1101/2020.05.26.116434v1

The purpose of this tool is to provide an automatic scoring of (resting state) MEEG epochs 
based on the power spectrum density 

The main improvements in the package are:
1) INPUT:  eeg_structure in the eeglab format (**.set**) with channels location 
2) OUTPUT:
 
2a) score value not only averaged across channels but also for each channel (averaged across epochs)
          ![scorepochs](https://github.com/aavitale/scorepochs_set/blob/main/fig1_scorepoch_xchannel.jpg)
2b) **multichannel scroll and topoplot** (for a single epoch) 
   of the outlier channels in the low and high frequencies
   ![multichanscroll](https://github.com/aavitale/scorepochs_set/blob/main/fig2_multichannel_scroll_topoplot.jpg)
   
-main function: scorEpochs_set.m
-folders required in the path: plot_multichan_master; freezeColors
