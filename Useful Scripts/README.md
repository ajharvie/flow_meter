Some MATLAB scripts used for producing plots, and simulating the flow meter. 
1) "DFT_Signals_Fig4.m", was used to produce figure 4 in the main paper. 
2) "FFT_vs_SDFT.m" produces figure S1. It requires "sdft_plug.m" to work. 
3) "Signal_Generator.m" was used to produce the simulated signals from the supporting information.
4) "Simulator.m" does what it says on the tin - use the signals generated with "Signal_Generator.m". The first 1500 output samples (consistent with the first 15s, or the first full sDFT window) are not output, as these are related to the initial "filling" of the sDFT.
