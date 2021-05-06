# BinauralModeling
Functions to model binaural sensitivity based on the output of a phenomenological auditory nerve model

As a front end, we use the UR_EAR phenomenological auditory nerve (AN) model available at https://www.urmc.rochester.edu/labs/carney.aspx

We have modified the low-pass filter of the inner hair cell in order to limit phase locking, thus providing a better match to binaural sensitivity data in human listeners. Specifically, in the "model_IHC_BEZ2018.c" source code (line 358), the low-pass cutoff was reduced to 1250 Hz (from 3000 Hz), and the mex functions were recompiled.

The UR_EAR functions "model_IHC_BEZ2018.m" and "model_Synapse_BEZ2018.m" along with their source code are included in this repository, but please cite the relevant papers provided at the link above if this code is to be used, in particular, Bruce et al. (2018) https://doi.org/10.1016/j.heares.2017.12.016

Also included is a stimulus file ('BUG_T16_1_3.wav') from the BUG corpus described in Kidd et al. (2008) https://asa.scitation.org/doi/10.1121/1.2998980. Please cite if this stimulus is to be used.

Also included is a centrality weighting function written by Arturo Moncada-Torres. This function contains a number of centrailty weight options that can be applied to the SCC. In our implementation, we use the centrailty weight described by Stern & Shear (1996) https://doi.org/10.1121/1.417937

The functions "genBEZpsth.m" and"genBEZpsth_stochasic.m" generate spike train outputs from the BEZ auditory nerve model with momdified IHC filter

The script "getITDthresholds_PureTones.m" generates predicted ITD sensitivity to a set of pure tone stimuli. Stimulus parameters were chosen to follow Brughera et al. (2013)
The script "getITDthresholds_NBnoise.m" generates predicted ITD sensitivity to a pair of narrowband noise stimuli (500 Hz and 4kHz center frequencies). These stimuli have been used in a number of binaural studies, (e.g. Spencer et al., 2016), and parameter values were chosen to match this literature.
The script "getITDthresholds_RustleNoise.m" generates predicted ITD sensitivity to a set of octave-band rustle noise stimuli. Stimulus parameters were chosen to match a dataset collected in our lab that has yet to be published, but an interim summary was presented at the 2020 Binaural Bash.

ITD sensitivity is simulated in a similar fashion for these different stimuli, in a pipeline modeled after Moncada-Torres et al. (2018).
