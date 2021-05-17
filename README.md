# BinauralModeling
Functions to model binaural sensitivity based on the output of a phenomenological auditory nerve model

As a front end, we use the UR_EAR phenomenological auditory nerve (AN) model available at https://www.urmc.rochester.edu/labs/carney.aspx

We have modified the low-pass filter of the inner hair cell in order to limit phase locking, thus providing a better match to binaural sensitivity data in human listeners. Specifically, in the "model_IHC_BEZ2018.c" source code (line 358), the low-pass cutoff was reduced to 1250 Hz (from 3000 Hz), and the mex functions were recompiled.

The UR_EAR functions "model_IHC_BEZ2018.m" and "model_Synapse_BEZ2018.m" along with their source code are included in this repository, but please cite the relevant papers provided at the link above if this code is to be used, in particular, Bruce et al. (2018) https://doi.org/10.1016/j.heares.2017.12.016

Also included is a centrality weighting function written by Arturo Moncada-Torres. This function contains a number of centrailty weight options that can be applied to the SCC. In our implementation, we use the centrailty weight described by Stern & Shear (1996) https://doi.org/10.1121/1.417937

Also included is a function "LSOmodelCOC.m", written by Go Ashida and modified by Jonas Klug, distributed using the Apache License. This LSO model is described by Ashida et al. (2017) https://doi.org/10.1371/journal.pcbi.1005903, and was also used by Klug et al. (2020) https://doi.org/10.1121/10.0001602

Also included is a stimulus file ('BUG_T16_1_3.wav') from the BUG corpus described in Kidd et al. (2008) https://asa.scitation.org/doi/10.1121/1.2998980. Please cite if this stimulus is to be used.

The functions "genBEZpsth.m" and"genBEZpsth_stochasic.m" generate spike train outputs from the BEZ auditory nerve model with modified IHC filter

# modeling ITD sensitivity
The function "simITDthresholds.m" simulates ITD thresholds for a given stimulus type (e.g. pure tone, narrowband noise), a given model type ('MSO' vs 'LSO'), and a given set of stimulus parameters. ITD sensitivity is simulated using a pipeline modeled after Moncada-Torres et al. (2018).

Also included as standalone scripts:
The script "simITDthresholds_NBnoise.m" simulates ITD sensitivity for a pair of narrowband noise stimuli (500 Hz and 4kHz center frequencies). These stimuli have been used in a number of binaural studies, (e.g. Spencer et al., 2016), and parameter values were chosen to match this literature.
The script "simITDthresholds_RustleNoise.m" simulates ITD sensitivity for a set of octave-band rustle noise stimuli. Stimulus parameters were chosen to match a dataset collected in our lab that has yet to be published, but an interim summary was presented at the 2020 Binaural Bash.

# modeling temporal weighting functions
the function "simANwghts.m" generates a simulated TWF for the speech token 'BUG_T16_1_3.wav'. See Baltzell et al. (under review) for details

# modeling hearing impairment
The function "simHCloss.m" returns cohc and cihc values corresponding to a particular audiometric hearing loss
This function calls the function "getCHCfun.m", which returns a surface of audiometric thresholds corresponding to a pair of ohc/ihc values, for a particular audimetric frequency
