# BinauralModeling
Functions to model binaural sensitivity based on the output of a phenomenological auditory nerve model

As a front end, we use the UR_EAR phenomenological auditory nerve (AN) model available at https://www.urmc.rochester.edu/labs/carney.aspx

We have modified the low-pass filter of the inner hair cell in order to limit phase locking, thus providing a better match to binaural sensitivity data in human listeners. Specifically, in the "model_IHC_BEZ2018.c" source code (line 358), the low-pass cutoff was reduced to 1250 Hz (from 3000 Hz), and the mex functions were recompiled.

The two functions "model_IHC_BEZ2018.m" and "model_Synapse_BEZ2018.m" are included in this repository, but please cite the relevant papers provided at the link above if this code is to be used, in particular, Bruce et al. (2018) https://doi.org/10.1016/j.heares.2017.12.016

Also included is a stimulus file ('BUG_T16_1_3.wav') from the BUG corpus described in Kidd et al. (2008) https://asa.scitation.org/doi/10.1121/1.2998980. Please cite if this stimulus is to be used.

Also included is a centrality weighting function written by Arturo Moncada-Torres. This function contains a number of centrailty weight options that can be applied to the SCC. In our implementation, we use the centrailty weight described by Stern & Shear (1996) https://doi.org/10.1121/1.417937
