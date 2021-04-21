# BinauralModeling
Functions to model binaural sensitivity based on the output of a phenomenological auditory nerve model

As a front end, we use the UR_EAR phenomenological auditory nerve (AN) model available at https://www.urmc.rochester.edu/labs/carney.aspx

We have modified the low-pass filter of the inner hair cell in order to limit phase locking, thus providing a better match to binaural sensitivity data in human listeners. Specifically, in the "model_IHC_BEZ2018.c" source code (line 358), the low-pass cutoff was reduced to 1250 Hz (from 3000 Hz), and the .mat function was recompiled.

The two functions "model_IHC_BEZ2018.m" and "model_Synapse_BEZ2018.m" are included in this repository, but please cite the relevant papers provided at the link above if this code is to be used.
