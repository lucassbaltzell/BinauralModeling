function [dB_detect] = getDetectThresh(an_cf,stim_cf,ANpar,OHC,IHC)
%get detection threshold for an AN cf (an_cf) at a stimulus frequency
%(stim_cf). Detection threshold is defined as a firing rate greater than 10
%spikes/s above spontaneous rate
%INPUTS
%an_cf = center frequency of auditory nerve fibers
%stim_cf = stimulus frequency
%ANpar = parameters for auditory nerve
%OHC = cohc parameter in AN model
%IHC = cihc parameter in AN model
%OUPUTS
%dB_detect = detection threshold in dB SPL

%created by Luke Baltzell 05/13/21

if nargin == 2 || isempty(ANpar) == 1
    ANpar.nrep = 1;       %number of repetitions for the psth
    ANpar.cohc = 1;       %OHC scaling factor: 1 is normal OHC function; 0 is complete OHC dysfunction
    ANpar.cihc = 1;       %IHC scaling factor: 1 is normal IHC function; 0 is complete IHC dysfunction
    ANpar.species = 2;    %human with BM tuning from Shera et al. (PNAS 2002)
    ANpar.noiseType = 1;  %1 for variable fGn and 0 for fixed (frozen) fGn
    ANpar.implnt = 0;     %implnt is for "approxiate" or "actual" implementation of the power-law functions: "0" for approx. and "1" for actual implementation
    ANpar.spont.mu = [0.1 4 70]; %[low_sr med_sr high_sr] row vector
    ANpar.spont.sd = [0.1 4 30]; %[low_sr med_sr high_sr] row vector
    ANpar.spont.lims = [1e-3 0.2; 0.2 18; 18 180]; % [low_sr; med_sr; high_sr]
    ANpar.tabs_lims = [0.0002085 0.0006915];  %tabs is the range of absolute refractory period in s
    ANpar.trel_lims = [0.000131 0.000894];  %trel is the range of relative refractory period in s
end

if nargin == 5
    ANpar.cohc = OHC;
    ANpar.cihc = IHC;
end

%% define probe stimuli (Zhang et al., 2001; figure 3)
fs = 100000;
dur = 0.05;
t = [1/fs:1/fs:dur];
tc = 0.0025;
stim = sin(2*pi*stim_cf*t);
stim = rampdamp(stim,tc,fs);
stim0 = 0*stim;

%% get detection threshold
Nfibers = 200;
dB = [0:5:90];
flg = 0;
ind = 0;
while flg == 0
    ind = ind+1;
    [psth,x] = genBEZpsth_stochastic(stim,fs,dB(ind),an_cf,Nfibers,ANpar);
    rate(ind) = (sum(sum(squeeze(psth)))/Nfibers)/x(end);
    psth0 = genBEZpsth_stochastic(stim0,fs,dB(ind),an_cf,Nfibers,ANpar);
    rate0(ind) = (sum(sum(squeeze(psth0)))/Nfibers)/x(end);
    if rate(ind) - rate0(ind) > 10
        flg = 1;
    elseif ind == length(dB)
        flg = 1;
    end
end
dB_detect = dB(ind);

end
