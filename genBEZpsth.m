function [psth,x] = genBEZpsth(stim,stim_fs,dB,CFs,nfibers,ANpar)
%generates post-stimulus spike histograms for a vector of center
%frequencies and for a set of fibers corresponding to each center
%frequency. In this function we call the Bruce, Erfani & Zilany model 
%described in Bruce et al. (2018): https://doi.org/10.1016/j.heares.2017.12.016
%AN parameters are fixed and given in ANpar

%stim: stimulus, could be mono or stereo
%stim_fs: sampling rate of stimulus
%dB: desired intensity of stimulus in dB. For stereo stimulus, this can
%have two values (left and right)
%CFs: vector of center frequencies
%nfibers: number of fibers per CF
%ANpar: AN parameters

%psth: time-by-CFs-by-nfibers-by-nch output
%x: x-axis for psth in seconds

%created by Luke Baltzell 04/28/21

if nargin == 5
    ANpar.nrep = 1;       %number of repetitions for the psth
    ANpar.cohc = 1;       %OHC scaling factor: 1 is normal OHC function; 0 is complete OHC dysfunction
    ANpar.cihc = 1;       %IHC scaling factor: 1 is normal IHC function; 0 is complete IHC dysfunction
    ANpar.species = 2;    %human with BM tuning from Shera et al. (PNAS 2002)
    ANpar.noiseType = 1;  %1 for variable fGn and 0 for fixed (frozen) fGn
    ANpar.implnt = 0;     %implnt is for "approxiate" or "actual" implementation of the power-law functions: "0" for approx. and "1" for actual implementation
    ANpar.spont = 60;     %spontanteous firing rate of fiber (60 spikes/sec is relatively high)
    ANpar.tabs = 0.0007;  %tabs is the absolute refractory period in s
    ANpar.trel = 0.0006;  %trel is the baselines mean relative refractory period in s
end

%make sure stimuli are column vectors
[nch,dim] = min(size(stim));
if dim == 1
    stim = stim';
end

%determine if stereo and make sure length(dB) matches nch
if nch == 2 
    sflg = 1;
    if length(dB) == 1
        dB(2) = dB;
    end
else
    sflg = 0;
end

%resample to 100000 Hz if necessary
fs = 100000;
dt = 1/fs;
if stim_fs ~= fs
    stim = resample(stim,fs,stim_fs);
end
Ts = length(stim)/fs;
reptime = Ts+0.005; %time between stimulus repetitions in seconds (or duration of single repetition)
D = round(reptime*fs); %duration of psth
x = [1:D]/fs; %x-axis for psth in seconds

%convert to Pa from dB
p0 = 0.00002;
Pa = p0*10.^(dB/20);
if rms(stim) == 0
    pin = stim;
else
    pin = stim.*(Pa./rms(stim));
end

if sflg == 1
    psth = zeros(D,length(CFs),nfibers,2);
    for f = 1:length(CFs)
        for n = 1:nfibers
            for i = 1:2 %left/right
                vihc = model_IHC_BEZ2018(pin(:,i)',CFs(f),ANpar.nrep,dt,reptime,ANpar.cohc,ANpar.cihc,ANpar.species);
                psth(:,f,n,i) = model_Synapse_BEZ2018(vihc,CFs(f),ANpar.nrep,dt,ANpar.noiseType,ANpar.implnt,ANpar.spont,ANpar.tabs,ANpar.trel);
            end
        end
    end
else
    psth = zeros(D,length(CFs),nfibers);
    for n = 1:nfibers
        for f = 1:length(CFs)
            vihc = model_IHC_BEZ2018(pin',CFs(f),ANpar.nrep,dt,reptime,ANpar.cohc,ANpar.cihc,ANpar.species);
            psth(:,f,n) = model_Synapse_BEZ2018(vihc,CFs(f),ANpar.nrep,dt,ANpar.noiseType,ANpar.implnt,ANpar.spont,ANpar.tabs,ANpar.trel);
        end
    end
end

end