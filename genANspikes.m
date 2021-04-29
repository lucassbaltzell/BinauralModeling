function psth = genANspikes(stim,stim_fs,dB,CFs,nfibers)
%generates post-stimulus spike histograms for a vector of center
%frequencies and for a set of fibers corresponding to each center
%frequency. In this function we call the Bruce, Erfani & Zilany model 
%described in Bruce et al. (2018): https://doi.org/10.1016/j.heares.2017.12.016
%AN parameters are fixed

%psth: time-by-CFs-by-nfibers-by-nch output
%stim: stimulus, could be mono or stereo
%stim_fs: sampling rate of stimulus
%dB: desired intensity of stimulus in dB. For stereo stimulus, this can
%have two values (left and right)
%CFs: vector of center frequencies
%nfibers: number of fibers per CF

%created by Luke Baltzell 04/23/21. This function has been replaced by
%genBEZpsth, which allows AN parameters to be given as input

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

%convert to Pa from dB
p0 = 0.00002;
Pa = p0*10.^(dB/20);
pin = stim.*(Pa./rms(stim));

%set AN parameters
nrep = 1;       %number of repetitions for the psth
reptime = Ts+0.005; %time between stimulus repetitions in seconds (or duration of single repetition)
cohc = 1;       %OHC scaling factor: 1 is normal OHC function; 0 is complete OHC dysfunction
cihc = 1;       %IHC scaling factor: 1 is normal IHC function; 0 is complete IHC dysfunction
species = 2;    %human with BM tuning from Shera et al. (PNAS 2002)

noiseType = 1;  %1 for variable fGn and 0 for fixed (frozen) fGn
implnt = 0;     %implnt is for "approxiate" or "actual" implementation of the power-law functions: "0" for approx. and "1" for actual implementation
spont = 50;     %spontanteous firing rate of fiber (50 spikes/sec is relatively high)
tabs = 0.0007;  %tabs is the absolute refractory period in s
trel = 0.0006;  %trel is the baselines mean relative refractory period in s

D = round(reptime*fs);

if sflg == 1
    psth = zeros(D,length(CFs),nfibers,2);
    for f = 1:length(CFs)
        for n = 1:nfibers
            for i = 1:2 %left/right
                vihc = model_IHC_BEZ2018(pin(:,i)',CFs(f),nrep,dt,reptime,cohc,cihc,species);
                psth(:,f,n,i) = model_Synapse_BEZ2018(vihc,CFs(f),nrep,dt,noiseType,implnt,spont,tabs,trel);
            end
        end
    end
else
    psth = zeros(D,nfibers,length(CFs));
    for n = 1:nfibers
        for f = 1:length(CFs)
            vihc = model_IHC_BEZ2018(pin',CFs(f),nrep,dt,reptime,cohc,cihc,species);
            psth(:,f,n) = model_Synapse_BEZ2018(vihc,CFs(f),nrep,dt,noiseType,implnt,spont,tabs,trel);
        end
    end
end

end