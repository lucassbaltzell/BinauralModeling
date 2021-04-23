function psth = genANspikes_stochastic(stim,stim_fs,dB,CFs,nfibers,srflg)
%generates post-stimulus spike histograms for a vector of center
%frequencies and for a set of fibers corresponding to each center
%frequency. In this function we call the Bruce, Erfani & Zilany (BEZ) model 
%described in Bruce et al. (2018): https://doi.org/10.1016/j.heares.2017.12.016
%In this function, some AN parameters are stochastic, providing a more 
%realistic simulation of a population AN response. Following Bruce et al.
%(2018), we draw a fully correlated pair of tabs and trel values for each
%fiber, and following Klug et al. (2020), who use the same AN front end for
%their LSO-inspired binarual model (https://doi.org/10.1121/10.0001602), we
%let these draws be fully correlated across ears for binarual input

%psth: time-by-CFs-by-nfibers-by-nch output
%stim: stimulus, could be mono or stereo
%stim_fs: sampling rate of stimulus
%dB: desired intensity of stimulus in dB. For stereo stimulus, this can
%have two values (left and right)
%CFs: vector of center frequencies
%nfibers: number of fibers per CF
%srflg: indicates spontaneous rate preference for fibers. "0" indicates mix,
%"1" indicates low, "2" indicates medium, and "3" indicates high. Following
%Liberman (1978), we use 16% low, 23% med, and 61% high spontaneous rate
%for srflg = 0

%created by Luke Baltzell 04/23/21

if nargin == 5
    srflg = 0;
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
if srflg == 1
    spont.mu = 0.1; %mean spontanteous firing rate (sr) of fiber
    spont.sd = 0.1; %standard deviation of sr
    spont.lims = [1e-3 0.2]; %range of acceptable sr values
elseif srflg == 2
    spont.mu = 4;
    spont.sd = 4;
    spont.lims = [0.2 18];
elseif srflg == 3
    spont.mu = 70;
    spont.sd = 30;
    spont.lims = [18 180];
else
    spont.mu = [0.1 4 70];
    spont.sd = [0.1 4 30];
    spont.lims = [1e-3 0.2; 0.2 18; 18 180];
end
tabs_lims = [0.0002085 0.0006915];  %tabs is the range of absolute refractory period in s
trel_lims = [0.000131 0.000894];  %trel is the range of relative refractory period in s

D = round(reptime*fs); %number of samples for the psth
tinds = rand(length(CFs),nfibers); %generate random values for tabs and trel
if srflg ~= 0
    sr = normrnd(spont.mu,spont.sd,length(CFs),nfibers);
    sr(sr < spont.lims(1)) = spont.lims(1);
    sr(sr > spont.lims(2)) = spont.lims(2);
else
    N = length(CFs)*nfibers;
    %low spont
    sr_low = normrnd(spont.mu(1),spont.sd(1),ceil(N*0.16),1);
    sr_low(sr_low < spont.lims(1,1)) = spont.lims(1,1);
    sr_low(sr_low > spont.lims(1,2)) = spont.lims(1,2);
    %med spont
    sr_med = normrnd(spont.mu(2),spont.sd(2),ceil(N*0.23),1);
    sr_med(sr_med < spont.lims(2,1)) = spont.lims(2,1);
    sr_med(sr_med > spont.lims(2,2)) = spont.lims(2,2);
    %high spont
    sr_high = normrnd(spont.mu(3),spont.sd(3),ceil(N*0.61),1);
    sr_high(sr_high < spont.lims(3,1)) = spont.lims(3,1);
    sr_high(sr_high > spont.lims(3,2)) = spont.lims(3,2);
    %concatenate and reshape
    sr = cat(1,sr_low,sr_med,sr_high);  %concatenate
    sr = sr(randperm(length(sr)));      %randomize
    sr = sr(1:N);                       %resize
    sr = reshape(sr,length(CFs),nfibers);%reshape
end
if sflg == 1
    psth = zeros(D,length(CFs),nfibers,2);
    for f = 1:length(CFs)
        for n = 1:nfibers
            for i = 1:2 %left/right
                tabs = tinds(f,n)*(tabs_lims(2)-tabs_lims(1))+tabs_lims(1);
                trel = tinds(f,n)*(trel_lims(2)-trel_lims(1))+trel_lims(1);
                vihc = model_IHC_BEZ2018(pin(:,i)',CFs(f),nrep,dt,reptime,cohc,cihc,species);
                psth(:,f,n,i) = model_Synapse_BEZ2018(vihc,CFs(f),nrep,dt,noiseType,implnt,sr(f,n),tabs,trel);
            end
        end
    end
else
    psth = zeros(D,nfibers,length(CFs));
    for n = 1:nfibers
        for f = 1:length(CFs)
            tabs = tinds(f,n)*(tabs_lims(2)-tabs_lims(1))+tabs_lims(1);
            trel = tinds(f,n)*(trel_lims(2)-trel_lims(1))+trel_lims(1);
            vihc = model_IHC_BEZ2018(pin',CFs(f),nrep,dt,reptime,cohc,cihc,species);
            psth(:,f,n) = model_Synapse_BEZ2018(vihc,CFs(f),nrep,dt,noiseType,implnt,spont,tabs,trel);
        end
    end
end

end