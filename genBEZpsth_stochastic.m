function [psth,x] = genBEZpsth_stochastic(stim,stim_fs,dB,CFs,nfibers,ANpar)
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
%ANpar: AN parameters. By default, ANpar calls for a mix of spontaneous
%rates and a distribution of refractory periods.

%psth: time-by-CFs-by-nfibers-by-nch output
%x: x-axis for psth

%created by Luke Baltzell 04/28/21

if nargin == 5
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

if range([length(ANpar.spont.mu) length(ANpar.spont.sd) length(ANpar.spont.lims)]) ~= 0
    error('make sure spont parameters are consistent')
end

%check if ANpar.spont calls for a single sr type or for a mix, and get "sr"
%distribution
if length(ANpar.spont.mu) == 1
    sr = normrnd(ANpar.spont.mu,ANpar.spont.sd,length(CFs),nfibers);
    sr(sr < ANpar.spont.lims(1)) = ANpar.spont.lims(1);
    sr(sr > ANpar.spont.lims(2)) = ANpar.spont.lims(2);
elseif size(ANpar.spont.mu,2) == 3
    N = length(CFs)*nfibers;
    %low spont
    sr_low = normrnd(ANpar.spont.mu(1),ANpar.spont.sd(1),ceil(N*0.16),1);
    sr_low(sr_low < ANpar.spont.lims(1,1)) = ANpar.spont.lims(1,1);
    sr_low(sr_low > ANpar.spont.lims(1,2)) = ANpar.spont.lims(1,2);
    %med spont
    sr_med = normrnd(ANpar.spont.mu(2),ANpar.spont.sd(2),ceil(N*0.23),1);
    sr_med(sr_med < ANpar.spont.lims(2,1)) = ANpar.spont.lims(2,1);
    sr_med(sr_med > ANpar.spont.lims(2,2)) = ANpar.spont.lims(2,2);
    %high spont
    sr_high = normrnd(ANpar.spont.mu(3),ANpar.spont.sd(3),ceil(N*0.61),1);
    sr_high(sr_high < ANpar.spont.lims(3,1)) = ANpar.spont.lims(3,1);
    sr_high(sr_high > ANpar.spont.lims(3,2)) = ANpar.spont.lims(3,2);
    %concatenate and reshape
    sr = cat(1,sr_low,sr_med,sr_high);  %concatenate
    sr = sr(randperm(length(sr)));      %randomize
    sr = sr(1:N);                       %resize
    sr = reshape(sr,length(CFs),nfibers);%reshape
end

%get tabs and trel distribution
tinds = rand(length(CFs),nfibers); %generate random values for tabs and trel
tabs = tinds*(ANpar.tabs_lims(2)-ANpar.tabs_lims(1))+ANpar.tabs_lims(1);
trel = tinds*(ANpar.trel_lims(2)-ANpar.trel_lims(1))+ANpar.trel_lims(1);

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
x = [1:D]/fs; %x-axis for psth

%convert to Pa from dB
p0 = 0.00002;
Pa = p0*10.^(dB/20);
pin = stim.*(Pa./rms(stim));

if sflg == 1
    psth = zeros(D,length(CFs),nfibers,2);
    for f = 1:length(CFs)
        for n = 1:nfibers
            for i = 1:2 %left/right
                vihc = model_IHC_BEZ2018(pin(:,i)',CFs(f),ANpar.nrep,dt,reptime,ANpar.cohc,ANpar.cihc,ANpar.species);
                psth(:,f,n,i) = model_Synapse_BEZ2018(vihc,CFs(f),ANpar.nrep,dt,ANpar.noiseType,ANpar.implnt,sr(f,n),tabs(f,n),trel(f,n));
            end
        end
    end
else
    psth = zeros(D,length(CFs),nfibers);
    for n = 1:nfibers
        for f = 1:length(CFs)
            vihc = model_IHC_BEZ2018(pin',CFs(f),ANpar.nrep,dt,reptime,ANpar.cohc,ANpar.cihc,ANpar.species);
            psth(:,f,n) = model_Synapse_BEZ2018(vihc,CFs(f),ANpar.nrep,dt,ANpar.noiseType,ANpar.implnt,sr(f,n),tabs(f,n),trel(f,n));
        end
    end
end

end