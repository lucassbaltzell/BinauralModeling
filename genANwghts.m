%This script simulates temporal weighting functions for "two". This
%simulation is based on simulated spike trains at the output of the BEZ
%phenomenolgical auditory nerve (AN) model. See Baltzell et al. (under 
%review) for further details

%load stimulus
[stim,stim_fs] = audioread('BUG_T16_1_3.wav'); %"two"

%define AN parameters
fs = 100000; %sampling rate required for AN model
cf_n = 40; %number of AN center frequencies
CFs = logspace(log10(200),log10(10000),cf_n); %define center frequencies between 200 and 10000 Hz
dB = 70; %desired intensity of stimulus
Nfibers = 25; %desired number of fibers for each CF
Ntrials = 250; %desired number of trials over which to obtain simulated lateralization judgemnets

%resample stimulus prior to spatialization
stim = resample(stim,fs,stim_fs);
dur = length(stim)/fs;

%generate lateralized stimuli
winlen = 0.03; %duration of temporal bin
nbins = floor(dur./winlen);
dly_rng = 0.3; %range of ITDs in ms (+- 150 us)
tStim = cell(1,Ntrials);
tDlys = zeros(nbins,Ntrials);
for n = 1:Ntrials
    dlys = dly_rng*rand(nbins,1) - dly_rng/2;  
    y = makestim_TWF(stim,fs,winlen,nbins,dlys);
    y = cat(1,y,zeros(0.015*fs,2)); %add 15 ms zero-pad to allow auditory nerve to fire after stimulus has ended
    tDlys(:,n) = dlys;
    tStim{1,n} = y;
end

binwidth_t = 20; %in s (20 us)
binwidth = binwidth_t*(fs/1e6); %samples
for t = 1:Ntrials
    tstim = tStim{1,t};
    psth = genANspikes_stochastic(tstim,fs,[dB dB],CFs,Nfibers,0);
    SCCwght = zeros(size(psth,1),length(CFs));
    for f = 1:cf_n
        xl = squeeze(psth(:,f,:,1));
        xr = squeeze(psth(:,f,:,2));
        [SCC,xax] = getSCC(xl,xr,binwidth,fs);
        SCCwght(:,f) = centralityWeighting1D(xax,SCC,CFs(f),'stern','pdf',5);
    end
    SCCall = mean(SCCwght,2)'; %average over frequency
    centroid(t) = sum(SCCall.*xax*binwidth_t)/sum(SCCall*binwidth_t); %get centroid
end

%run linear regression to obtain weights
y = centroid';
x = tDlys';
mdl = fitlm(x,y);
Beta = mdl.Coefficients{2:end,1};
SE = mdl.Coefficients{2:end,2};
Beta_norm = Beta/sum(Beta);
SE_norm = SE/sum(Beta);

%plot weights
binlims = [0:winlen:winlen*15]*1e3;
xbin = mean([binlims(1:end-1); binlims(2:end)]);
figure
p1 = errorbar(xbin,Beta_norm,SE_norm,'ko','markerfacecolor','k','markersize',8);
xticks(binlims)
xlim([binlims(1)-(winlen*1e3/2) binlims(end)+((winlen*1e3/2))])
xlabel('time (ms)','fontsize',12)
ylim([0 0.25])
ylabel('normalized weight (a.u.)','fontsize',12)
title('TWF for "two"','fontsize',14)
set(gca,'linewidth',2)