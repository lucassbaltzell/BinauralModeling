%This script integrates SCCs over frequency before computing centroid
pth = pwd;

%simulate two thresholds
fs = 100000;
[stim,stim_fs] = audioread('MBUG50K_T16_1_3_F_two.wav');
% [stim,stim_fs] = audioread('MBUG50K_T16_5_3_F_six.wav');
% [stim,stim_fs] = audioread('MBUG50K_T16_6_3_F_eight.wav');
% [stim,stim_fs] = audioread('MBUG50K_T16_7_3_F_nine.wav');
stim = resample(stim,fs,stim_fs);
dur = length(stim)/fs;
t = [1/fs:1/fs:dur];
% tc = 0.025;
dB = 70;

% %list of delay anchor points (including 0us reference)
% dlys = [0 10 20 40 80 160 320 640 1280]*1e-6;
cf_n = 40;
cfvec = logspace(log10(200),log10(10000),cf_n);

addpath('/Users/lukebaltzell/Documents/UR_EAR_2020b')
Npsth = 35;

Ntrials = 250;
% Ntrials = Ntrials/2;
%gen trial stim
params.tardB = 70;
params.refdB = 133;
params.stim = stim;
params.winlen = 0.03;
% params.winlen = 0.4;
params.n = floor(dur./params.winlen);
params.fs = fs;
sigma = 0.3;
tStim = cell(1,Ntrials);
tDlys = zeros(params.n,Ntrials);
for n = 1:Ntrials
    dlys = sigma*rand(params.n,1) - sigma/2;
    y = makestim_TWF(params,dlys);
    y = cat(1,y,zeros(0.03*fs,2));
    tDlys(:,n) = dlys;
    tStim{1,n} = y;
end

cd('/Users/lukebaltzell/Documents/UR_EAR_2020b/Luke_ITDmodelling')
% [~,D] = genANspikes(randn(1,dur*fs),fs,dB,500); %dummy spike gen to get params
[~,D] = genANspikes(randn(1,length(y)),fs,dB,500); %dummy spike gen to get params
binwidth_t = 20; %in us
binwidth = binwidth_t*(fs/1e6); %samples
edges = [-D:binwidth:D];
xax = edges(1:end-1)*(1e6/fs); %us

% ANset = cell(1,Nsets);
for t = 1:Ntrials
    tstim = tStim{1,t};
%     ANspikes = zeros(Npsth,D,2,cf_n);
%     for n = 1:Npsth
%         %generate AN spikes
%         ANspikes(n,:,:,:) = genANspikes(tstim,fs,dB,cfvec,1);
%     end
%     SCC_wght = zeros(cf_n,D);
%     for f = 1:cf_n
%         AN_L = squeeze(ANspikes(:,:,1,f));
%         AN_R = squeeze(ANspikes(:,:,2,f));
%         SCC = genSCC(AN_L,AN_R,Npsth,D,binwidth);
%         SCC_wght(f,:) = centralityWeighting1D(xax,SCC,cfvec(f),'stern','pdf',5);
%     end
%     SCC_wght = mean(SCC_wght);
%     centroid(t) = sum(SCC_wght.*xax*binwidth_t)/sum(SCC_wght*binwidth_t);
    centroid(t) = computeSCC(tstim,fs,dB,cfvec,Npsth);
end
cd(pth)
AN_150us.y = centroid';
AN_150us.x = -tDlys';
save('AN_nine_150us','AN_150us')