function [centroid,SCCall,xax] = computeSCC(stim,stim_fs,dB,CFs,Npsth)

[n,m] = size(stim);
if n < 2 || m < 2
    error('Input should be stereo')
end

if m > 2
    stim = stim';
end

%resample
fs = 100000;
dt = 1/fs;
stim = resample(stim,fs,stim_fs);
Ts = length(stim)/fs;

%convert to Pa
p0 = 0.00002;
Pa = p0*10^(dB/20);

pinL = stim(:,1)'.*(Pa/rms(stim(:,1)));
pinR = stim(:,2)'.*(Pa/rms(stim(:,2)));

nrep = 1;       %number of repetitions for the psth
reptime = Ts+0.005; %time between stimulus repetitions in seconds
cohc = 1;       %OHC scaling factor: 1 is normal OHC function; 0 is complete OHC dysfunction
cihc = 1;       %IHC scaling factor: 1 is normal IHC function; 0 is complete IHC dysfunction
species = 2;    %human with BM tuning from Shera et al. (PNAS 2002)

noiseType = 1;  %1 for variable fGn and 0 for fixed (frozen) fGn
implnt = 0;     %implnt is for "approxiate" or "actual" implementation of the power-law functions: "0" for approx. and "1" for actual implementation
spont = 50;
tabs = 0.0007;  %tabs is the absolute refractory period in s
trel = 0.0006;  %trel is the baselines mean relative refractory period in s

D = round(reptime*fs);
binwidth_t = 20; %in us
binwidth = binwidth_t*(fs/1e6); %samples
edges = [-D:binwidth:D];
xax = edges(1:end-1)*(1e6/fs); %us

% SCCall = zeros(length(edges)-1,iters);
psth_L = zeros(D,Npsth);
psth_R = zeros(D,Npsth); 
SCCwght = zeros(D,length(CFs));
for f = 1:length(CFs)
    for i = 1:Npsth
        vihc_L = model_IHC_BEZ2018(pinL,CFs(f),nrep,dt,reptime,cohc,cihc,species);
        vihc_R = model_IHC_BEZ2018(pinR,CFs(f),nrep,dt,reptime,cohc,cihc,species);
        psth_L(:,i) = model_Synapse_BEZ2018(vihc_L,CFs(f),nrep,dt,noiseType,implnt,spont,tabs,trel);
        psth_R(:,i) = model_Synapse_BEZ2018(vihc_R,CFs(f),nrep,dt,noiseType,implnt,spont,tabs,trel);
    end
    SCC = SCCfun(psth_L,psth_R,Npsth,binwidth,D,edges);
    SCCwght(:,f) = centralityWeighting1D(xax, SCC, CFs(f), 'stern', 'pdf', 5);
end
SCCall = mean(SCCwght,2)';
centroid = sum(SCCall.*xax*binwidth_t)/sum(SCCall*binwidth_t); %get centroid
end
% % pd = fitdist(ITDdist','Normal');
% pd = fitdist(CENTdist','Normal');
% % pkwght = mean(SCCmax);
% SCCmean = mean(SCCall,2); %in microseconds (use histogram to plot)
% x = [-1000:20:1000];
% y = pdf(pd,[-1000:20:1000]);
% y = y/rms(y);
% end


function SCC = SCCfun(psth1,psth2,N,binwidth,D,edges)
%computes the SCC in one direction. Specifically, it computes the 
%difference (in samples) between each spike in psth1 and all spikes in 
%psth2. If the spikes in psth2 are delayed relative to the spikes in psth1
%(i.e. left-leading), the SCC should be shifted in the negative direction.
%See Louage et al. (2004) for details, and please cite if this code is used
%for publication: http://www.physiology.org/doi/10.1152/jn.00816.2003

%psth1: reference psth (post-stimulus time histogram)
%psth2: psth referenced to psth1
%number of fibers in each psth (should be the same)
r1 = mean(mean(psth1));
r2 = mean(mean(psth2));
SCC = [];
for n = 1:N
    for nn = 1:N
        st1 = find(psth1(:,n) == 1);
        st2 = find(psth2(:,nn) == 1);
        for i = 1:length(st1)
            temp(i,:) = st1(i)-st2;
        end
        temprs = reshape(temp,length(st1)*length(st2),1);
        temp = [];
        SCC = cat(1,SCC,temprs);
    end
end
%SCC = cat(1,SCC,-SCC);
%norm = N^2*r1*r2*binwidth*D*2;
norm = N^2*r1*r2*binwidth*D; %removed "*2" to account for asymmetric SCC
histvec = histcounts(SCC,edges);
SCC = histvec./norm;
end