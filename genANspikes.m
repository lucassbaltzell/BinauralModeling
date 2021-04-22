function [psth,D] = genANspikes(stim,stim_fs,dB,CFs,sflg)
%changed to mono option sflg LB 11/17/20
if nargin == 4
    sflg = 0;
end
% [n,m] = size(stim);
% if n > 1 || m > 1
%     error('Input should be mono')
% end

%resample
fs = 100000;
dt = 1/fs;
if stim_fs ~= fs
    stim = resample(stim,fs,stim_fs);
end
Ts = length(stim)/fs;

%convert to Pa
p0 = 0.00002;
Pa = p0*10^(dB/20);

pin = stim.*(Pa/rms(stim));

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

if sflg == 1
    psth = zeros(D,2,length(CFs));
    for f = 1:length(CFs)
        for i = 1:2 %left/right
            vihc = model_IHC_BEZ2018(pin,CFs(f),nrep,dt,reptime,cohc,cihc,species);
            psth(:,i,f) = model_Synapse_BEZ2018(vihc,CFs(f),nrep,dt,noiseType,implnt,spont,tabs,trel);
        end
    end
else
    psth = zeros(D,length(CFs));
    for f = 1:length(CFs)
        vihc = model_IHC_BEZ2018(pin,CFs(f),nrep,dt,reptime,cohc,cihc,species);
        psth(:,f) = model_Synapse_BEZ2018(vihc,CFs(f),nrep,dt,noiseType,implnt,spont,tabs,trel);
    end
end