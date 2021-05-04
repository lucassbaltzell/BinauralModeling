function y = genRustle(dur,fs,gpwidth)
%This function creates a "rustle" noise token, with a sparseness determined
%by gpwdth. See Ewert et al. (2012) for details: https://doi.org/10.1007/s10162-011-0303-2
%OUTPUT
%y = rustle noise token
%INPUT
%dur: desired duration of noise token
%fs: sampling rate of noise token
%gpwidth: mean gap width between noise samples in microseconds

%created by Luke Baltzell, modified 05/03/21

if nargin == 2
    gpwidth = 5800;
end

noise = randn(dur*fs,1); %generate base noise token
gps = 1e-6*(gpwidth*rand(1,dur*fs)); %generate base gap durations in seconds
nind = 0; %index for noise samples
gind = 0; %index for gap samples
sflg = 0; %stop flag
while sflg == 0
    nind = nind+1;
    gind = gind+1;
    rnoise(gind,1) = noise(nind);
    gpsamps = round(gps(nind)*fs);
    if gpsamps ~= 0
        rnoise = cat(1,rnoise,zeros(gpsamps,1));
        gind = gind+gpsamps;
    end
    if length(rnoise) >= length(noise)
        sflg = 1;
    end
end

y = rnoise(1:length(noise));

end