function y = genRustle(dur,fs,gpwidth)
%gpwidth in microseconds

noise = randn(dur*fs,1);
gps = 1e-6*(gpwidth*rand(1,dur*fs)); %in seconds
% gps = 1e-6*(gpwidth*rand(1,dur*fs)); %in seconds
% temp = round(gps*fs);
% save('gpsamps','temp')
ind = 0;
n = 0;
sflg = 0;
% for n = 1:dur*fs
while sflg == 0
    n = n+1;
    ind = ind+1;
    rnoise(ind,1) = noise(n);
    gpsamps = round(gps(n)*fs);
    if gpsamps ~= 0
        rnoise = cat(1,rnoise,zeros(gpsamps,1));
        ind = ind+gpsamps;
    end
    if length(rnoise) >= length(noise)
        sflg = 1;
    end
end

y = rnoise(1:length(noise));

end






