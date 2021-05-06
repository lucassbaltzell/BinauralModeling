function [y,flims] = genNBrustle(dur,fs,gpwidth,cf,bw,tflg,nord)
%This function creates a "rustle" noise token, with a sparseness determined
%by gpwdth, and filtered according to cf, bw, and tflg
%OUTPUT
%y = rustle noise token
%INPUT
%dur = desired duration of noise token
%fs = sampling rate of noise token
%gpwidth = mean gap width between noise samples in microseconds
%cf = center frequency of noise
%bw = bandwidth in octave units (e.g. 1/3)
%tflag = type of filter. 0 for brickwall frequency domain filter, 1 for
%time comain filter of order nord

%created by Luke Baltzell on 05/04/21

if nargin == 5
    tflg = 0;
end

zpad = 0.01*fs; %zero-pad to allow for high-order bandpass filtering
rnoise = genRustle(dur,fs,gpwidth);
rnoise = cat(1,zeros(zpad,1),rnoise);

flims = [cf/2^(1/(2/bw)) cf*2^(1/(2/bw))];
if tflg == 0
    %brickwall frequency domain filter
    noise_fft = fft(rnoise);
    nfft = length(noise_fft);
    noise_fft_shft = fftshift(noise_fft);
    f = (1/dur)*[1:1:floor(nfft/2)];
    pbnd = zeros(1,length(f));
    lc_ind = round(dur*flims(1));
    uc_ind = round(dur*flims(2));
    pbnd(lc_ind:uc_ind) = 1;
    pbnd_all = cat(2,fliplr(pbnd(1:end-1)),0,pbnd);
    y_fft = ifftshift(noise_fft_shft.*pbnd_all);
    y = ifft(y_fft,'symmetric')';
elseif tflg == 1
    N = round(nord/2); %forwards-backwards (0-phase) filter so halve nord
    Nq = fs/2;
    Wn = [flims(1) flims(2)]/Nq;
    [z,p,k] = butter(N,Wn);
    [sos,g] = zp2sos(z,p,k);
    y = filtfilt(sos,g,rnoise);
    y = y(zpad+1:end);
end

end