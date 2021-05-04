function [y,flims] = genNBnoise(dur,fs,cf,bw,tflg)
%This function computes narrowband noise with center frequency 'cf' and
%bandwidth 'bw' in octave units (e.g. 1/3).
%OUTPUTS
%y = narrowband noise
%flims = lower and upper cutoff in frequency corresponding to bw
%INPUTS
%dur = desired duration of noise in seconds
%fs = sampling rate
%cf = center frequency of noise
%bw = bandwidth in octave units
%tflg = flag determining whether to use a brickwall frequency-domain filter
%or a time domain filter. If 0, frequency-domain filter, and if 1, time
%domain filter

%created by Luke Baltzell, modified 05/03/21

if nargin == 4
    tflg = 0;
end

noise = randn(1,dur*fs);
flims = [cf/2^(1/(2/bw)) cf*2^(1/(2/bw))];
if tflg == 0
    %brickwall frequency domain filter
    noise_fft = fft(noise);
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
    %16th-order time domain filter
    N = 8; %forwards-backwards (0-phase) filter so double order
    Nq = fs/2;
    Wn = [flims(1) flims(2)]/Nq;
    [z,p,k] = butter(N,Wn);
    [sos,g] = zp2sos(z,p,k);
    y = filtfilt(sos,g,noise);
end

end