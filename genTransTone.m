function y = genTransTone(fs,dur,cf,mf,tc,phi)
% This function returns a "transposed" tone, following Bernstein &
% Trahiotis (2002): https://doi.org/10.1121/1.1497620
%OUPUTS
%y: Transposed tone
%INPUTS
%fs: sampling rate
%dur: duration of SAM tone
%cf: carrier frequency
%mod: modulation frequency
%tc: duration of ramp/damp
%phi: desired starting phase of sine waves. The first element should refer
%to the phase of the carrier, and the second element to the modulator

%created by Luke Baltzell, modified 04/30/21

if nargin == 4
    tc = 0;
    phi = [0 0];
elseif nargin == 5
    phi = [0 0];
end

t = [1/fs:1/fs:dur];
car = sin(2*pi*cf*t + phi(1));
tmod = sin(2*pi*mf*t + phi(2));
tmod(tmod < 0) = 0; %hlf-wave rectify modulator

%apply boxcar filter to modulator in frequency domain
tmod_fft = fft(tmod);
nfft = length(tmod_fft);
tmod_fft_shft = fftshift(tmod_fft);
f = (1/dur)*[1:1:floor(nfft/2)]; %vector of positive frequencies (no dc)
pbnd = zeros(1,length(f));
flim = 2000; %boxcar low-pass filter at 2 kHz
uc_ind = round(dur*flim);
pbnd(1:uc_ind) = 1;
pbnd_all = cat(2,fliplr(pbnd(1:end-1)),1,pbnd); %concatenate neg freqs with dc with pos freqs
tmod_fft_filt = ifftshift(tmod_fft_shft.*pbnd_all);
tmod_filt = ifft(tmod_fft_filt,'symmetric')';

amtone = car.*tmod_filt';
if tc ~= 0
    amtone = rampdamp(amtone,tc,fs);
end
y = amtone./rms(amtone);

end