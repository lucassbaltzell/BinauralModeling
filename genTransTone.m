function y = genTransTone(fs,dur,cf,mf)

t = [1/fs:1/fs:dur];

car = sin(2*pi*cf*t);

tmod = sin(2*pi*mf*t);
tmod(tmod < 0) = 0;

tmod_fft = fft(tmod);
nfft = length(tmod_fft);
tmod_fft_shft = fftshift(tmod_fft);
if mod(nfft,2) == 1
    error('Change stimulus duration')
end

f = (1/dur)*[1:1:nfft/2];
pbnd = zeros(1,length(f));
flim = 2000;
% lc_ind = round(dur*flims(1));
uc_ind = round(dur*flim);
pbnd(1:uc_ind) = 1;
pbnd_all = cat(2,fliplr(pbnd(1:end-1)),1,pbnd);
tmod_fft_filt = ifftshift(tmod_fft_shft.*pbnd_all);
tmod_filt = ifft(tmod_fft_filt,'symmetric')';

amt = car.*tmod_filt';
y = amt./rms(amt);

end