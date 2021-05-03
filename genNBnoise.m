function [y,flims] = genNBnoise(dur,fs,cf,bw,tflg)

if nargin == 4
    tflg = 0;
end

noise = randn(1,dur*fs);

if tflg == 0
    noise_fft = fft(noise);
    nfft = length(noise_fft);
    noise_fft_shft = fftshift(noise_fft);
    if mod(nfft,2) == 1
        error('Change stimulus duration')
    end
    f = (1/dur)*[1:1:nfft/2];
    pbnd = zeros(1,length(f));
    flims = [cf/2^(1/(2*bw)) cf*2^(1/(2*bw))];
    lc_ind = round(dur*flims(1));
    uc_ind = round(dur*flims(2));
    pbnd(lc_ind:uc_ind) = 1;
    pbnd_all = cat(2,fliplr(pbnd(1:end-1)),0,pbnd);
    y_fft = ifftshift(noise_fft_shft.*pbnd_all);
    y = ifft(y_fft,'symmetric')';
    
elseif tflg == 1
    %write time domain filter
end

end