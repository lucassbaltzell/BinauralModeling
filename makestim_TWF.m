function y = makestim_TWF(stim,fs,winlen,n,D)

stim = stim/rms(stim); %normalize stim
stim = stim(1:round(winlen*n*fs)); %resize corresponding to number of bins

%generate lateralized copies of stimulus
X = zeros(length(stim),2,length(D));
for i = 1:length(D)
    if D(i) <= 0
        %delay right channel relative to left
        X(:,1,i) = stim;
        X(:,2,i) = FFTdelay(stim,abs(D(i))*1e-3,fs); %convert ms to s and apply delay
    else
        %delay left channel relative to right
        X(:,1,i) = FFTdelay(stim,abs(D(i))*1e-3,fs);
        X(:,2,i) = stim;
    end
end

%concatenate bins
tc = 0.005; %crossfade duration (seconds)
y = squeeze(X(:,:,1));
for i = 1:length(D)-1
    s0 = round((winlen*i)*fs);
    y = crossfade(y,squeeze(X(:,:,i+1)),s0,fs,tc);
end

end