function [y,ref] = makestim_TWF(params,D)

dBdwn = params.tardB - params.refdB;
stim = params.stim/rms(params.stim);
stim = stim(1:round(params.winlen*params.n*params.fs));

% if params.token == 1
%     % s0 = 6069;
%     s0 = round(0.1214*params.fs);
% elseif params.token == 2
%     % s0 = 7990
%     s0 = round(0.1598*params.fs);
% elseif params.token == 3
%     s0 = round(0.34*params.fs);
% end
tc = 0.005; %seconds

X = zeros(length(stim),2,length(D));
for i = 1:length(D)
    if D(i) <= 0
        X(:,1,i) = FFTdelay(stim,abs(D(i))/1000,params.fs); %convert ms to s
        X(:,2,i) = stim;
    else
        X(:,1,i) = stim;
        X(:,2,i) = FFTdelay(stim,abs(D(i))/1000,params.fs);
    end
end
ref(:,1) = stim;
ref(:,2) = stim;

% x1 = RW_bandfilter(stim,params.fs,D(:,1));
% x2 = RW_bandfilter(stim,params.fs,D(:,2));

y = squeeze(X(:,:,1));
for i = 1:length(D)-1
    s0 = round((params.winlen*i)*params.fs);
    y = crossfade(y,squeeze(X(:,:,i+1)),s0,params.fs,tc);
end

y = y*(10^(dBdwn/20));
ref = ref*(10^(dBdwn/20));
end