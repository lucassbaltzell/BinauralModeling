function [vs,CFs] = getVS(stim,stim_fs,dB,tf,CFs,pflg)
%calculates vector strength at the output of the BEZ auditory nerve model
%for a periodic stimulus (tone or transposed tone)

%stim: input stimulus
%stim_fs: sampling rate of stimulus
%dB: desired intensity of stimulus
%tf: target frequency (frequency against which vector strength is computed)
%CFs: range of AN center frequencies over which vector strength is computed
%pflg: plot figure if 1, do not plot if 0

%vs: vector of vector strengths for each CF
%CFs: returns CFs to plot vs if CFs not specified in argin

%created by Luke Baltzell, modified 04/28/21

if nargin == 4
    CFs = logspace(log10(200),log10(10000),60);
    pflg = 0
elseif nargin == 5
    pflg = 0
end
    

%check if stereo input
[n,m] = size(stim);
if n > 1 && m > 1
    error('Input should be mono')
end

%change to row vector
if n > 1
    stim = stim';
end
clear n m

%set AN parameters
nfibers = 1;
ANpar.nrep = 500;    %number of repetitions for the psth
ANpar.cohc = 1;       %OHC scaling factor: 1 is normal OHC function; 0 is complete OHC dysfunction
ANpar.cihc = 1;       %IHC scaling factor: 1 is normal IHC function; 0 is complete IHC dysfunction
ANpar.species = 2;    %human with BM tuning from Shera et al. (PNAS 2002)
ANpar.noiseType = 1;  %1 for variable fGn and 0 for fixed (frozen) fGn
ANpar.implnt = 0;     %implnt is for "approxiate" or "actual" implementation of the power-law functions: "0" for approx. and "1" for actual implementation
ANpar.spont = 60;     %spontanteous firing rate of fiber (60 spikes/sec is relatively high)
ANpar.tabs = 0.0007;  %tabs is the absolute refractory period in s
ANpar.trel = 0.0006;  %trel is the baselines mean relative refractory period in s

%generate psth
[PSTH,x] = genBEZpsth(stim,stim_fs,dB,CFs,nfibers,ANpar);

%get vector strength
omega = 2*pi*tf; %define angular velocity (cycles/second) against which to compare spike times
for f = 1:length(CFs)
    psth = squeeze(PSTH(:,f));
    spike_ind = find(psth > 0); %find all bins with spikes > 0
    n_ind = 0;
    for i = 1:length(spike_ind)
        t = x(spike_ind(i)); %time of spike in seconds
        N = psth(spike_ind(i)); %number of spikes at time t
        for n = 1:N
            n_ind = n_ind + 1;
            y(n_ind) = exp(1i*omega*t); %calculate phase of omega with respect to t, for N spikes
        end
    end   
    vs(f) = abs(sum(y)/n_ind); %normalize phase vector sum and get magnitude
end

if pflg == 1
    figure
    semilogx(CFs,vs,'bo-','linewidth',1.5)
    xlabel('CF (Hz)')
    ylabel('vector strength')
    set(gca,'fontsize',12)
    title(['reference: ' num2str(tf) ' Hz'])
end

end
