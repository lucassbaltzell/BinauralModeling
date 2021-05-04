%This script simulates ITD sensitivity for a pair of low-frequency and 
%high-frequency narrowband noises. This simulation is based on simulated 
%spike trains at the output of the BEZ phenomenolgical auditory nerve (AN) 
%model. This method roughly follows the procedure described by Moncada-
%Torres et al. (2018) https://doi.org/10.1121/1.5051322. The results of
%this simulation can be compared with the results of Spencer et al. (2016):
%https://doi.org/10.1121/1.4962444, among others.

%In this script, information is integrated across multiple CFs, and
%simulated ITD thresholds reflect an average sensitivity across multiple
%CFs spanning the bandwidth of the noise stimuli.

%For the sake of efficiency, rather than simulate new spike trains for each
%iteration, a set of spike trains are defined (Nsets), and are sampled with
%replacement on each iteration (Niters). Since stimuli are randomly drawn
%for each set of Nfibers in Nsets, on each iteration, a set of Nfibers
%corresponding to a particular stimulus is drawn randomly with replacement
%from Nsets

%Also for the sake of efficiency, spike trains are generated using monaural
%stimuli, and ITDs will be applied by delaying the spike trains. The high
%sample rate (100 kHz) allows for a delay resolution of 10 us, and so
%target ITD values (dlys) must be integer multiples of 10 us.

%created by Luke Baltzell for presentation at Binaural Bash 2020. Modified 
%by Luke Baltzell on 05/03/21

%define narrowband noise paramters following Bernstein 
fs = 100000;
dur = 0.3; %500 ms
tc = 0.015; %100 ms ramp
t = [1/fs:1/fs:dur];
dB = 65;
bw = 1/3; %3rd octave
cfvec_stim = [500 4000]; %center frequencies for noise stimuli
ncf_stim = length(cfvec_stim);
dlys = [0 10 20 40 80 160 320 640 1280]*1e-6; %set of delays over which to calculate sensitivity
cfvec_AN = logspace(log10(150),log10(10000),60); %generate dense vector of center frequencies
% cf_n = length(cfvec_AN);

%calculate AN cfs for each stimulus cf
bw_ext = 1/8; %extend bandwidth by 8th octave
for f = 1:ncf_stim
    [~,flims] = genNBnoise(dur,fs,cfvec_stim(f),bw);
    flims_ext = [flims(1)/2^(1/(2/bw_ext)) flims(2)*2^(1/(2/bw_ext))];
    [~,ind] = find(flims_ext(1) < cfvec_AN & flims_ext(2) > cfvec_AN);
    tmp = cfvec_AN(ind);
    if f > 1
        if length(tmp) ~= length(cf_fib)
            error('Manually adjust fiber count to be the same across noise cfs')
        end
    end
    cf_fib(f,:) = tmp;
end
ncf_fib = size(cf_fib,2);

%set AN parameters
Nfibers = 35; %desired number of fibers for each CF
Nsets = 150; %number of sets to draw from 250
Niters = 100; %number of iterations, drawn with replacement from Nsets
binwidth_t = 20; %in microseconds (20 us)
binwidth = binwidth_t*(fs/1e6); %samples

%generate Nsets
Ds = (dur+0.005)*fs; %add 5 ms zero pad (mirroring genANspikes) 
ANsets = cell(ncf_stim,Nsets);
for s = 1:ncf_stim
%     ANcf = zeros(Ds,cfvec_AN,Nfibers,2); %obtain left and right spike train for each draw
    for n = 1:Nsets
        [tstim,flims] = genNBnoise(dur,fs,cfvec_stim(s),bw);
        tstim = rampdamp(tstim,tc,fs);
        tstim(2,:) = tstim; %stereo stimulus to generate L/R pairs of spikes for each stimulus
        psth = genBEZpsth_stochastic(tstim,fs,dB,cf_fib(s,:),Nfibers);
        ANsets{s,n} = squeeze(psth);
    end
end

%obtain centroids
centroid = zeros(ncf_stim,length(dlys),Niters);
for s = 1:ncf_stim
%     ANcf = reshape(ANsets{s,f},Ds,Nfibers*Nsets);
    for n = 1:Niters
        nind = randi(Nsets); %draw single stereo set of spike trains
        for d = 1:length(dlys)
            samp_dly = round(dlys(d)*fs);
            for f = 1:ncf_fib
                xl = cat(1,zeros(samp_dly,Nfibers),squeeze(ANsets{s,nind}(1:end-samp_dly,f,:,1)));
                xr = squeeze(ANsets{s,nind}(:,f,:,2));
                [SCC,xax] = getSCC(xl,xr,binwidth,fs);
                SCCwght(f,:) = centralityWeighting1D(xax,SCC,cf_fib(s,f),'stern','pdf',5);
            end
            SCCall = mean(SCCwght,1);
            centroid(s,d,n) = sum(SCCall.*xax*binwidth_t)/sum(SCCall*binwidth_t); %get centroid
        end
    end
end

%obtain distributions
dist_stim = cell(ncf_stim,length(dlys));
for s = 1:ncf_stim
    for d = 1:length(dlys)
        dist_stim{s,d} = fitdist(squeeze(centroid(s,d,:)),'Normal');
    end
end

%obtain d-prime
sigma_int = 5; %microseonds
for s = 1:ncf_stim
    for d = 2:length(dlys)
        dprime(s,d-1) = (dist_stim{s,d}.mu - dist_stim{s,1}.mu)/...
            sqrt(sigma_int^2+(0.5*(dist_stim{s,d}.sigma^2 + dist_stim{s,1}.sigma^2)));
    end
end

%plot interpolated sensitivity and obtain ITD threshold
DPT = sqrt(2)*norminv(0.71); %2-down/1-up (71% correct)
x = log(abs(dlys(2:end)*1e6));
figure
for f = 1:ncf_stim
    y = abs(dprime(f,:));
    [fity,fitx,y_anch,thresh(f)] = fitInterp(x,y,DPT,0);
%     [fity,fitx,y_anch,thresh(f)] = fitSigmoid(x,y,DPT,0);
    subplot(1,ncf_stim,f)
    plot(fitx,fity,'b','linewidth',1.5) %plot interpolated function
    hold on
    plot(x,y_anch,'ko','markerfacecolor','k','markersize',10) %plot anchor points
    ylim([0 5])
    xline(log(thresh(f)),'r--','linewidth',1.5);
    xticks(x);
    for i = 1:length(x)
        xlbls{1,i} = num2str(exp(x(i)));
    end
    xticklabels(xlbls);
    xlabel('ITD \mus','fontsize',12)
    ylabel('d-prime','fontsize',12)
    title([num2str(cfvec_stim(f)) ' Hz'],'fontsize',14)
end
clear xlbls

%plot thresholds as a function of CF
thresh_p1 = thresh(thresh ~= max(exp(x)));
thresh_p2 = thresh(thresh == max(exp(x))); %identify all thresholds that equal the maximum possible ITD
figure
xf = [1:length(thresh_p1)];
semilogy(xf,thresh_p1,'bs--','markerfacecolor','b','markersize',12,...
    'linewidth',2);
hold on
for i = 1:length(thresh_p2)
    semilogy(length(xf)+i,2000,'bs','markerfacecolor','b','markersize',12,...
        'linewidth',2);
end
xlim([0.5 ncf_stim+0.5])
xticks(1:length(thresh))
for i = 1:length(thresh)
    xlbls{1,i} = num2str(cfvec_stim(i));
end
xticklabels(xlbls);
ylim([10 2000])
yline(1000,'r:','linewidth',2);
yticks([10,100,1000,2000])
yticklabels({'10','100','1000','u.m.'});
xlabel('Center Frequency (Hz)','fontsize',14)
ylabel('ITD threshold (\mus)','fontsize',14)
set(gca,'fontsize',12,'linewidth',1.5)
title('Narrowband Noise')
