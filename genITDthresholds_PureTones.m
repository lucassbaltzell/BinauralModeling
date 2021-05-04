%This script simulates ITD sensitivity for a set of pure tones. This
%simulation is based on simulated spike trains at the output of the BEZ
%phenomenolgical auditory nerve (AN) model. This method roughly follows the
%procedure described by Moncada-Torres et al. (2018)
%https://doi.org/10.1121/1.5051322, although we obtain somewhat different 
%results. The results of this simulation can be compared with the pure-tone
%ITD thresholds in human listeners obtained by Brughera et al. (2013)

%In this script, information is not integrated across multiple CFs, and
%simulated ITD thresholds reflect a single set of fibers with CFs
%corresponding to the frequency of the pure tone

%For the sake of efficiency, spike trains are generated using monaural
%stimuli, and ITDs will be applied by delaying the spike trains. The high
%sample rate (100 kHz) allows for a delay resolution of 10 us, and so
%target ITD values (dlys) must be integer multiples of 10 us.

%%%%ONLY FOR PURE TONE SIMULATIONS%%%%
%Also for the sake of efficiency (improved sampling), rather than simulate 
%new spike trains for each iteration over which the distribution is built,
%a set of spike trains are defined (Nsets), and are sampled with
%replacement on each iteration (Niters). Because pure tone stimuli are the 
%same for each set of Nfibers in Nset, a pool of Nfibers*Nsets is defined, 
%and on each iteration, a set of Nfibers is drawn for each ear from the 
%pool. This reduces the number of simulated spike trains required to
%acheive random sampling representative of the population
%%%%_____________________________%%%%

%Delays are limited to 320 us becasue for pure tones, the effective ITD
%decreases if the delay exceeds a duration corresponding to half of the
%period. In practice, this decrease in sensitivity occurs as the delay
%approaches half of the period (e.g. Sayers et al., 1964), and seems to be
%maximal around T/4, so for a 1.5 kHz pure tone, we expect sensitivity to
%decrease at delays greater than 167 microseconds. For this reason, the
%fitInterp and fitSigmoid functions carry the maximum sensitivity obtained
%at any delay through to all larger delays.

%created by Luke Baltzell for presentation at Binaural Bash 2020. Modified 
%by Luke Baltzell on 04/27/21

%define pure tone parameters (following Brughera et al., 2013)
fs = 100000;
dur = 0.5; %500 ms
tc = 0.1; %100 ms ramp
t = [1/fs:1/fs:dur];
dB = 70;
dlys = [0 10 20 40 80 160 320]*1e-6; %set of delays over which to calculate sensitivity
cfvec = [250 500 750 1000 1250 1500];
cf_n = length(cfvec);
% cfvec = logspace(log10(200),log10(10000),cf_n);

%generate pure tones
Stim = cell(1,cf_n);
for f = 1:cf_n
    pt = sin(2*pi*cfvec(f)*t);
    Stim{1,f} = rampdamp(pt,tc,fs);
end

%set AN parameters
Nfibers = 35; %desired number of fibers for each CF
Nsets = 250; %number of sets to draw from
Niters = 100; %number of iterations, drawn with replacement from Nsets
binwidth_t = 20; %in microseconds (20 us)
binwidth = binwidth_t*(fs/1e6); %samples

%generate Nsets
Ds = (dur+0.005)*fs; %add 5 ms zero pad (mirroring genANspikes) 
ANsets = cell(1,cf_n); 
for f = 1:cf_n
    tstim = Stim{1,f};
    ANcf = zeros(Ds,Nfibers,Nsets); %mono stimuli for spike generation
    for n = 1:Nsets
        psth = genBEZpsth_stochastic(tstim,fs,dB,cfvec(f),Nfibers);
        ANcf(:,:,n) = squeeze(psth);
    end
    ANsets{1,f} = ANcf;
end

%obtain centroids
centroid = zeros(cf_n,length(dlys),Niters);
for f = 1:cf_n
    ANcf = reshape(ANsets{1,f},Ds,Nfibers*Nsets);
    for n = 1:Niters
        %becasue stimuli are mono, a different spike train index can be
        %used for each ear
        inds = randi(Nfibers*Nsets,Nfibers,2); %generate random indexes (Nfibers for each ear)
        for d = 1:length(dlys)
            samp_dly = round(dlys(d)*fs);
            xl = cat(1,zeros(samp_dly,Nfibers),ANcf(1:end-samp_dly,inds(:,1)));
            xr = ANcf(:,inds(:,2));
            [SCC,xax] = getSCC(xl,xr,binwidth,fs);
            SCCwght = centralityWeighting1D(xax,SCC,cfvec(f),'stern','pdf',5);
            centroid(f,d,n) = sum(SCCwght.*xax*binwidth_t)/sum(SCCwght*binwidth_t); %get centroid
        end
    end
end

%obtain distributions
dist_pt = cell(cf_n,length(dlys));
for f = 1:cf_n
    for d = 1:length(dlys)
        dist_pt{f,d} = fitdist(squeeze(centroid(f,d,:)),'Normal');
    end
end

%obtain d-prime
sigma_int = 5; %microseonds
for f = 1:cf_n
    for d = 2:length(dlys)
        dprime(f,d-1) = (dist_pt{f,d}.mu - dist_pt{f,1}.mu)/...
            sqrt(sigma_int^2+(0.5*(dist_pt{f,d}.sigma^2 + dist_pt{f,1}.sigma^2)));
    end
end

%plot interpolated sensitivity and obtain ITD threshold
DPT = sqrt(2)*norminv(0.8); %3-down/1-up (80% correct)
x = log(abs(dlys(2:end)*1e6));
figure
for f = 1:cf_n
    y = abs(dprime(f,:));
    [fity,fitx,y_anch,thresh(f)] = fitInterp(x,y,DPT,0);
%     [fity,fitx,y_anch,thresh(f)] = fitSigmoid(x,y,DPT,0);
    subplot(1,cf_n,f)
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
    title([num2str(cfvec(f)) ' Hz'],'fontsize',14)
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
xlim([0.5 cf_n+0.5])
xticks(1:length(thresh))
for i = 1:length(thresh)
    xlbls{1,i} = num2str(cfvec(i));
end
xticklabels(xlbls);
ylim([10 2000])
yline(1000,'r:','linewidth',2);
yticks([10,100,1000,2000])
yticklabels({'10','100','1000','u.m.'});
xlabel('Frequency (Hz)','fontsize',14)
ylabel('ITD threshold (\mus)','fontsize',14)
set(gca,'fontsize',12,'linewidth',1.5)
title('Pure Tones')