function [thresh,cf_fib] = genITDthresholds(stimtype,mdltype,stimpar)
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

%For the sake of efficiency, spike trains are generated using monaural
%stimuli, and ITDs will be applied by delaying the spike trains. The high
%sample rate (100 kHz) allows for a delay resolution of 10 us, and so
%target ITD values (dlys) must be integer multiples of 10 us.

%created by Luke Baltzell for presentation at Binaural Bash 2020. Modified 
%by Luke Baltzell on 05/03/21

% mdltype = 'MSO';
 
if nargin == 2
    stimpar.fs = 100000;
    if strcmp(stimtype,'pureTone') == 1
        %define pure tone paramters following Brughera et al. (2013)
        stimpar.stim_cfs = [250 500 750 1000 1250 1500]; %stimulus center frequencies
        stimpar.dur = 0.5;
        stimpar.tc = 0.1;
        stimpar.dB = 70;
        stimpar.t = [1/stimpar.fs:1/stimpar.fs:stimpar.dur];
        stimpar.adptrck = 3; %3-down/1-up
        stimpar.afc = 2; %2AFC
    elseif strcmp(stimtype,'nbNoise') == 1
        %define narrowband noise paramters following Spencer et al. (2016)
        stimpar.stim_cfs = [500 4000]; %stimulus center frequencies
        stimpar.dur = 0.3;
        stimpar.tc = 0.015;
        stimpar.dB = 65;
        stimpar.bw = 1/3; %3rd octave
        stimpar.tflg = 0;
        stimpar.adptrck = 2; %2-down/1-up
        stimpar.afc = 2; %2AFC
    elseif strcmp(stimtype,'rustleNoise') == 1
        %define rustle noise paramters following Baltzell et al. (in preparation)
        stimpar.stim_cfs = [500 4000]; %stimulus center frequencies
        stimpar.dur = 0.5;
        stimpar.tc = 0.025;
        stimpar.gpwidth = 3620; %gapwidth of rustle in microseconds
        stimpar.dB = 40;
        stimpar.bw = 1; %3rd octave
        stimpar.tflg = 1;
        stimpar.nord = 16;
        stimpar.adptrck = 3; %2-down/1-up
        stimpar.afc = 2; %2AFC
    end
end

ncf_stim = length(stimpar.stim_cfs);
dlys = [0 10 20 40 80 160 320 640 1280]*1e-6; %set of (positive) delays over which to calculate sensitivity
% cfvec_AN = logspace(log10(125),log10(10000),70); %generate dense vector of center frequencies
cfvec_AN = getGreenwoodCF([0:0.01:1]); %generate 100 AN CFs over the length of the cochlea 

%% calculate AN cfs for each stimulus cf
bw_ext = 1/8; %extend bandwidth by 8th octave
eqflg = 1; %get equal numbers of fibers for each stimulus
cf_fib = getCFfib(cfvec_AN,stimpar.stim_cfs,bw_ext,stimtype,stimpar,eqflg);

%% determine binaural type
if strcmp(mdltype,'MSO') == 1
    Nfibers = 28; %desired number of fibers for each CF
    sigma_int = 5; %microseonds
elseif strcmp(mdltype,'LSO') == 1
    lso.Tref = 1.6e-3;
    lso.ThEx = 3;
    lso.WiEx = 1.1e-3;
    lso.gIn = 2;
    lso.WiIn = 3.1e-3;
    fibersPerNeuron_ipsi_ex = 20;
    fibersPerNeuron_contra_in = 8;
    Nfibers = fibersPerNeuron_ipsi_ex + fibersPerNeuron_contra_in;
    sigma_int = 2; %microseonds
else
    error('select model type (MSO or LSO)')
end

%% set AN parameters
Nsets = 100; %number of sets to stimulate
binwidth_t = 20; %in microseconds (20 us)
binwidth = binwidth_t*(stimpar.fs/1e6); %samples
%genBEZpsth_stochastic calls default AN parameters. To change these
%parameters provide an additional argumment ANpar to function call

%% generate psth and obtain laterality estimate
% Ds = (dur+0.005)*stimpar.fs; %add 5 ms zero pad (mirroring genANspikes) 
lat_est = zeros(ncf_stim,length(dlys),Nsets);
for s = 1:ncf_stim
    for n = 1:Nsets
        if strcmp(stimtype,'nbNoise') == 1
            tstim = genNBnoise(stimpar.dur,stimpar.fs,stimpar.stim_cfs(s),...
                stimpar.bw,stimpar.tflg);
            tstim = rampdamp(tstim,stimpar.tc,stimpar.fs);
            tstim(2,:) = tstim; %stereo stimulus to generate L/R pairs of spikes for each stimulus
        elseif strcmp(stimtype,'rustleNoise') == 1
            tstim = genNBrustle(stimpar.dur,stimpar.fs,stimpar.gpwidth,...
                stimpar.stim_cfs(s),stimpar.bw,stimpar.tflg,stimpar.nord);
            tstim = rampdamp(tstim,stimpar.tc,stimpar.fs);
            tstim(2,:) = tstim; %stereo stimulus to generate L/R pairs of spikes for each stimulus
        elseif strcmp(stimtype,'pureTone') == 1
            phi = (2*pi)*rand(1) - pi; %random phase
            tstim = sin(2*pi*stimpar.stim_cfs(s)*stimpar.t + phi);
            tstim = rampdamp(tstim,stimpar.tc,stimpar.fs);
            tstim(2,:) = tstim;
        end
        psth = genBEZpsth_stochastic(tstim,stimpar.fs,stimpar.dB,cf_fib{1,s},Nfibers);
        for d = 1:length(dlys)
            for f = 1:length(cf_fib{1,s})
                if strcmp(mdltype,'MSO') == 1
                    xl = spikedelay(squeeze(psth(:,f,:,1)),dlys(d),stimpar.fs);
                    xr = squeeze(psth(:,f,:,2));
                    %get centrality weighted SCCs
                    [SCC,xax] = getSCC(xl,xr,binwidth,stimpar.fs);
                    SCCwght(f,:) = centralityWeighting1D(xax,SCC,cf_fib{1,s}(f),'stern','pdf',5);
                elseif strcmp(mdltype,'LSO') == 1
                    %get input to left hemifield
                    xl_ipsi_ex = spikedelay(squeeze(psth(:,f,1:fibersPerNeuron_ipsi_ex,1)),dlys(d),stimpar.fs); %apply delay to left channel
                    xl_ipsi_ex = sum(xl_ipsi_ex,2); %sum over ipsi_ex fibers
                    xr_contra_in = sum(squeeze(psth(:,f,fibersPerNeuron_ipsi_ex+1:end,2)),2); %sum over contra_in fibers
                    %get input to right hemifield
                    xr_ipsi_ex = sum(squeeze(psth(:,f,1:fibersPerNeuron_ipsi_ex,2)),2); %sum over ipsi_ex fibers
                    xl_contra_in = spikedelay(squeeze(psth(:,f,fibersPerNeuron_ipsi_ex+1:end,1)),dlys(d),stimpar.fs); %apply delay to left channel
                    xl_contra_in = sum(xl_contra_in,2); %sum over contra in fibers             
                    %gemerate coincidence detection spike trains
                    cdOut_Lhem = LSOmodelCOC(xl_ipsi_ex,xr_contra_in,lso,1/stimpar.fs);
                    cdOut_Rhem = LSOmodelCOC(xr_ipsi_ex,xl_contra_in,lso,1/stimpar.fs);
                    rate_Lhem = sum(cdOut_Lhem,2)/dur;
                    rate_Rhem = sum(cdOut_Rhem,2)/dur;
                    %get spike rate difference across hemifields
                    rate_dif(f) = rate_Lhem - rate_Rhem;
                end
            end
            if strcmp(mdltype,'MSO') == 1
                SCCall = mean(SCCwght,1); %average over cf_fib
                lat_est(s,d,n) = sum(SCCall.*xax*binwidth_t)/sum(SCCall*binwidth_t); %get centroid
            elseif strcmp(mdltype,'LSO') == 1
                lat_est(s,d,n) = mean(rate_dif); %get mean rate difference
            end
        end
    end
end

%% obtain distributions
dist_stim = cell(ncf_stim,length(dlys));
for s = 1:ncf_stim
    for d = 1:length(dlys)
        dist_stim{s,d} = fitdist(squeeze(lat_est(s,d,:)),'Normal');
    end
end

%% obtain d-prime
for s = 1:ncf_stim
    for d = 2:length(dlys)
        dprime(s,d-1) = (dist_stim{s,d}.mu - dist_stim{s,1}.mu)/...
            sqrt(sigma_int^2+(0.5*(dist_stim{s,d}.sigma^2 + dist_stim{s,1}.sigma^2)));
    end
end

%% get d-prime threshold
if stimpar.adptrck == 2
    pc = 0.71; %2-down/1-up
elseif stimpar.adptrck == 3
    pc = 0.78; %3-down/1-up
else
    error('adaptive track not supported, please add yours')
end
%compute d-prime from pc following Green & Dai (1991) P&P
if stimpar.afc == 2
    DPT = sqrt(2)*norminv(pc) + 0; %2AFC
elseif stimpar.afc == 3
    DPT = sqrt(1.7)*norminv(pc) + 0.56; %3AFC
else
    error('mAFC not supported, please add yours')
end

%% plot interpolated sensitivity and obtain ITD threshold
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
    title([num2str(stimpar.stim_cfs(f)) ' Hz'],'fontsize',14)
end
clear xlbls

%% plot thresholds as a function of CF
ind_p1 = find(thresh ~= max(exp(x))); %identify all thresholds below the maximum possible ITD
thresh_p1 = thresh(ind_p1);
ind_p2 = find(thresh == max(exp(x))); %identify all thresholds that equal the maximum possible ITD
thresh_p2 = thresh(ind_p2);

figure
if isempty(thresh_p2)
    semilogy(ind_p1,thresh_p1,'bs--','markerfacecolor','b','markersize',12,...
        'linewidth',2);
else
    semilogy(ind_p1,thresh_p1,'bs','markerfacecolor','b','markersize',12,...
        'linewidth',2);
    hold on
    semilogy(ind_p2,2000,'bs','markerfacecolor','b','markersize',12,...
        'linewidth',2);
end
xlim([0.75 ncf_stim+0.25])
xticks(1:length(thresh))
for i = 1:length(thresh)
    xlbls{1,i} = num2str(stimpar.stim_cfs(i));
end
xticklabels(xlbls);
ylim([10 2000])
yline(1000,'r:','linewidth',2);
yticks([10,100,1000,2000])
yticklabels({'10','100','1000','u.m.'});
xlabel('Center Frequency (Hz)','fontsize',14)
ylabel('ITD threshold (\mus)','fontsize',14)
set(gca,'fontsize',12,'linewidth',1.5)
if strcmp(mdltype,'MSO') == 1
    title('Narrowband Noise (MSO)')
elseif strcmp(mdltype,'LSO') == 1
    title('Narrowband Noise (LSO)')
end

end
