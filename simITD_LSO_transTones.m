pth = pwd;

%simulate Monaghan et al. (2015)
fs = 100000;
dur = 0.3;
tc = 0.02;
dB = 75;

dlys = [0 10 20 40 80 160 320 640 1280]*1e-6;
% dlys = [-4000:250:4000]*1e-6;
cfvec = logspace(log10(100),log10(10000),60);
% cf_fib = [3200 4000 4700];


addpath('/Users/lukebaltzell/Documents/UR_EAR_2020b')
addpath('/Users/lukebaltzell/Documents/UR_EAR_2020b/Luke_ITDmodelling')
addpath('/Users/lukebaltzell/Documents/UR_EAR_2020b/Luke_ITDmodelling/ashida_lso')
% addpath('/Users/lukebaltzell/Documents/JASA-04864R2_CODE/IBiDT/model_stage/periphery/zilany_model/zilany_2018')

%generate SAM
p_car = 4000;
p_mods = [32 64 128 256 512 1024];

Rate_Dif = cell(1,length(p_mods));
for mf = 1:length(p_mods)
tone = genTransTone(fs,dur,p_car,p_mods(mf));
tone = rampdamp(tone,tc,fs);
flims = [p_car-p_mods(mf) p_car+p_mods(mf)];
flims_extend = [flims(1)/(2^(1/(2*3))) flims(2)*(2^(1/(2*3)))]; %extend by 3rd oct
%get lower cutoff
[~,tmpind] = find(flims_extend(1) > cfvec);
tmp_cf_lc = cfvec(tmpind(end));
flims_extend_low = flims_extend(1)/(2^(1/(2*8))); %include if within 8th octave
if tmp_cf_lc > flims_extend_low
    cf_lc_ind = tmpind(end);
else
    cf_lc_ind = tmpind(end)+1;
end
%get upper cutoff
[~,tmpind] = find(flims_extend(2) < cfvec);
tmp_cf_uc = cfvec(tmpind(1));
flims_extend_hi = flims_extend(2)*(2^(1/(2*8))); %include if within 8th octave
if tmp_cf_uc < flims_extend_hi
    cf_uc_ind = tmpind(1);
else
    cf_uc_ind = tmpind(1)-1;
end
cf_fib = cfvec(cf_lc_ind:cf_uc_ind);
n_cf = length(cf_fib);

%set LSO parameters
% lso.Tref = 1.6e-3;    %1.6; % [ms] 
lso.Tref = 1.6e-3;    %1.6; % [ms] ... LB mod
lso.ThEx = 3;         %8;   % threshold 
% lso.WiEx = 1.1e-3;    %0.8; % [ms] coincidence window
lso.WiEx = 1.1e-3;    %0.8; % [ms] coincidence window
lso.gIn = 2;          %2;   % threshold increase by inhibition 
% lso.WiIn = 3.1e-3;    %1.6; % [ms] inhibition window
lso.WiIn = 3.1e-3;    %1.6; % [ms] inhibition window... LB mod
fibersPerNeuron_ipsi_ex = 20;
fibersPerNeuron_contra_in = 8;
% fibersPerNeuron_ipsi_ex = 30;
% fibersPerNeuron_contra_in = 12;
Npsth = fibersPerNeuron_ipsi_ex + fibersPerNeuron_contra_in;

%% generate AN spikes
DT = 1/fs;
Ts = dur;

Nset = 200;
[~,D] = genANspikes(randn(1,round(dur*fs)),fs,dB,500,0); %dummy spike gen to get params
binwidth_t = 20; %in us
binwidth = binwidth_t*(fs/1e6); %samples
edges = [-D:binwidth:D];
xax = edges(1:end-1)*(1e6/fs); %us

%generate spikes
ANspikes = zeros(Npsth,D,2,n_cf);
for n = 1:Nset
    for f = 1:n_cf
        %generate AN spikes
        ANspikes(n,:,:,f) = genANspikes(tone,fs,dB,cf_fib(f),1);
    end
end

iters = 200;
for n = 1:iters
    %draw AN spikes
    inds = randperm(Nset,Npsth);
    inds_ex = inds(1:fibersPerNeuron_ipsi_ex);
    inds_in = inds(fibersPerNeuron_ipsi_ex+1:end);
    for d = 1:length(dlys)
        for f = 1:n_cf
            %Apply Delays
            if dlys(d) >= 0
                AN_ipsi_ex_Lhem = spikedelay(squeeze(ANspikes(inds_ex,:,1,f)),dlys(d),fs);
                AN_ipsi_ex_Lhem = sum(AN_ipsi_ex_Lhem,1);
                AN_ipsi_ex_Rhem = sum(squeeze(ANspikes(inds_ex,:,2,f)),1);
            
                AN_contra_in_Lhem = sum(squeeze(ANspikes(inds_in,:,2,f)),1);  
                AN_contra_in_Rhem = spikedelay(squeeze(ANspikes(inds_in,:,1,f)),dlys(d),fs);
                AN_contra_in_Rhem = sum(AN_contra_in_Rhem,1);
            else
                AN_ipsi_ex_Lhem = sum(squeeze(ANspikes(inds_ex,:,1,f)),1);
                AN_ipsi_ex_Rhem = spikedelay(squeeze(ANspikes(inds_ex,:,2,f)),abs(dlys(d)),fs);
                AN_ipsi_ex_Rhem = sum(AN_ipsi_ex_Rhem,1);
                 
                AN_contra_in_Lhem = spikedelay(squeeze(ANspikes(inds_in,:,2,f)),abs(dlys(d)),fs);
                AN_contra_in_Lhem = sum(AN_contra_in_Lhem,1);
                AN_contra_in_Rhem = sum(squeeze(ANspikes(inds_in,:,1,f)),1);
            end
            
            
            %Generate CD
            cdOut_Lhem = LSOmodelCOC(AN_ipsi_ex_Lhem,AN_contra_in_Lhem,lso,DT);
            cdOut_Rhem = LSOmodelCOC(AN_ipsi_ex_Rhem,AN_contra_in_Rhem,lso,DT);
            rate_Lhem(n,d,f) = sum(cdOut_Lhem,2)/dur;
            rate_Rhem(n,d,f) = sum(cdOut_Rhem,2)/dur;
            rate_dif(n,d,f) = rate_Lhem(n,d,f) - rate_Rhem(n,d,f);
        end
    end
end
Rate_Dif{1,mf} = rate_dif;
end

% %average over cf
% rate_dif_mu = mean(rate_dif,3);

%calculate d-prime for each cf
dist = cell(length(p_mods),length(dlys));
for mf = 1:length(p_mods)
    rate_dif_mu = mean(Rate_Dif{1,mf},3);
    for d = 1:length(dlys)
        dist{mf,d} = fitdist(rate_dif_mu(:,d),'Normal');
    end
end

sigma_int = 2; %delta spikes/sec
for mf = 1:length(p_mods)
    for d = 2:length(dlys)
        dprime(mf,d-1) = (dist{mf,d}.mu - dist{mf,1}.mu)/...
            sqrt(sigma_int^2+(0.5*(dist{mf,d}.sigma^2 + dist{mf,1}.sigma^2)));
    end
end

%obtain threshold
DPT = 1.2; %80% correct (Monaghan)
x = log(abs(dlys(2:end)*1e6));
h = figure;
for mf = 1:length(p_mods)
    subplot(1,length(p_mods),mf)
    y = abs(dprime(mf,:));
    %         [~,~,fsig,fsx,yco,itd_thresh(i)] = fitsigmoid(x,y,0,DPT);
    [fsig,fsx,yco,itd_thresh(mf)] = fitInterp(x,y,0,DPT);
    plot(fsx,fsig,'b','linewidth',1.5)
    hold on
    plot(x,yco,'ko','markerfacecolor','k','markersize',10)
    ylim([0 5])
    xline(log(itd_thresh(mf)),'r--','linewidth',1.5);
    xticks(x);
    xlbls = exp(x);
    xticklabels({num2str(xlbls(1)),num2str(xlbls(2)),num2str(xlbls(3)),...
        num2str(xlbls(4)),num2str(xlbls(5)),num2str(xlbls(6)),...
        num2str(xlbls(7)),num2str(xlbls(8))});
    xlabel('ITD \mus')
    ylabel('d-prime')
    title(['4 kHz Transposed Tone (' num2str(p_mods(mf)) ' Hz)'])
    set(gca,'fontsize',14)
end

itd_thresh(itd_thresh == max(xlbls)) = 2000;
Monaghan.thresholds = [36,50,62,2000;...
    47,40,150,2000;...
    110,142,560,2000;...
    185,395,2000,2000;...
    220,270,2000,2000];
Monaghan.MFs = [128 256 512 1024];
BT.thresholds_mu = [105,85,76,105];
BT.MFs = [32,64,128,256];

h = figure;
[~,ind] = find(itd_thresh == 2000);
y = itd_thresh(1:min(ind)-1);
x = p_mods(1:length(y));
h1 = loglog(x,y,'bs--','markerfacecolor','b','markersize',14,...
    'linewidth',2);
hold on
if length(ind) == 1
    loglog(p_mods(end),2000,'bs','markerfacecolor','b','markersize',14,...
        'linewidth',2);
elseif length(ind) == 2
    loglog(p_mods(end-1),2000,'bs','markerfacecolor','b','markersize',14,...
        'linewidth',2);
end

h2 = loglog(BT.MFs,BT.thresholds_mu,'k^:','markerfacecolor','k','markersize',11,...
    'linewidth',1.5);

for i = 1:size(Monaghan.thresholds,1)
    y = Monaghan.thresholds(i,:);
    x = Monaghan.MFs;
    h3 = loglog(x,y,'o','color',[0.5 0.5 0.5],'markerfacecolor',[0.9 0.9 0.9],...
        'markersize',10,'linewidth',1.5);
end

xticks(p_mods)
ylim([30 2000])
yticks([100 1000 2000])
yticklabels({'100','1000','u.m.'})
yline(1000,':','color',[0.5 0 0],'linewidth',2);
xlabel('Modulation Frequency (Hz)')
ylabel('ITD threshold')
legend([h1 h2 h3],{'model predictions','data from Bernstein & Trahiotis','data from Monaghan et al.'},...
    'location','northwest')
title('4 kHz carrier')
set(gca,'fontsize',16,'linewidth',2)

LSO.cf_fib = cf_fib;
LSO.params = lso;
LSO.dur = dur;
LSO.dB = dB;
LSO.sigma_int = sigma_int;
LSO.rate_dif = rate_dif;
LSO.rate_Lhem = rate_Lhem;
LSO.rate_Rhem = rate_Rhem;
LSO.dlys = dlys;
LSO.n_ipsi_ex = fibersPerNeuron_ipsi_ex;
LSO.n_contra_in = fibersPerNeuron_contra_in;
cd('LSO')
savefig(h,'LSO_TransTone_threshold')
save('LSO_TransTone_threshold','LSO')
cd(pth)