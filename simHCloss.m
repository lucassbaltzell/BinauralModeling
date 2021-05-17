function [Out,In] = simHCloss(cfs,HIthresh,cohcIn,loadflg)
%Computes OHC/IHC parameters corresponding to detection thresholds in
%HIthresh. This function takes outer hair cell paramters in cohcIn as
%a starting point, and finds cich necessary to acheive desired HIthresh
%INPUTS
%cfs = audiometric frequencies to be evaluated
%HIthresh = detection thresholds in dB HL to be evaluated
%cohcIn = Outer hair cell parameters to be evaluated. By default, the input
%[0 1] is used, which finds corresponding cihc parameters for full outer
%hair cell loss (cohc = 0) and full outer hair cell health (cohc = 1)
%loadflg = 1 if CHCfun has already been generated, 0 if not
%OUTPUTs
%Out = output parameters (cfs-by-cohcIn-by-HIthresh).
%Out.cohc = output cohc paramters required for loss in HIthresh (may be
%different form input)
%Out.cihc = cihc parameters required to acheive HIthresh given desired cohc
%paramters
%Out.thresh = thresholds corresponding to cohc/cihc paramters

%created by Luke Baltzell 05/17/21

if nargin == 2
    cohcIn = [0 1];
    loadflg = 1;
end

if loadflg == 1
    load('CHCfun');
else
    CHCfun = getCHCfun(cfs,2,1);
end

%remove non-monotonicites (smooth)
for c = 1:size(CHCfun.thresh,1)
    for r = 1:size(CHCfun.thresh,2)
        [rv,ind] = min(CHCfun.thresh(c,r,:));
        if ind ~= size(CHCfun.thresh,3)
            CHCfun.thresh(c,r,ind+1:end) = rv;
        end
    end
end

tol = 5; %tolerance in dB
Out.cihc = zeros(length(cfs),length(cohcIn),length(HIthresh));
Out.cohc = zeros(length(cfs),length(cohcIn),length(HIthresh));
Out.thresh = zeros(length(cfs),length(cohcIn),length(HIthresh));
for f = 1:length(cfs)
    cfind = find(cfs(f) == CHCfun.cfs);
    if isempty(cfind)
        error('CHCfun does not contain desired cf')
    end
    thresh = squeeze(CHCfun.thresh(cfind,:,:));
    for c = 1:length(cohcIn)
        [~,ohc_row] = sort(abs(cohcIn(c) - CHCfun.CHC)); %find row corresponding to desired OHC loss
        for lv = 1:length(HIthresh)
            flg = 0;
            fi = 0;
            while flg == 0
                fi = fi+1;
                threshOHC = thresh(ohc_row(fi),:);
                [epsilon,IHCind] = min(abs(threshOHC - HIthresh(lv)));
                IHCind = find(threshOHC == threshOHC(IHCind),1,'last');
                if epsilon <= tol
                    Out.cihc(f,c,lv) = CHCfun.CHC(IHCind);
                    Out.cohc(f,c,lv) = CHCfun.CHC(ohc_row(fi));
                    Out.thresh(f,c,lv) = threshOHC(IHCind);
                    flg = 1;
                end
            end
        end
    end
end

In.thresh = HIthresh;
In.cfs = cfs;
In.cohc = cohcIn;

if loadflg == 0
    save('CHCfun','CHCfun')
end
end
    
    