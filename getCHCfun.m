function CHCfun = getCHCfun(cfs,type,pflg)
%Returns surface describing the detection threshold for each pair
%(cihc,cohc). Used to compute detection thresholds corresponding to amounts
%of hair cell damage with certain OHC/IHC proportions
%INPUTS
%cfs = set of frequencies to compute surface for. an_cfs and stim_cfs will
%refer to the same element in cfs
%type = if 1, obtain detection thresholds separately for OHC and IHC loss.
%If 2, obtain detection thresholds for joint OHC/IHC loss
%pflg = plot flag
%OUPUTS
%threshq = cfsxOHCxIHC matrix (each column refers to a given CIHC value 

%created by Luke Baltzell 05/14/21

%set resolutions
CHC = [0:0.1:1];
CHCq = [0:0.01:1];

%obtain type 1 thresholds
if type == 1
    thresh = zeros(2,length(CHC));
    threshq = zeros(length(cfs),2,length(CHCq));
    for f = 1:length(cfs)
        for n = 1:length(CHC)
            dBo = getDetectThresh(cfs(f),cfs(f),[],CHC(n),1);
            dBi = getDetectThresh(cfs(f),cfs(f),[],1,CHC(n));
            thresh(1,n) = SPLtoHL(dBo,cfs(f));
            thresh(2,n) = SPLtoHL(dBi,cfs(f));
        end
        threshq(f,1,:) = interp1(CHC,thresh(1,:),CHCq);
        threshq(f,2,:) = interp1(CHC,thresh(2,:),CHCq);
    end
    
    if pltflg == 1
        figure
        for i = 1:length(cfs)
            subplot(1,length(cfs),i)
            plot(CHC,thresh(1,:),'k-','linewidth',2)
            hold on
            plot(CHC,thresh(1,:),'b--','linewidth',2)
            xlabel('cohc/cihc')
            ylabel('detection threshold (dB HL)')
            legend('OHC loss (cihc = 1)','IHC loss (cohc = 1)')
            title(num2str(cfs(i)))
        end
    end

% obtain type 2 thresholds
elseif type == 2
    thresh = zeros(length(CHC),length(CHC));
    threshq = zeros(length(cfs),length(CHCq),length(CHCq));
    for f = 1:length(cfs)
        for n = 1:length(CHC)
            for m = 1:length(CHC)
                dB = getDetectThresh(cfs(f),cfs(f),[],CHC(n),CHC(m));
                thresh(n,m) = SPLtoHL(dB,cfs(f));
            end
        end
        threshq(f,:,:) = interp2(CHC',CHC,thresh,CHCq',CHCq);
    end
    
    if pflg == 1
        figure
        for i = 1:length(cfs)
            subplot(1,length(cfs),i)
            surf(CHCq,CHCq,squeeze(threshq(i,:,:)))
            view(2)
            xlabel('CIHC')
            ylabel('COHC')
            c = colorbar;
            c.Label.String = 'Threshold (dB HL)';
            title(num2str(cfs(i)))
        end
    end
end

CHCfun.thresh = threshq;
CHCfun.CHC = CHCq;
CHCfun.cfs = cfs;
end