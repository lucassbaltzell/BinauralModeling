function [SCC,x_us] = getSCC(psth1,psth2,binwidth,fs)
%computes the SCC in one direction. Specifically, it computes the 
%difference (in samples) between each spike in psth1 and all spikes in 
%psth2. If the spikes in psth2 are delayed relative to the spikes in psth1
%(i.e. left-leading), the SCC should be shifted in the negative direction.
%See Louage et al. (2004) for details, and please cite if this code is used
%for publication: http://www.physiology.org/doi/10.1152/jn.00816.2003

%SCC: shuffled cross-correlogram for inputs psth1 and psth2
%x_us: x-axis for plotting SCC in microseconds
%psth1: reference psth (post-stimulus time histogram). Each AN fiber should
%ocupy a column (time-by-fibers)
%psth2: psth referenced to psth1
%binwidth: desired width of bin in samples over which coincidence is detected
%fs: sampling rate of psth

if nargin == 3 && nargout == 2
    error('sampling rate required to return x_us')
end

[D,N1] = size(psth1);
[D2,N2] = size(psth2);

%make sure psths are same length
if D-D2 ~= 0
    error('psth1 and psth2 need to have the same duration')
end

r1 = mean(mean(psth1));
r2 = mean(mean(psth2));
SCC = [];
temp = [];
for n = 1:N1
    for nn = 1:N2
        st1 = find(psth1(:,n) == 1);
        st2 = find(psth2(:,nn) == 1);
        for i = 1:length(st1)
            temp(i,:) = st1(i)-st2;
        end
        temprs = reshape(temp,length(st1)*length(st2),1);
        temp = [];
        SCC = cat(1,SCC,temprs);
    end
end
%SCC = cat(1,SCC,-SCC);
norm = N1*N2*r1*r2*binwidth*D; %removed "*2" to account for asymmetric SCC
edges = [-D:binwidth:D]; %get left-sided edges for histogram
histvec = histcounts(SCC,edges);
SCC = histvec./norm;
x_us = edges(1:end-1)*(1e6/fs); %us
end