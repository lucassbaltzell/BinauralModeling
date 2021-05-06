function spOut = spikedelay(spIn,dly,fs)

%can take matrix of spikes in the form fiber-by-psth
if size(spIn,1) > size(spIn,2)
    spIn = spIn';
end

spOut = cat(2,zeros(size(spIn,1),round(dly*fs)),spIn);
spOut = spOut(:,1:size(spIn,2));

end