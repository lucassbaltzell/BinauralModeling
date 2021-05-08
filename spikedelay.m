function spOut = spikedelay(spIn,dly,fs)
%applies delay to spike trains. For a 100 kHz sampling rate, delay step
%size is 10 us
%OUTPUTS
%spOUT = delayed spike trains
%INPUTS
%spIn = input spike trains (psth-by-fiber)
%dly = time delay
%fs = sampling rate

%created by Luke Baltzell, modified 05/08/21

[len_psth,nfibers] = size(spIn);
if nfibers > len_psth
    disp('input should be psth by fiber')
end

samp_dly = round(dly*fs);
spOut = cat(1,zeros(samp_dly,nfibers),squeeze(spIn(1:end-samp_dly,:)));

end