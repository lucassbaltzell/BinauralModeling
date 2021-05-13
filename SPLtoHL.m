function dBout = SPLtoHL(dBin,fin)
%returns dB HL for dB SPL input at frequency f
%INPUTS
%dBin = input dB (SPL)
%fin = frequency of tone (Hz)
%OUTPUTS
%dB_HL = output dB (HL)

%created by Luke Baltzell 05/13/21

%HL curve
fqs = [125 160 200 250 315 400 500 630 750 800 1000 ...
    1250 1500 1600 2000 2500 3000 3150 4000 5000 6000 6300 8000];
dB_SPL = [45 38.5 32.5 27 22 17 13.5 10.5 9 8.5 7.5 ...
    7.5 7.5 8 9 10.5 11.5 11.5 12 11 16 21 15.5];

[remain,ind] = min(abs(fin-fqs));
if remain == 0
    dBout = dBin - dB_SPL(ind);
else
    fqs_interp = [fqs(1):5:fqs(end)];
    dB_SPL_interp = interp1(fqs,dB_SPL,fqs_interp,'pchip');
    [~,ind] = min(abs(fin-fqs_interp));
    dBout = dBin - dB_SPL_interp(ind);
end

end