function y = genSAM(fs,dur,cf,mod,tc,phi)
% This function returns a sinusoidally amplitude-modulated tone
%OUPUTS
%y: SAM tone
%INPUTS
%fs: sampling rate
%dur: duration of SAM tone
%cf: carrier frequency
%mod: modulation frequency
%tc: duration of ramp/damp
%phi: desired starting phase of sine waves. The first element should refer
%to the phase of the carrier, and the second element to the modulator

%created by Luke Baltzell, modified 04/30/21

if nargin == 4
    tc = 0;
    phi = [0 0];
elseif nargin == 5
    phi = [0 0];
end

t = [1/fs:1/fs:dur];
car = sin(2*pi*cf*t + phi(1));
mod = 0.5+0.5*sin(2*pi*mod*t + phi(2));

amtone = car.*mod;
if tc ~= 0
    amtone = rampdamp(amtone,tc,fs);
end
y = amtone./rms(amtone);

end