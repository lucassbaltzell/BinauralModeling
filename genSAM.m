function y = genSAM(fs,dur,cf,mod)

t = [1/fs:1/fs:dur];

car = sin(2*pi*cf*t);
mod = 0.5+0.5*sin(2*pi*mod*t);

amt = car.*mod;
y = amt./rms(amt);

end