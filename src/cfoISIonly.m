function y = cfoISIonly(u, ppm, phi, fs)

% time vector
t = (0:length(u)-1)./fs;

% REMOVE distortion
df = ppm*1e-6*fs;
y = u.*exp(-1j.* (2*pi*df.*t + phi));
