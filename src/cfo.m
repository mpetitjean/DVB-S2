function y = cfo(u, df, phi, fs)

% time vector
t = (0:length(u)-1)./fs;

% Add distortion
temp = exp(1j.* (2*pi*df .* t + phi));
y = u .* temp;