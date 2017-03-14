function m = oversampling(m_base, f_m, f_s)
% INPUTS :
% - m_base : base message
% - f_m : frequency of the base message
% - f_s : sampling frequency
% OUTPUT :
% - m : oversampled message
% Oversamples the message m_base from frenquency f_m to frenquency f_s

N = length(m_base);
r = f_s/f_m; % oversampling factor

m = zeros(1,N*r);
for i=1:N
    m((i-1)*r+1:i*r) = m_base(i); % copies r times the ith value of m_bases
end