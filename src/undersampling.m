function m = undersampling(m_sampled, f_m, f_s)
% INPUTS :
% - m_base : base message
% - f_m : frequency of the base message
% - f_s : sampling frequency
% OUTPUT :
% - m : oversampled message
% Undersamples the message m_base from frenquency f_s to frenquency f_m

N = length(m_sampled);
r = f_s/f_m; % oversampling factor
m = zeros(1,N/r);
for i=1:N/r
    m(i) = sum(m_sampled((i-1)*r+1:i*r))/r;
end