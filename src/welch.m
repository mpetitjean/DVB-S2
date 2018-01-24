function [psd,f] = welch(x,t,L,D)


%%%
%   Evalue la PSD d un signal par la methode de Welch
%%%

N = floor((length(x)-L)/D+1);
Xt = [];

for nn = 1:N,
    Xt = [Xt,x((nn-1)*D+1:(nn-1)*D+L).'];
end;

W = 1-([0:L-1]-(L-1)/2).^2*4/(L+1)^2;
Wt = repmat(W',1,N);
Yf = fftshift(fft(Xt.*Wt,L,1),1);

psd = sum(abs(Yf).^2,2)'/N;
psd = psd/max(psd);
f = [-L/2:L/2-1]/L.*(1/(t(2)-t(1)));

end