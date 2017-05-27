function [shiftoa, dftoa] = frame(y, a, N, K, f_sym)
% y: receiver samples (after downsampling) (row vector) 
% a: Pilot (row vector, length N)
% f_sym : symbol frequency
nmax = 200;
D = zeros(K,nmax);
y = conj(y);
for k = 1:K
    for l = k:N-1
        for n = 1:nmax
	 		y1 = y(n+l).*a(l+1); % length(y1) = length(y)
	        y2 = conj(y(n+l-k).*a(l-k+1)); % length(y2) = length(y)
	        D(k,n) = D(k,n) + y1.* y2;
	%         y1 = zeros(size(y));
	%         y2 = y1;
	%         y1(1:end-l) = y(l+1:end);
	%         y2(1:end-(l-k)) = y(l-k+1:end);
	%         D(k+1,:) = D(k+1,:) + conj(y1).*a(l+1).*conj(conj(y2).*a(l-k+1));
        end
    end
    D(k,:) = D(k,:)/(N-k);
end
% D = D(:,1:end-N+1);
[~, shiftoa] = max(sum(abs(D)));
Y = f_sym*angle(D(:,shiftoa))./(2*pi*(1:K).');
dftoa = -sum(Y)/K;
% cfoHat = 0;
% for k=1:K
%         cfoHat = cfoHat + angle(D(k,shiftoa))/(2*pi*k/f_sym);
% end
%     cfoHat = -cfoHat/K;
%     dftoa = cfoHat;
