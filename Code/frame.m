function [shiftoa, dftoa] = frame(y, a, N, K, f_sym)

D = zeros(K,length(y));
for k = 1:K-1
    coef = 1/(N-k);
    for l = k:N-1
        
%         y1 = zeros(size(y));
%         y2 = y1;
%         y1(1:end-l) = y(l+1:end);
%         y2(1:end-(l-k)) = y(l-k+1:end);
%         D(k+1,:) = D(k+1,:) + conj(y1).*a(l+1).*conj(conj(y2).*a(l-k+1));
    end
    D(k+1,:) = coef * D(k+1,:);
end
D = D(:,1:end-N+1);
[~, shiftoa] = max(sum(abs(D)));
D = f_sym*angle(D(:,shiftoa))./(2*pi*(1:K).');
dftoa = -sum(D)/K;
