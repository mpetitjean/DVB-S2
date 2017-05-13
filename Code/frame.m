function [shiftoa, dftoa] = frame(y, a, N, K, f_sym)

D = zeros(K,length(y));
for k = 1:K
    coef = 1/(N-k);
    for ll = k:N-1
        rightshift = circshift(y,[0 ll+1]);
        rightshift(1:ll+1) = 0;
        rightshift2 = circshift(y,[0 ll-k+1]);
        rightshift2(1:ll-k+1) = 0;
        D(k,:) = D(k,:) + conj(rightshift).*a(ll+1).*conj(conj(rightshift2).*a(ll-k+1));
    end
    D(k,:) = coef * D(k,:);
end
[~, shiftoa] = max(sum(abs(D)));
D = f_sym*angle(D(:,shiftoa))./(2*pi*(1:K).');
dftoa = -sum(D)/K;
