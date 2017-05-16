function [nHat, cfoHat] = frameAcq(pilotSymbol, modulated, N, K,FM)
    D = zeros(K, length(modulated));
    for k=1:K-1
        if k==5
            disp('')
        end
        for l = k:N-1
            next = padto((conj(modulated(l+1:end)).*pilotSymbol(l+1)), [length(modulated),1]);
            next = next.*padto(conj(conj(modulated(l-k+1:end)).*pilotSymbol(l-k+1)), [length(modulated),1]);
            next = next.';
            D(k+1,:) = D(k+1,:) + next;
        end
        D(k+1,:) = D(k+1,:)/(N-k);
    end
    D = D(:,1:end-N+1);
    [~,nHat] = max(sum(abs(D)));

Y = FM*angle(D(:,nHat))./(2*pi*(1:K).');
cfoHat = -sum(Y)/K;
%     cfoHat = 0;
%     for k=1:K
%         cfoHat = cfoHat + angle(D(k,nHat))/(2*pi*k/FM);
%     end
%     cfoHat = -cfoHat/K;
end
