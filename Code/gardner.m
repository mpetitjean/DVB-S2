function [corrected, epsilon] = gardner(samples,k,ratio)
epsilon = zeros(1,ceil(length(samples)/ratio));
corrected = zeros(size(epsilon));
corrected(1) = samples(1);
for n=1:length(epsilon)-1
    interpolate = interp1(1:ratio+1,samples(ratio*(n-1)+1:ratio*n+1),[ratio/2+1 ratio+1]-epsilon(n),'pchip');
    corrected(n+1) = interpolate(2);
    epsilon(n+1) = epsilon(n) + 2*k*real(interpolate(1)*(conj(corrected(n+1)) - conj(corrected(n))));
end
epsilon = epsilon./ratio;