function corrected = gardner(samples,k)
epsilon = zeros(1,length(samples)/2);
corrected = zeros(size(epsilon));
corrected(1) = samples(1);
for n=2:length(epsilon)-1
    interpolate = interp1(1:3,samples(2*(n-1)+1:2*n+1),[1.5 3]-epsilon(n));
    corrected(n+1) = interpolate(2);
    epsilon(n+1) = epsilon(n) * 2*k*real(interpolate(1)*(conj(corrected(n))-conj(corrected(n-1))));
end