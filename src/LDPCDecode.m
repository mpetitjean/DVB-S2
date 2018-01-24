function bits_info = LDPCDecode(coded_bits, H, maxit) 

% Decoder
sizeC = max(size(H));
sizeF = min(size(H));
ratio = sizeC/sizeF;
l_bits = length(coded_bits);
bits_info = zeros(1, l_bits/ratio);
% take bloc per bloc
for startIndex = 1:sizeC:l_bits
    y = coded_bits(startIndex:startIndex+sizeC-1);
    m = y;
    iter = 0;
    syndrom = mod(H*m',2);
    % Iterate on tanner graph while not converged
    while  any(syndrom) && (iter < maxit)
        iter = iter + 1;
        m = tanner(y, syndrom, m, H);
        syndrom = mod(H*m',2);
    end
    % Drop parity bits
    bits_info(1+(startIndex-1)/ratio:(startIndex-1)/ratio+sizeF) = ...
        m(end-sizeF+1:end);
end
