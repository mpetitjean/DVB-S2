function bits_info = LDPCDecode(coded_bits, H, maxit) 

% Decoder
sizeC = max(size(H));
sizeF = min(size(H));
ratio = sizeC/sizeF;
l_bits = length(coded_bits);
bits_info = zeros(1, l_bits/ratio);
% take bloc per bloc
for startIndex = 1:sizeC:l_bits
    bloc_in = coded_bits(startIndex:startIndex+sizeC-1);
    iter = 1;
    bloc_o = tanner(bloc_in, H);

    % Iterate on tanner graph while not converged
    while  ~isequal(bloc_in,bloc_o) && (iter < maxit)
        bloc_in = bloc_o;
        iter = iter + 1;
        bloc_o = tanner(bloc_in, H);
    end
    % Drop parity bits
    bits_info(1+(startIndex-1)/ratio:(startIndex-1)/ratio+sizeF) = ...
        bloc_o(end-sizeF+1:end);
end
