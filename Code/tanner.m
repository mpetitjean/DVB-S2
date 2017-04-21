function bloc_o = tanner(bloc_in,H)
%{
    Perform one iteration of hard decoding for the tanner graph
%}
    % Compute error pattern
    syndrom = mod(H*bloc_in',2);
    bloc_o =  (sum(xor(syndrom,bloc_in).*H) > sum(H)/2);