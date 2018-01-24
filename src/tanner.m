function bloc_o = tanner(y, syndrom, m,H)
%{
    Perform one iteration of hard decoding for the tanner graph
%}
    % Compute error pattern
bloc_o = (sum(xor(syndrom,m).*H) + y) > sum(H)/2;