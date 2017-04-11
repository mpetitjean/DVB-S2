function bloc_o = tanner(bloc_in,H)
%{
    Perform one iteration of hard decoding for the tanner graph
%}
    % Compute error pattern
    error = mod(H*bloc_in',2);
    vote = zeros(size(bloc_in));
    % 
    for i=1:length(error)
        vote = vote + (xor(error(i),bloc_in)).*H(i,:);
    end
    bloc_o = (vote > sum(H,1)/2);
