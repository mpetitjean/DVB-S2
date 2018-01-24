function [coded_bits, H] = smallLDPCEncode(info_bits, H0)
%{
    Does a 1/2 code rate LDPC encoder on the basis of a hard-coded 
    parity check matrix that must first be arranged in order to make the
    identy matrix appear.
%}
 
Isize = min(size(H0));

if ( mod(length(info_bits),Isize) ~= 0)
    error('Bad number of bits to use H0');
end

% Arrange the matrix so that H = [I P']
combili = mod(inv(H0(:,1:Isize)),2);
Ptransp = mod(combili*H0(:,Isize+1:end),2);
P = Ptransp';
Itest = mod(combili*H0(:,1:Isize),2);
I = eye(Isize);

if ( Itest ~= I)
    error('Could not arrange H0')
end

H = [I Ptransp];
G = [P I];

% Code the information: the encoder simply consists in dividing the bit 
% stream into blocks of bits and by modulo-2 multiplying each block with 
% the generator matrix

ratio = max(size(H0))/Isize;
coded_bits = zeros(1, ratio*length(info_bits));
for blkbegin = 1:Isize:length(info_bits)
    coded_bits(ratio*(blkbegin-1)+1:ratio*(blkbegin-1)+ratio*Isize) = ...
        mod(info_bits(blkbegin:blkbegin+Isize-1)*G,2);
end










