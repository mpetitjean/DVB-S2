function bits = text2bits(text_tx)

%%%
% transforms a string in binary sequence
%%%

text_tx_int = uint8(double(text_tx)) ; % 8 bits/character
text_tx_bin = de2bi(text_tx_int,8) ;   % decimal to binary
nbchar = size(text_tx_bin,1);
bits = double(reshape(text_tx_bin’,1,nbchar*8));
