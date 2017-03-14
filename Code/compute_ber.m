function compute_ber(bits_tx, bits_rx, num)

N = length(bits_tx);
for i = 1:num
    for k = 1:N
        if ( bits_tx(k) ~= bits_rx(i,k));
            ber(i) = ber(i) + 1;
        end
    end
end