function ber = compute_ber(bits_tx, bits_rx, num)

N = length(bits_tx);
ber = zeros(1, num);
ber(1,:) = sum(bits_tx ~= bits_rx,2);
ber = ber./N;