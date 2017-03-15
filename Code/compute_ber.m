function ber = compute_ber(bits_tx, bits_rx, num)

N = length(bits_tx);
ber = zeros(1, num);
for i = 1:num
    for k = 1:N
        if ( bits_tx(k) ~= bits_rx(i,k))
            ber(1,i) = ber(1,i) + 1;
            if (i ~= 1)
                disp('ça bug plus')
            end
        end
    end
end

%ber = ber./N;