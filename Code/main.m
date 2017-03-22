clear all;
close all;

% Simulation parameters
f_sym 	= 1e6;
f_samp	= 4e6;

Nbps = 6;
modulation = 'qam';
N = 100*Nbps;
taps = 101;
rolloff = 0.3;

ratio_min = 1;     % Different E_b/N0 values (dB)
step = 1;
ratio_max = 1;  

% Message generation
bits = randi(2,1,N)-1;      % random bits generation

% Mapping
symb_tx = mapping(bits,Nbps,modulation).';
E_b = 1/(2*f_samp*length(symb_tx))*(trapz(abs(symb_tx).^2));

% Upsampling
message_symb = oversampling(symb_tx, f_sym, f_samp);

% Nyquist
nyquist_impulse = nyquist(taps, rolloff, f_samp, f_sym);
message_symb_n = conv(message_symb, nyquist_impulse);

% Add noise
%message_noisy = noise(message_symb_n, ratio_min, step, ratio_max, f_samp, E_b);
message_noisy = message_symb_n;

% Nyquist
num = size(message_noisy,1);
message_noisy_n = zeros(num, length(message_noisy)+length(nyquist_impulse)-1);
for i = 1:num
	message_noisy_n(i,:) = conv(message_noisy(i,:), fliplr(nyquist_impulse));
end

% Drop meaningless samples
message_noisy_n = message_noisy_n(:,taps:end-(taps-1));

% Downsampling
symb_rx = zeros(num, length(symb_tx));
for i = 1:num
    symb_rx(i,:) = undersampling(message_noisy_n(i,:), f_sym, f_samp);
end
figure;
plot(symb_rx(1,:), 'x');
grid on;

% Demapping
bits_rx = zeros(num, length(bits));
for i = 1:num
    bits_rx(i,:) = demapping(symb_rx(i,:).',Nbps,modulation);
end

% Compute BER and plot
ber = compute_ber(bits, bits_rx, num);
figure;
semilogy(ratio_min:step:ratio_max,ber, '-o');
xlabel('Ratio $E_b/N_0$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('BER (log scale)', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
