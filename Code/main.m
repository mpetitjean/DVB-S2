clear;
close all;

% Simulation parameters
f_sym 	= 1e6;
f_samp	= 8e8;

Nbps = 6;
modulation = 'qam';
N = 1000*Nbps;

ratio_min = 1;     % Different E_b/N0 values
step = 10;
ratio_max = 251;  

% Message generation
bits = randi(2,1,N)-1;      % random bits generation

% Mapping
symb_tx = mapping(bits,Nbps,modulation)';
% figure;
% plot(symb_tx, 'x');
% hold on;
E_b  = 1/(2*f_samp*length(symb_tx))*(trapz(abs(symb_tx).^2));

% Upsampling
message_symb = oversampling(symb_tx, f_sym, f_samp);
% Nyquist

% Add noise
message_noisy = noise(message_symb, ratio_min, step, ratio_max, f_samp, E_b);

% Nyquist

% Downsampling
num = size(message_noisy,1);
symb_rx = zeros(num, length(symb_tx));
for i = 1:num
    symb_rx(i,:) = undersampling(message_noisy(i,:), f_sym, f_samp);
end
figure;
plot(symb_rx(2,:), 'x');
grid on;

% Demapping
bits_rx = zeros(num, length(bits));
for i = 1:num
    bits_rx(i,:) = demapping(symb_rx(i,:)',Nbps,modulation);
end

% Compute BER and plot
ber = compute_ber(bits, bits_rx, num);
figure;
semilogy(10*log10(ratio_min:step:ratio_max),ber, '-o');
xlabel('Ratio $E_b/N_0$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('BER (log scale)', 'Interpreter', 'latex', 'FontSize', 12);
grid on;
