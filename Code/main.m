clear;
close all;

% Simulation parameters
f_sym 	= 500;
f_samp	= 20e3;

Nbps = 6;
modulation = 'qam';
N = 10000*Nbps;

ratio_min = 1;     % Different E_b/N0 values
step = 10;
ratio_max = 551;  

% Message generation
bits = randi(2,1,N)-1;      % random bits generation
%message = 'Hello Cedric';
%bits = text2bits(message);

% Mapping
symb_tx = mapping(bits,Nbps,modulation)';

% Upsampling
message_symb = oversampling(symb_tx, f_sym, f_samp);
% Nyquist

% Add noise
message_noisy = noise(message_symb, ratio_min, step, ratio_max, f_samp);

% Nyquist

% Downsampling
num = size(message_noisy,1);
symb_rx = zeros(num, length(symb_tx));
for i = 1:num
    symb_rx(i,:) = undersampling(message_noisy(i,:), f_sym, f_samp);
end

% Demapping
bits_rx = zeros(num, length(bits));
for i = 1:num
    bits_rx(i,:) = demapping(symb_rx(i,:)',Nbps,modulation);
end

% Compute BER and plot
ber = compute_ber(bits, bits_rx, num);
semilogy(10*log10(ratio_min:step:ratio_max),ber);
xlabel('Ratio $E_b/N_0$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('BER (log scale)', 'Interpreter', 'latex', 'FontSize', 12);
grid on;



    