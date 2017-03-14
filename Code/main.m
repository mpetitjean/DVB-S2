clear;
close all;

% Simulation parameters
f_sym 	= 500;
f_samp	= 20e3;

Nbps = 6;
modulation = 'qam';
N = 100*Nbps;

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
message_noisy = noise(message_symb, 1, 1, 20, f_samp);

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
compute_ber(bits, bits_rx, num);



    