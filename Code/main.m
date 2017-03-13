clear;
close all;

% Simulation parameters
Nbps = 6;
modulation = 'qam';

% Message generation
message = 'Hello Cedric';
message_bits = text2bits(message);

% Mapping
symb_tx = mapping(message_bits,Nbps,modulation);

% Add noise
message_noisy = noise(message_bits, 1, 10, 151);

