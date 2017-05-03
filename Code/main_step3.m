% clear;
% close all;
function ber = main_step1(Nbps,precision, ratio_min, step, ratio_max)
% Simulation parameters
f_sym 	= 2e6;
f_samp	= 8e6;

%Nbps = 2;
%precision = 1e6;
N = precision*Nbps;
taps = 101;
rolloff = 0.3;

%ratio_min = -5;     % Different E_b/N0 values (dB)
%step = 1;
%ratio_max = 15;

% Message generation
bits = randi([0 1],[1 N]);     % random bits generation
 
% Mapping
if(Nbps > 1)
    symb_tx = mapping(bits,Nbps,'qam').';
else
    symb_tx = mapping(bits,Nbps,'pam').';
end
    
% 
% figure;
% plot(symb_tx(1,:), 'x');
% title('Tx');
% grid on;

% Upsampling
message_symb = upsample(symb_tx, f_samp/f_sym);

% Nyquist
nyquist_impulse = nyquist(taps, rolloff, f_samp, f_sym);
message_symb_n = conv(message_symb, nyquist_impulse);

E_b = 1/(2*f_samp*N)*(trapz(abs(message_symb_n).^2));

% Add noise
message_noisy = noise(message_symb_n, ratio_min, step, ratio_max, f_samp, E_b);
%message_noisy = message_symb_n;

% Nyquist
num = size(message_noisy,1);
message_noisy_n = zeros(num, length(message_noisy)+length(nyquist_impulse)-1);
normalization = max(real(conv(nyquist_impulse, nyquist_impulse)));
parfor i = 1:num
	message_noisy_n(i,:) = conv(message_noisy(i,:), fliplr(nyquist_impulse))./normalization;
end

% Drop meaningless samples
message_noisy_n = message_noisy_n(:,taps:end-(taps-1));

% Downsampling
symb_rx = zeros(num, length(symb_tx));
parfor i = 1:num
    %symb_rx(i,:) = undersampling(message_noisy_n(i,:), f_sym, f_samp);
    symb_rx(i,:) = downsample(message_noisy_n(i,:), f_samp/f_sym);
end

% figure;
% plot(symb_rx(1,:), 'x');
% title('Rx');
% grid on;

% Demapping
bits_rx = zeros(num, length(bits));
parfor i = 1:num
    if(Nbps > 1)
        bits_rx(i,:) = demapping(symb_rx(i,:).',Nbps,'qam');
    else
        bits_rx(i,:) = demapping(real(symb_rx(i,:)).',Nbps,'pam');
    end
end

% Compute BER and plot
ber = compute_ber(bits, bits_rx, num);
% figure;
% semilogy(ratio_min:step:ratio_max,ber, 'o');
% xlabel('Ratio $E_b/N_0$', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('BER (log scale)', 'Interpreter', 'latex', 'FontSize', 12);
% grid on;
% hold on;
% load('ber_th_Nbps2.mat');
% semilogy(ebno4QAM,ber4QAM, '-');
% legend('Simulation', 'Theory')