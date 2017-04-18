clear;
close all;

% Simulation parameters
f_sym 	= 2e6;
f_samp	= 8e6;
Nbps = 6;
taps = 101;
rolloff = 0.3;
ratio_min = 5;     % Different E_b/N0 values (dB)
step = 1;
ratio_max = 5;
maxit = 2;
info_blksize = 128;
code_rate = 1/2;
code_blksize = 128/code_rate;
N = info_blksize*Nbps*code_rate*100;

% Message generation
info_bits = randi(2,1,N)-1;      % random bits generation
disp('Bits generated')


% LDPC coding
% Create initial parity check matrix
H0 = makeLdpc(info_blksize,code_blksize,0,1,3); 

% Encode information bits
[check_bits, H] = makeParityChk(info_bits(1:info_blksize)', H0, 0); 
coded_bits = [check_bits' info_bits(1:info_blksize)];
for start = info_blksize+1:info_blksize:N 
 [check_bits, ~] = makeParityChk(info_bits(start:start+info_blksize-1)', H0, 0);
 coded_bits = [coded_bits check_bits' info_bits(start:start+info_blksize-1)];
end
disp('LDPC encoding done')

% Mapping
if(Nbps > 1)
    symb_tx = mapping(coded_bits,Nbps,'qam').';
else
    symb_tx = mapping(coded_bits,Nbps,'pam').';
end
disp('Mapping done')

% figure;
% plot(symb_tx(1,:), 'x');
% title('Tx');
% grid on;

% Upsampling
message_symb = upsample(symb_tx, f_samp/f_sym);
disp('Upsample done')

% Nyquist
nyquist_impulse = nyquist(taps, rolloff, f_samp, f_sym);
message_symb_n = conv(message_symb, nyquist_impulse);
disp('First RRC filter done')


E_b = 1/(2*f_samp*N)*(trapz(abs(message_symb_n).^2));
% Add noise
message_noisy = noise(message_symb_n, ratio_min, step, ratio_max, f_samp, E_b);
%message_noisy = message_symb_n;
disp('Noise added')

% Nyquist
num = size(message_noisy,1);
message_noisy_n = zeros(num, length(message_noisy)+length(nyquist_impulse)-1);
normalization = max(real(conv(nyquist_impulse, nyquist_impulse)));
parfor i = 1:num
	message_noisy_n(i,:) = conv(message_noisy(i,:), fliplr(nyquist_impulse))./normalization;
end
disp('Second RRC filter done')

% Drop meaningless samples
message_noisy_n = message_noisy_n(:,taps:end-(taps-1));

% Downsampling
symb_rx = zeros(num, length(symb_tx));
parfor i = 1:num
    %symb_rx(i,:) = undersampling(message_noisy_n(i,:), f_sym, f_samp);
    symb_rx(i,:) = downsample(message_noisy_n(i,:), f_samp/f_sym);
end
disp('Downsampling done')
% figure;
% plot(symb_rx(1,:), 'x');
% title('Rx');
% grid on;

% Demapping
bits_rx = zeros(num, length(coded_bits));
parfor i = 1:num
    if(Nbps > 1)
        bits_rx(i,:) = demapping(symb_rx(i,:).',Nbps,'qam');
    else
        bits_rx(i,:) = demapping(real(symb_rx(i,:)).',Nbps,'pam');
    end
end
disp('Demapping done')

% LDPC decoding
bits_info = LDPCDecode(bits_rx, H, maxit);
disp('LDPC decoding done')

% Compute BER and plot
ber = compute_ber(info_bits, bits_info, num);
disp('End')
%figure;
%semilogy(ratio_min:step:ratio_max,ber, '-o');
% xlabel('Ratio $E_b/N_0$', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('BER (log scale)', 'Interpreter', 'latex', 'FontSize', 12);
% grid on;
