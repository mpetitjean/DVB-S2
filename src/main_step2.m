% clear;
% close all;
function ber = main_step2(Nbps,precision, ratio_min, step, ratio_max, code_rate, maxit)
% Simulation parameters
f_sym 	= 2e6;
f_samp	= 8e6;
% Nbps = 6;
taps = 101;
rolloff = 0.3;
% ratio_min = -5;     % Different E_b/N0 values (dB)
% step = 1;
% ratio_max = 15;
% maxit = 2;
info_blksize = 128;
% code_rate = 1/2;
code_blksize = info_blksize/code_rate;
N = info_blksize*Nbps*precision;


% Message generation
info_bits = logical(randi([0 1],[1 N]));      % random bits generation
disp('Bits generated')


% LDPC coding
% Create initial parity check matrix

H = makeLdpc(info_blksize,code_blksize,0,1,3);
info_bits = reshape(info_bits, info_blksize,N/info_blksize);
% Encode information bits
[check_bits, H] = makeParityChk(info_bits, H, 1);

coded_bits = vertcat(check_bits, info_bits);
clear check_bits;
coded_bits = coded_bits(:).';
info_bits = info_bits(:).';
disp('LDPC encoding done')
H = logical(H);
%info_bits = logical(info_bits);
% Mapping
if(Nbps > 1)
    symb_tx = mapping(coded_bits,Nbps,'qam').';
    mode = 'qam';
else
    symb_tx = mapping(coded_bits,Nbps,'pam').';
    mode = 'pam';
end
disp('Mapping done')
clear coded_bits;
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
clear message_symb symb_tx;

E_b = 1/(2*f_samp*N/code_rate)*(trapz(abs(message_symb_n).^2));
% Add noise
message_noisy = noise(message_symb_n, ratio_min, step, ratio_max, f_samp, E_b);
%message_noisy = message_symb_n;
disp('Noise added')
clear message_symb_n;
% Nyquist
num = size(message_noisy,1);
normalization = max(real(conv(nyquist_impulse, nyquist_impulse)));
message_noisy_n = zeros(num,length(message_noisy)+length(nyquist_impulse)-1);

parfor i=1:num
message_noisy_n(i,:) = conv(message_noisy(i,:), fliplr(nyquist_impulse))./normalization;
end

disp('Second RRC filter done')
clear message_noisy nyquist_impulse;
% Drop meaningless samples
message_noisy_n = message_noisy_n(:,taps:end-(taps-1));

% Downsampling
symb_rx= downsample(message_noisy_n.', f_samp/f_sym).';
disp('Downsampling done')
clear message_noisy_n symb_tx;
% figure;
% plot(symb_rx(1,:), 'x');
% title('Rx');
% grid on;

% Demapping
bits_rx = zeros(size(symb_rx,1),size(symb_rx,2)*Nbps);
if strcmp(mode,'pam')
    symb_rx = real(symb_rx);
end

symb_rx = symb_rx.';
parfor i = 1:num
        bits_rx(i,:) = demapping(symb_rx(:,i),Nbps,mode);
end
disp('Demapping done')
clear symb_rx;
% LDPC decoding
bits_info = zeros(size(info_bits));
H = double(H);
parfor i = 1:num
bits_info(i,:) = LDPCDecode(bits_rx(i,:), H, maxit);
end
disp('LDPC decoding done')
clear bits_rx H;
% Compute BER and plot
ber = compute_ber(info_bits, bits_info, num);
disp('End')
% figure;
% semilogy(ratio_min:step:ratio_max,ber, '-o');
% xlabel('Ratio $E_b/N_0$', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('BER (log scale)', 'Interpreter', 'latex', 'FontSize', 12);
% grid on;
