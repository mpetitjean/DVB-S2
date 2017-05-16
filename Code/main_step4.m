% clear;
% close all;
function [shiftoa, dftoa] = main_step4(f_sym, f_samp, Nbps,precision, ratio_min, step, ratio_max,shift, k, df, phi,Nw,Kw)
% Simulation parameters
f_gardner = f_samp/2;
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
disp('Bits generated')

% Mapping
if(Nbps > 1)
    symb_tx = mapping(bits,Nbps,'qam').';
%     mode = 'qam';
else
    symb_tx = mapping(bits,Nbps,'pam').';
%     mode = 'qam';
end
disp('Mapping done')
% 
% figure;
% plot(symb_tx(1,:), 'x');
% title('Tx');
% grid on;

% Upsampling
message_symb = upsample(symb_tx, f_samp/f_sym);
Pilot = symb_tx(1:Nw/Nbps);
disp('Upsample done')
clear symb_tx

% Nyquist
nyquist_impulse = nyquist(taps, rolloff, f_samp, f_sym);
message_symb_n = conv(message_symb, nyquist_impulse);
disp('First RRC filter done')
clear message_symb

E_b = 1/(2*f_samp*N)*(trapz(abs(message_symb_n).^2));

% Add noise
message_noisy = noise(message_symb_n, ratio_min, step, ratio_max, f_samp, E_b);
disp('Noise added')
clear message_symb_n

% Add CFO
num = size(message_noisy,1);
parfor i = 1:num
    message_noisy(i,:) = cfo(message_noisy(i,:), df, phi, f_samp);
end

% Nyquist
message_noisy_n = zeros(num, length(message_noisy)+length(nyquist_impulse)-1);
normalization = max(real(conv(nyquist_impulse, nyquist_impulse)));

parfor i = 1:num
	message_noisy_n(i,:) = conv(message_noisy(i,:), fliplr(nyquist_impulse))./normalization;
end
disp('Second RRC filter done')
clear message_noisy normalization nyquist_impulse

% Drop meaningless samples
message_noisy_n = message_noisy_n(:,taps:end-(taps-1));

% Downsampling
%symb_rx= downsample(message_noisy_n(:,1+shift:end).', f_samp/f_sym).';
symb_rx= downsample(message_noisy_n(:,1+shift:end).', f_samp/f_gardner).';

% Gardner
corrected = zeros(num,precision);
epsilon = zeros(num,precision);
for ii=1:num
    [corrected(ii,:), epsilon(ii,:)] = gardner(symb_rx(ii,:),k,f_gardner/f_sym);
end
disp('Downsampling done')
clear message_noisy_n
shiftoa = zeros(1,num);
dftoa = shiftoa;
for ii = 1:num
[shiftoa(ii), dftoa(ii)] = frameAcq(Pilot,corrected(ii,:).', Nw/Nbps, Kw, f_sym);
end
% for ii=1:num
%     fd{ii} = corrected(ii, shiftoa(ii):end);
% end
%corrected = corrected(:,
%Manual perfect cfo correction after filtering
% parfor i= 1:num
%     dd{i} = cfo(fdi{i}, dftoa(i), -phi, f_samp).*exp(-1j*(taps-1)/2/f_samp*2*pi*df);
% end
% % figure;
% % plot(symb_rx(1,:), 'x');
% % title('Rx');
% % grid on;
% 
% % Demapping
% %bits_rx = zeros(num,length(bits)-Nw);
% if strcmp(mode,'pam')
%     dd = cellfun(@real,dd,'un',0);
% end
% dd = cellfun(@transpose,dd,'un',0);
% parfor i = 1:num
%         bits_rx{i} = demapping(dd(:,i),Nbps,mode);
% end
% disp('Demapping done')
% clear symb_rx corrected

% Compute BER and plot
%ber = compute_ber(bits, bits_rx, num);
% figure;
% semilogy(ratio_min:step:ratio_max,ber, 'o');
% xlabel('Ratio $E_b/N_0$', 'Interpreter', 'latex', 'FontSize', 12);
% ylabel('BER (log scale)', 'Interpreter', 'latex', 'FontSize', 12);
% grid on;
% hold on;
% load('ber_th_Nbps2.mat');
% semilogy(ebno4QAM,ber4QAM, '-');
% legend('Simulation', 'Theory')