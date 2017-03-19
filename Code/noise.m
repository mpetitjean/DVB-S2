function message_noisy = noise(message, ratio_min, step, ratio_max, f_samp, E_b)

% Compute a matrix, each line being the noisy message with different SNR
% 
% INPUTS:
% - message : noiseless binary message
% - ratio_min : lowest SNR
% - ratio_max : highest SNR
% - step : step from ratio min to max
%
% OUTPUTS:
% - matrix of noisy messages

% Various SNRs
ratio = ratio_min:step:ratio_max;
message_noisy = zeros(length(ratio),length(message));

figure;
for i = 1:length(ratio)
	% Compute N0 corresponding to the wanted SNR then add noise
	N0 = E_b/(10^(ratio(i)/10));
	noisePower = N0*f_samp*2;
	message_noisy(i,:) = sqrt(noisePower/2)*(randn(1,length(message)) + 1i*randn(1,length(message))) + message;
    semilogy(i,N0, '-x');
    hold on
end