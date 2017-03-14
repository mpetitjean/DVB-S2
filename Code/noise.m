function message_noisy = noise(message, ratio_min, step, ratio_max, f_samp, Nbps)

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

for i = 1:length(ratio)
	% Compute N0 corresponding to the wanted SNR then add noise
	N0 = 1/(2*length(message)*f_samp)*sum(message.^2)/ratio(i);
	sigma = sqrt(N0*f_samp);
	for k = 1:length(message)
		message_noisy(i,k) = sigma*(randn(1) + 1i*randn(1)) + message(k);
    end
end