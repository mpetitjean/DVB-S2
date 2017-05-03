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
% Compute N0 corresponding to the wanted SNR then add noise
N0 = E_b./(10.^(ratio/10));
noisePower = N0.'*f_samp*2;
message_noisy = sqrt(noisePower/2).*(randn(length(ratio),length(message)) + 1i*randn(length(ratio),length(message))) + message;

% brawgn()