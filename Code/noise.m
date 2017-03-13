function message_noisy = noise(message, ratio_min, step, ratio_max)

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

% Compute bit energy
E_b = sum(message.^2);

% Compute N0
ratio = ratio_min:step:ratio_max;
for i = 1:length(ratio)
	N0(i) = E_b/(10^(ratio(i)/10));
end

% Add noise




end
