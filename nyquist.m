%function  result_t = nyquist(taps, beta)
clear all; close all;
taps = 101;
beta = 0.3;
f_samp = 3e6;
f_symb = 1e6;
T = 1/f_symb;
result = zeros(1, taps);
stepo = 1 / taps * f_samp;
highf = stepo * (taps - 1)/2;
f = linspace(-highf, highf, taps);
lcorner = (1 - beta) / (2 * T);
rcorner = (1 + beta) / (2 * T);
for i = 1:1:taps
	if (abs(f(i)) < lcorner)
		result(i) = sqrt(T);
	elseif (abs(f(i)) > rcorner)
		result(i) = 0;
    else
        result(i) = sqrt(T / 2 * (1 + cos(pi * T / beta * (abs(f(i)) - lcorner))));
    end 
end
figure;
stem(result)
result_t = ifft(result);
result_t = [result_t(52:end) result_t(1:51)];
result_t = result_t./result_t(51);
figure;
t = (-(taps - 1) / 2 : (taps - 1) / 2)./f_samp;
plot(t,result_t)
