%function  result_t = nyquist(taps, beta)
taps = 51;
beta = 0.5;
f_samp = 8e8;
f_symb = 1e6;
T = 1/f_symb;
result = zeros(1, taps);
stepo = 1 / taps * f_symb;
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
figure;
plot(t,result_t)
t = (-(taps - 1) / 2 : (taps - 1) / 2)./f_samp;

