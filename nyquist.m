function  result_t = nyquist(taps, beta)
T = 1e-6;
result = zeros(taps,1);
stepo = 1 / taps * 8e8;
highf = stepo * (taps - 1)/2;
f = linspace(-highf, highf, taps);
lcorner = (1 - beta) / (2 * T);
rcorner = (1 + beta) / (2 * T);
for i = 1:1:taps
	if (abs(f(i)) < lcorner)
		result(i) = sqrt(T);
	elseif (abs(f(i)) > lcorner)
		result(i) = 0;
	else
		result(i) = sqrt(T / 2 * (1 + cos(pi * T / beta * (abs(f(i)) - lcorner))));
	end
end
stem(result)
result_t = ifftshift(result);
t = (-(taps - 1) / 2 : (taps - 1) / 2) / 8e8;

