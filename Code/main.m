clear; close all
Nbps = 4;
precision = 1e4;
ratio_min = -5;
ratio_max = 15;
code_rate = 1/2;
step = 1;
fm 	= 2e6;
fs	= 200e6;
a = ratio_min:step:ratio_max;
%% 
for maxit=[1 3 5]
    ber = (main_step2(Nbps,precision, ratio_min, step, ratio_max, code_rate, maxit) + main_step2(Nbps,precision, ratio_min, step, ratio_max, code_rate, maxit))./2;
    csvwrite(['ber_Nbps' num2str(Nbps) '(LDPC_it' num2str(maxit) ').csv'],[a.' ber.']);
end

%% 
figure
c=1;
for shift=[0 2 5 10 20]
   ber = (samplingshift(Nbps,precision, ratio_min, step, ratio_max, shift) + samplingshift(Nbps,precision, ratio_min, step, ratio_max, shift))./2;
   semilogy(a, ber)
   leg{c} = ['t0 = ' num2str(shift*fm/fs) ' T'];
   c = c+1;
   hold on
end
legend(leg);