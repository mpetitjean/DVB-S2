clear; close all
Nbps = 2;
precision = 1e4;
ratio_min = -5;
ratio_max = 15;
code_rate = 1/2;
step = 1;
fm 	= 2e6;
fs	= 200e6;
FC = 2e9;
df = [0 10 50].*(1e-6*FC);
phi = 0;
f_sym 	= 1e6;
f_samp	= 16e6;
shift = round(0.45*f_samp/f_sym);
k = 0.1;
Nexp = 25;
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

%% 
figure
hold on
k = [0.02 0.05 0.1];
precision = 1e3;
means = zeros(length(k), precision);
stdv = zeros(length(k), precision);
step = 20;
for ii= 1:length(k)
    for jj = 1:Nexp
        [~, epsilon] = samplingshift(f_sym, f_samp, Nbps,precision, 25, 1, 25, shift, k(ii));
        samplingerror = (shift - epsilon*f_samp/f_sym)/f_samp;
        means(ii,:) = means(ii,:) + samplingerror;
        stdv(ii,:) = stdv(ii,:) + samplingerror.^2;
    end
    means(ii,:) = means(ii,:)/Nexp;
    stdv(ii,:) = sqrt(stdv(ii,:)/Nexp - means(ii,:).^2);
    current = plot(1:step:length(means),means(ii,1:step:end));
    c = get(current, 'color');
    plot(1:step:length(means),means(ii,1:step:end) + stdv(ii,1:step:end),':','Color', c)
    plot(1:step:length(means),means(ii,1:step:end) - stdv(ii,1:step:end),':','Color',c)
end

%% Degueu apd ici

figure
hold on
precision = 1e3;
means = zeros(length(df), precision);
stdv = zeros(length(df), precision);
step = 20;
for ii= 1:length(df)
    for jj = 1:Nexp
        [~, epsilon] = gardnerCFO(f_sym, f_samp, Nbps,precision, 25, 1, 25, shift, k,df(ii),phi);
        samplingerror = (shift - epsilon*f_samp/f_sym)/f_samp;
        means(ii,:) = means(ii,:) + samplingerror;
        stdv(ii,:) = stdv(ii,:) + samplingerror.^2;
    end
    means(ii,:) = means(ii,:)/Nexp;
    stdv(ii,:) = sqrt(stdv(ii,:)/Nexp - means(ii,:).^2);
    current = plot(1:step:length(means),means(ii,1:step:end));
    c = get(current, 'color');
    plot(1:step:length(means),means(ii,1:step:end) + stdv(ii,1:step:end),':','Color', c)
    plot(1:step:length(means),means(ii,1:step:end) - stdv(ii,1:step:end),':','Color',c)
end

%%

figure
hold on
precision = 1e3;
means = zeros(length(df), precision);
stdv = zeros(length(df), precision);
step = 20;
Nw = 40;
Kw = 8;
% for ii= 1:length(df)
%     for jj = 1:Nexp
        epsilon = main_step4(f_sym, f_samp, Nbps,precision, 5, 1, 25,shift, k, df(1), phi,Nw,Kw);
%         samplingerror = (shift - epsilon*f_samp/f_sym)/f_samp;
%         means(ii,:) = means(ii,:) + samplingerror;
%         stdv(ii,:) = stdv(ii,:) + samplingerror.^2;
%     end
%     means(ii,:) = means(ii,:)/Nexp;
%     stdv(ii,:) = sqrt(stdv(ii,:)/Nexp - means(ii,:).^2);
%     current = plot(1:step:length(means),means(ii,1:step:end));
%     c = get(current, 'color');
%     plot(1:step:length(means),means(ii,1:step:end) + stdv(ii,1:step:end),':','Color', c)
%     plot(1:step:length(means),means(ii,1:step:end) - stdv(ii,1:step:end),':','Color',c)
% end

