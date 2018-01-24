clear;
close all;

i = 1;
min = 0;
max = 100;
FC = 2e9;
step = 25;
for df = [0:25:100].*(1e-6*FC)
    ber(i,:) = (berCFO(4, 1e4, -5, 1, 15, df, 0) + berCFO(4, 1e4, -5, 1, 15, df, 0) + berCFO(4, 1e4, -5, 1, 15, df, 0) + berCFO(4, 1e4, -5, 1, 15, df, 0) )./4;
    disp(i)
    i = i+1;
end

figure;
df = [0:25:100].*(1e-6*FC);
for i = 1:length(df)
    % Je sais Ced c'est degueu
    semilogy(-5:1:15, ber(i,:))
    hold on;
    legendInfo{i} = ['ppm = ' num2str(df(i)/FC*1e6)];
end
legendInfo{i+1} = 'Ideal';
load('ber_th_Nbps4.mat');
plot(ebno16QAM,ber16QAM);
legend(legendInfo);
grid on;
