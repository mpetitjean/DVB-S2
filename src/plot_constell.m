clear;
close all;

load('symb_Nbps2.mat');
figure;
subplot(1,2,1)
plot(symb_rx, 'x')
hold on
plot(symb_tx, '*')
title('Nbps = 2')

load('symb_Nbps6.mat');
subplot(1,2,2)
plot(symb_rx, 'x')
hold on
plot(symb_tx, '*')
title('Nbps = 6')