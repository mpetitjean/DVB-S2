clear;
close all;

i = 1;
min = 10;
max = 100;
step = 10;
for ppm = min:step:max
    ber(i,:) = berCFO(4, 1e4, -5, 1, 15, ppm, 0);
    disp(i)
    i = i+1;
end

figure;
ppm = min:step:max;
for i = 1:length(ppm)
    % Je sais Ced c'est degueu
    semilogy(-5:15, ber(i,:))
    hold on;
    legendInfo{i} = ['ppm = ' num2str(ppm(i))];
end
legend(legendInfo);