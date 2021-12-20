% Matlab File: plotn.m
% This program is to be used after impn.m to plot impulse response functions
% Order of variables is: Date K A G Lambda C N Y Y/N I R


time = 1:nir;

subplot(221), plot(time,IR(:,3),time,IR(:,4),time,IR(:,2))
title('A, G, K')
xlabel('time')
legend('A','G','K')

subplot(222), plot(time,IR(:,6),time,IR(:,10),time,IR(:,8)) 
title('C, I, Y ')
xlabel('time')
legend('C','I','Y')

subplot(223), plot(time,IR(:,7),time,IR(:,9))
title('Labor, Wages')
xlabel('time')

subplot(224), plot(time,IR(:,11))
title('Interest Rate')
xlabel('time')
legend('r')

clear time

