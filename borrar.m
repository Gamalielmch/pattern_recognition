clc
close all
format long
x = [0.99 1.02 1.15 1.29 1.46 1.36 0.87 1.23 1.55 1.40 1.19 1.15 0.98...
    1.01 1.11 1.20 1.26 1.32 1.43 0.95];
y = [90.1 89.05 91.43 93.74 96.73 94.45 87.59 91.77 99.42 93.65 93.54...
    92.52 90.56 89.54 89.85 90.39 93.25 93.41 94.98 87.33];
N = length(y);
Sumx = sum(x);
Sumx2 = sum(x.^2);
Sumy = sum(y);
Sumy2 = sum(y.^2);
Sumxy = sum(x.*y);
a = (N*Sumxy-Sumy*Sumx)/(N*Sumx2- Sumx^2);
b = (Sumy*Sumx2-Sumx*Sumxy)/(N*Sumx2- Sumx2^2);
Y = a*x+b;
disp(a);
disp(b);
figure(1),plot(x,y,'o'),hold on,
xlabel('Eje x: Nivel de hdrocarburos (%)'),
ylabel('Eje y: Pureza de oxígeno (%)'),grid,
title('Datos de niveles de oxígeno de hidrocarburos')
plot(x,Y)