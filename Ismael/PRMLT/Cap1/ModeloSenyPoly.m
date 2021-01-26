% Modelo polinomial vs sin(2pix)
format long
close all
clear all

x = [0:0.001:1];
N = max(size(x));
TrueMo = sin(2*pi*x);
SIGMA = 0.2;
tn = TrueMo + normrnd(0,SIGMA,1,N);

figure(1), plot(x,TrueMo),title('Modelo Senoidal y con Ruido')
hold on, plot(x,tn,'o r')

[ejex,ejey] = hist(tn,20);
figure(2),bar(ejey,ejex/(N*.5),'r'), axis([-7 7 0 3]),
title('Histograma de t_n como V.A. Normal/Gaussiana?')
meanX = mean(tn)
varX = cov(tn)
momSkew = skewness(tn)
momkurt = kurtosis(tn)

% Calculos de sumatorias para el sist. de ecuaciones normales: polinomio de
% orden 2 : p_2 (x)

E1a0 = max(size(x));
E1a1 = sum(x);
E1a2 = sum(x.^2);
E2a0 = E1a1;
E2a1 = E1a2;
E2a2 = sum(x.^3);
E3a0 = E2a1;
E3a1 = E2a2;
E3a2 = sum(x.^4);

y1 = sum(tn);
y2 = sum(tn.*x);
y3 = sum(tn.*x.^2);

% Matriz del sistema (Regresor) lineal de ecuaciones para Minimos Cuadrados

A = [E1a0 E1a1 E1a2; E2a0 E2a1 E2a2; E3a0 E3a1 E3a2];
y = [y1; y2; y3];

% Estimacion de Minimos Cuadrados para los a2, a1, y a0

aest = A\y

f  = aest(3)*x.^2 + aest(2)*x + aest(1);

% Error cuadratico prom:

SumEx = sum((tn - f).^2)/E1a0

% 
% Dibujo de la comparacion entre la salida original y la salida estimada 
% 

figure,plot(x,TrueMo),hold on,plot(x,tn,'g*'),hold on,plot(x,f,'r'),title('Salida con ruido (asteriscos verdes) y Salida Estimada (llena rojo), Verdadeo (azul)');

% Calculos de sumatorias para el sist. de ecuaciones normales: polinomio de
% orden 3: p_3 (x)

E1a0 = max(size(x));
E1a1 = sum(x);
E1a2 = sum(x.^2);
E1a3 = sum(x.^3);
E2a0 = E1a1;
E2a1 = E1a2;
E2a2 = E1a3;
E2a3 = sum(x.^4);
E3a0 = E2a1;
E3a1 = E2a2;
E3a2 = E2a3;
E3a3 = sum(x.^5);
E4a0 = E3a1;
E4a1 = E3a2;
E4a2 = E3a3;
E4a3 = sum(x.^6);

y1 = sum(tn);
y2 = sum(tn.*x);
y3 = sum(tn.*x.^2);
y4 = sum(tn.*x.^3);

% Matriz del sistema (Regresor) lineal de ecuaciones para Minimos Cuadrados

A = [E1a0 E1a1 E1a2 E1a3; E2a0 E2a1 E2a2 E2a3; E3a0 E3a1 E3a2 E3a3; E4a0 E4a1 E4a2 E4a3];
y = [y1; y2; y3; y4];

% Estimacion de Minimos Cuadrados para los a3, a2, a1, y a0

aest = A\y

f  = aest(4)*x.^3 + aest(3)*x.^2 + aest(2)*x + aest(1);

% Error cuadratico prom:

SumEx = sum((tn - f).^2)/E1a0

% 
% Dibujo de la comparacion entre la salida original y la salida estimada 
% 

figure,plot(x,TrueMo),hold on,plot(x,tn,'g*'),hold on,plot(x,f,'r'),title('Salida con ruido (asteriscos verdes) y Salida Estimada (llena rojo), Veradero (azul)');