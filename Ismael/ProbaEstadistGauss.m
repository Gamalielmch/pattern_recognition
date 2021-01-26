%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Comprobación de algunos conceptos
%  Distribución Normal
clc
format long
close all
clear all

MU = 1;     % parametros de la distribución
SIGMA = 2;  % MU es el valor medio y Var = SIGMA^2
M = 1;      % dimensiones de la V.A.
N = 5000;
X = normrnd(MU,SIGMA,M,N);

figure(1),plot(X),title('V.A. obtenida a partir de una distribución Normal/Gaussiana')

% Trazado de hostograma, distancia inter bins
%
[ejex,ejey] = hist(X,20);
figure(2),bar(ejey,ejex/(N*.5),'r'), axis([-7 7 0 3]),
title('Histograma de una V.A. Normal/Gaussiana')
meanX = mean(X)
varX = cov(X)
momSkew = skewness(X)
momkurt = kurtosis(X)

y = [];
for i=1:10,
    X = normrnd(MU,SIGMA,M,N);
    y = [y X];
end
meanProcY = mean(y)
varProcY = cov(y)

figure(3),plot(y),title('Proceso estocástico Gaussiano Y')
% Trazado de histograma, distancia inter bins
%
[ejex,ejey] = hist(y,20);
figure(4), bar(ejey,ejex/(length(y)*.5),'g'), axis([-7 7 0 3]),
title('Histograma del proceso Gaussiano Y')

Y = normrnd(MU,SIGMA,2,N);
figure(5),plot(Y(1,:),Y(2,:),'*r'),title('Procesos estocásticos Gaussianos X e Y')
meanProcXY = mean(mean(Y))
varProcXY = cov(Y(1,:),Y(2,:))
figure(6),hist(mean(Y)),title('Histograma de la media del proceso X e Y')

figure(7),surf(Y),title('Procesos estocásticos Gaussianos X e Y')
figure(8),hist(Y),title('Histograma de los Procesos estocásticos Gaussianos X e Y')
