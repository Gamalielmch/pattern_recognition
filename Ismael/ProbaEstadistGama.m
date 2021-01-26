%
%  Comprobaci�n de algunos conceptos
%  Distribuci�n Gama
format long
close all
clear all

A = 1;  % par�metos de la distribuci�n
B = 2;
M = 1;
N = 500;
X = gamrnd(A,B,M,N);

figure,plot(X),title('V. A. muestreada de una distribuci�n Gamma')

% Trazado de hostograma, distancia inter bins
%
[ejex,ejey]= hist(X,50);
figure,bar(ejey,ejex/trapz(ejey,ejex),'r'), axis([-1 8 0 3]),
title('Histograma de una V. A. Gamma')
mX = mean(X)
varX = cov(X)
momSkew = skewness(X)
momkurt = kurtosis(X)