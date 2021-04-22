%
%  Comprobaci�n de algunos conceptos
%  Distribuci�n Beta
clc
format long
close all
clear all

A = 1;  % par�metos de la distribuci�n
B = 2;
M = 1;
N = 5000000;
X = betarnd(A,B,M,N);

figure,plot(X),title('V.A. muestreada de una distribuci�n Beta')

% Trazado de hostograma, distancia inter bins
%
[ejex,ejey] = hist(X,14);
figure,bar(ejey,ejex/trapz(ejey,ejex),'r'), axis([-3 3 0 3]),
title('Histograma de una V.A. Beta')
mX = mean(X)
varX = cov(X)
momSkew = skewness(X)
momkurt = kurtosis(X)