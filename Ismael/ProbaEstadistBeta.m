%
%  Comprobación de algunos conceptos
%  Distribución Beta
clc
format long
close all
clear all

A = 1;  % parámetos de la distribución
B = 2;
M = 1;
N = 5000000;
X = betarnd(A,B,M,N);

figure,plot(X),title('V.A. muestreada de una distribución Beta')

% Trazado de hostograma, distancia inter bins
%
[ejex,ejey] = hist(X,14);
figure,bar(ejey,ejex/trapz(ejey,ejex),'r'), axis([-3 3 0 3]),
title('Histograma de una V.A. Beta')
mX = mean(X)
varX = cov(X)
momSkew = skewness(X)
momkurt = kurtosis(X)