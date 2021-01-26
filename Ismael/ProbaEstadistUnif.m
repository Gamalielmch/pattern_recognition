%
%  Comprobación de algunos conceptos
%  Distribución uniforme
clc
format long
close all
clear all

A = -2;      % parametros de la distribución
B = 2;
M = 1;
N = 50000;
X = unifrnd(A,B,M,N);

figure(1), plot(X,'r'),title('V.A. obtenida a partir de una distribución Uniforme U(A,B)')

% Trazado de hostograma, distancia inter bins
%
[ejex,ejey] = hist(X,20);
figure(2), bar(ejey,ejex/trapz(ejey,ejex),'c'), axis([-3 3 0 3]),title('Histograma de una V.A. Uniforme')
meanX = mean(X)
varX = cov(X)