%% Script para demostración de la estimación no paramétrica por kernels JIRV
%
close all
clear all
path(path,'C:\Users\josei\Desktop\RecuperadosHPviejita\Kernelest');
n = 50;             % Longitud de la muestra N = n, comenzar con 50, luego 500,
                   % 5000, y 10000: ver que sucede
%err = randn(n,1);  % Sigue una distribución N(0,1)
x1 = 1/4 + 1/10*randn(n/2,1); 
x2 = 3/4 + 1/10*randn(n/2,1);
err = [1.5*x1 x2];  % Sigue una distribución P = 1.5*N(1/4,(1/10)^2)+N(3/4,(1/10)^2)
err = err(:);       % P Bimodal como en el libro de Bishop pág. 121

[x,y] = hist(err,50);  % Histograma con 50 Bins                        
%figure,bar(y,x/trapz(y,x),'y')                 
                            
hn = bwsnr(err)

errest = Kerestn(err,hn,length(err));  % Kernel exponencial o Gaussiano
l = linspace(min(err),max(err),length(errest));
figure(1),plot(l,errest),title('Kernel exponencial hn metodo bwsnr')
xlabel('Distribución de probabilidad'),ylabel('Amplitud de probabilidad')
                         
errestd = Kerestnd(err,hn,length(err));  % Kernel exponencial normalizado
l = linspace(min(err),max(err),length(errestd));
hold on, plot(l,errestd,'+ g')                          

hn2 = bwos(err)

errest2 = Kerestn(err,hn2,length(err));  % Kernel exponencial con hn2
l = linspace(min(err),max(err),length(errest2));
figure(2),plot(l,errest2,'k'),title('Kernel exponencial hn metodo bwos')
xlabel('Distribución de probabilidad'),ylabel('Amplitud de probabilidad')                          

hn3 = bwsjpib(err)
  
errest3 = Kerestn(err,hn3,length(err));   % Kernel exponencial con hn3
l = linspace(min(err),max(err),length(errest3));
figure(3),plot(l,errest3,'m'),title('Kernel exponencial hn metodo bwsjpib')
xlabel('Distribución de probabilidad'),ylabel('Amplitud de probabilidad')

[kde,lx,mker] = gpkde(err,0);        % Kernel exponencial con hn optimo                  
figure(4),plot(lx,kde),title('Kernel exponencial hn metodo optimizado')
xlabel('Distribución de probabilidad'),ylabel('Amplitud de probabilidad')

errest5 = Kerestco(err,40,length(err));  % Kernel cosenoidal
l = linspace(min(err),max(err),length(errest5));
figure(5),plot(l,errest5,'g'),title('Kernel cosenoidal con exponente 40')
xlabel('Distribución de probabilidad'),ylabel('Amplitud de probabilidad')


errest6 = KerestH(err,1,length(err));    % Kernel de Hilbert
l = linspace(min(err),max(err),length(errest6));
figure(6),plot(l,errest6,'y'),title('Kernel de Hilbert, no calcula hn')
xlabel('Distribución de probabilidad'),ylabel('Amplitud de probabilidad')

%------------------------------------
% Comparacion de diferentes kernels:
%------------------------------------

figure(7),bar(y,x/trapz(y,x),'r'),title('Comparación entre Kernels')
hold on,plot(lx,kde,'b'),
hold on,plot(l,errest2,'k')
hold on,plot(l,errest5,'g')
hold on,plot(l,errest6,'y')
xlabel('Distribución de probabilidad'),ylabel('Amplitud de probabilidad')