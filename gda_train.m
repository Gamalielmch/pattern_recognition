%%%%%%%%%%%%%%%%%%%%% LOAD and plot
% clear all


function [mu,sigma,phi]=gda_train()
load fisheriris

x=meas;
n_class=3; %múmero de clases
n_charac=size(x,2); %número de características
y=zeros(length(x),1);
%etiqueta setosa=0 versicolor=1 virginica=2
match  = ismember(species, 'versicolor');
y(match==1)=1;
match  = ismember(species, 'virginica');
y(match==1)=2;

%tamaño de muestra
m=length(y);

%media
media=mean(x);

%removiendo media
xo=x-media;

%aplicando pca para visualizar
[coeff,score,latent] = pca(xo);
x_setosa=score(y==0,:);
x_versi=score(y==1,:);
x_vergi=score(y==2,:);
%graficando PCA 1,2 y 3
Xc = score*coeff';
figure
hold off
plot3(x_setosa(:,1),x_setosa(:,2),x_setosa(:,3),'r+')
hold on
plot3(x_versi(:,1),x_versi(:,2),x_versi(:,3),'xm')
plot3(x_vergi(:,1),x_vergi(:,2),x_vergi(:,3),'ob')
%biplot(coeff(:,1:2),'scores',score(:,1:2),'varlabels',{'v_1','v_2','v_3','v_4'});

%%%%%%%%%%%%%%%%%%%%%%%% GDA
X_mu=[];
phi=zeros(n_class,1);
mu=zeros(n_class,n_charac);
for i=1:n_class
    phi(i)=length(find(y==i-1))/m;
    mu(i,:)=sum(x(y==i-1,:))/length(find(y==i-1));
    X_mu=[X_mu; x(y==i-1,:)-mu(i,:)];
end

sigma=(X_mu'*X_mu)/m;
end


