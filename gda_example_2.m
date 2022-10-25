%%%%%%%%%%%%%%%%%%%%% LOAD and plot
clear all
load fisheriris

x=meas;
nc=3; %múmero de clases
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
for i=1:nc-1
phi(i)=length(find(y==i-1))/m;
mu(i,:)=sum(x(y==i-1,:))/length(find(y==i-1));
X_mu=[X_mu; x(y==i-1,:)-mu(i,:)];
end

sigma=(X_mu'*X_mu)/m;
[x1,y1]=meshgrid(linspace(-3,8,100)',linspace(-3,8,100)');
X1=[x1(:),y1(:)];
z1=mvnpdf(X1,mu_1,sigma);
contour(x1,y1,reshape(z1,100,100),8)
hold on
z2=mvnpdf(X1,mu_0,sigma);
contour(x1,y1,reshape(z2,100,100),8)

z2=mvnpdf(X1,mu_0,sigma);
contour(x1,y1,reshape(z2,100,100),8)


