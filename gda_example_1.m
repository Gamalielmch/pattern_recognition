%%%%%%%%%%%%%%%%%%%%%%%% GDA
m=200;
n=2;
rp=mvnrnd([1,1],[1 0; 0 1],m/2);
rn=mvnrnd([4,4],[1 0; 0 1],m/2);
y=[ones(m/2,1); zeros(m/2,1)];
figure, hold off
plot3(rp(:,1),rp(:,2),y(1:m/2,1),'b+')
hold on
plot3(rn(:,1),rn(:,2),y(m/2+1:m,1),'ro')

phi=length(find(y==1))/m;
mu_0=sum(rn)/length(find(y==0));
mu_1=sum(rp)/length(find(y==1));
X=[rp;rn];
X_mu1=X(y==1,:)-mu_1;
X_mu2=X(y==0,:)-mu_0;
X_mu=[X_mu1;X_mu1];
sigma=(X_mu'*X_mu)/m;
[x1,y1]=meshgrid(linspace(-3,8,100)',linspace(-3,8,100)');
X1=[x1(:),y1(:)];
z1=mvnpdf(X1,mu_1,sigma);
contour(x1,y1,reshape(z1,100,100),8)
hold on
z2=mvnpdf(X1,mu_0,sigma);
contour(x1,y1,reshape(z2,100,100),8)


