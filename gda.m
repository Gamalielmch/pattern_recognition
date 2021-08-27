function gda
load fisheriris
X=meas;
y=zeros(length(X),1);
%etiqueta setosa=0 versicolor=1 virginica=2
match  = ismember(species, 'versicolor');
y(match==1)=1;
match  = ismember(species, 'virginica');
y(match==1)=2;
[mu,sigma,phi]=gda_train(X,y);
x=[1,2,3,4];
probs=gda_test(x,mu,sigma,phi)
end


function [mu,sigma,phi]=gda_train(X,y)
n_class=length(unique(y)); %número de clases
n_charac=size(X,2); %número de características
m=length(y);
X_mu=[];
phi=zeros(n_class,1);
mu=zeros(n_class,n_charac);
for i=1:n_class
    phi(i)=length(find(y==i-1))/m;
    mu(i,:)=sum(X(y==i-1,:))/length(find(y==i-1));
    X_mu=[X_mu; X(y==i-1,:)-mu(i,:)];
end

sigma=(X_mu'*X_mu)/m;
end

function probs=gda_test(x,mu,sigma,phi)
probs=mvnpdf(x,mu,sigma).*phi;
probs=probs/sum(probs);
end



