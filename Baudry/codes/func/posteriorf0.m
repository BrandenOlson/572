function [f,f0]=posteriorf0(y,mu,S,p,K,n)

%POSTERIOR Calcule la matrice des proba a posteriori des composantes d'un
%mélange de gaussiennes, pour un échantillon.
%   POSTERIOR(y,mu,S,p) calcule la matrice n*K (n=taille de l'échantillon
%   y; K=nombre de composantes du mélange=nbre de lignes de mu) des proba a
%   posteriori que les observations y_i aient été tirées selon la kième
%   composante du mélange.
%   
%   POSTERIOR(y,mu,S,p,K) (resp. POSTERIOR(y,mu,S,p,K,n)) fait de même sans
%   recalculer K (resp. K et n).

if nargin==4
    K=size(mu,1);
    n=size(y,1);
elseif nargin==5
    n=size(y,1);
end

d=ones(1,K);

for k=1:K
    d(1,k)=det(S(:,:,k));
end

if any(d<10^(-300))
        
        f=NaN;
        
else

    post=zeros(n,K);

    for k=1:K
        post(:,k)=p(k)*mvnpdf(y,mu(k,:),S(:,:,k));
    end

    f=zeros(size(post));

    di=sum(post,2);

    f_zero=(di<10^(-20));
        
    f(~f_zero,:)=diag(di(~f_zero))^(-1)*post(~f_zero,:);
    f_inc=(post(f_zero,:)==repmat(max(post(f_zero,:)')',1,K));
    f(f_zero,:)=f_inc;
    f(sum(f')'>1.5,:)=repmat([1 zeros(1,K-1)],size(f(sum(f')'>1.5,:),1),1);
    
end

end
