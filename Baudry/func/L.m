function f=L(y,mu,S,p,K,n)

if nargin==4
    K=size(mu,1);
    n=size(y,1);
elseif nargin==5
    n=size(y,1);
end

d=1;

for k=1:K
    d=d*det(S(:,:,k));
end

if d<10^(-300)
        
        f=NaN;
        
else

    Lp=zeros(n,K);

    for k=1:K
        Lp(:,k)=p(k)*mvnpdf(y,mu(k,:),S(:,:,k));
    end

    f=sum(log(sum(Lp,2)));

end