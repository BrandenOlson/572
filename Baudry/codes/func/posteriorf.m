function f=posteriorf(y,mu,S,p,K,n,tm)

switch nargin
    case 6
        tm=min(n,50);
end

f=zeros(n,K);

N=floor(n/tm);

for i=1:N
    yN=y((i-1)*tm+1:i*tm,:);
    f((i-1)*tm+1:i*tm,:)=posteriorf0(yN,mu,S,p,K,tm);
end

if N<(n/tm)
    yN=y(N*tm+1:n,:);
    f(N*tm+1:n,:)=posteriorf0(yN,mu,S,p,K,n-N*tm);
end

end
