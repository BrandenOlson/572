function z=MAP_combi(y,mu,S,p,M,n)

% MAP_combi(y,mu,S,p,M,n), with y the data (n*d matrix), 
% mu the expectations of the Gaussian components (K*d matrix), S the 
% covariance matrices of the Gaussian components (d*d*K matrix), p the 
% mixing proportions (1*d vector) and M the combining matrix (K*K_combi 
% matrix of 0 and 1), returns the labels of the combined solution, computed
% by the maximum a posteriori method.

if nargin==5
    n=size(y,1);
end

K=size(M,2);
K_combi=size(M,1);

z=zeros(n,K_combi);
poster=zeros(n,K);

for k=1:K
    poster(:,k)=p(k)*mvnpdf(y,mu(k,:),S(:,:,k));
end
    
poster=(M*poster')';

z=(repmat(max(poster'),K_combi,1)'==poster);
