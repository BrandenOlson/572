function [mu,S,p]=extract_param(par,K,d)

switch nargin
    case 1
        K=size(par.mean,2);
        d=size(par.mean{1},2);
    case 2
        d=size(par.mean{1},2);
end

mu=zeros(K,d);
S=reshape(zeros(d*d*K,1),d,d,K);
p=zeros(1,K);

for k=1:K
    mu(k,:)=par.mean{k};
    S(:,:,k)=par.dispersion{k};
    p(k)=par.proportion(k);
end
