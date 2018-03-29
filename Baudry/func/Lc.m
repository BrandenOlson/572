function f=Lc(y,mu,S,p,K,n)

switch nargin

    case 4
        
        post=posteriorf(y,mu,S,p);

        f=L(y,mu,S,p)+sum(sum(post.*log_perso(post)));
    
    case 5
        
        post=posteriorf(y,mu,S,p,K);

        f=L(y,mu,S,p,K)+sum(sum(post.*log_perso(post)));

    case 6
        
        post=posteriorf(y,mu,S,p,K,n);

        f=L(y,mu,S,p,K,n)+sum(sum(post.*log_perso(post)));
        
end
