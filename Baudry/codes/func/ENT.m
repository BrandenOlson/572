function f=ENT(y,mu,S,p,K,n,tm)

switch nargin

    case 4
        tm=min(200,n);
        
        if size(y,1)>tm
            N=floor(size(y,1)/tm);
        else
            N=0;
        end
        
        f=zeros(N+1);
        
        for i=1:N
            
            post=posteriorf(y(((i-1)*tm+1):(i*tm),:),mu,S,p);

            f(i)=-sum(sum(post.*log_perso(post)));
        end
        
        if N<(size(y,1)/tm)
            post=posteriorf(y((N*tm+1):size(y,1),:),mu,S,p);
            f(N+1)=-sum(sum(post.*log_perso(post)));
        end
        
        f=sum(f);
        
        f=f(1);
        
    case 5
        
        tm=min(200,n);
        
        if size(y,1)>tm
            N=floor(size(y,1)/tm);
        else
            N=0;
        end
        
        f=zeros(N+1);
        
        for i=1:N
            
            post=posteriorf(y(((i-1)*tm+1):(i*tm),:),mu,S,p,K);

            f(i)=-sum(sum(post.*log_perso(post)));
        end
        
        if N<(size(y,1)/tm)
            post=posteriorf(y((N*tm+1):size(y,1),:),mu,S,p,K);
            f(N+1)=-sum(sum(post.*log_perso(post)));
        end
        
        f=sum(f);
        
        f=f(1);

    case 6
        
        tm=min(200,n);

        if n>tm
            N=floor(n/tm);
        else
            N=0;
        end
        
        f=zeros(N+1,1);
        
        for i=1:N
            
            post=posteriorf(y(((i-1)*tm+1):(i*tm),:),mu,S,p,K,tm);

            f(i)=-sum(sum(post.*log_perso(post)));
            
        end
        
        if N<(n/tm)
            post=posteriorf(y((N*tm+1):n,:),mu,S,p,K,n-N*tm);
            f(N+1)=-sum(sum(post.*log_perso(post)));
        end
   
        f=sum(f);
        f=f(1);
        
    case 7
        
        if n>tm
            N=floor(n/tm);
        else
            N=0;
        end
        
        f=zeros(N+1,1);
        
        for i=1:N
            
            post=posteriorf(y(((i-1)*tm+1):(i*tm),:),mu,S,p,K,tm);

            f(i)=-sum(sum(post.*log_perso(post)));
            
        end
        
        if N<(n/tm)
            post=posteriorf(y((N*tm+1):n,:),mu,S,p,K,n-N*tm);
            f(N+1)=-sum(sum(post.*log_perso(post)));
        end
   
        f=sum(f);
        f=f(1);
        
end
        
