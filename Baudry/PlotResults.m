addpath('func');

%% Plot the observed data

figure;

plot(exp.data.obs(:,1),exp.data.obs(:,2),'*b');

title(['Observed Data. n=' int2str(exp.cond.n)]);

%% Plot the BIC solution

figure;

map_contour(exp.data.obs,exp.res.BIC.mu,exp.res.BIC.S,exp.res.BIC.p,exp.res.BIC.K);

title(['BIC Solution. K=' int2str(exp.res.BIC.K) '. ENT=' int2str(exp.res.BIC.ent)]);
box on;

%% Plot the combined solutions

y=exp.data.obs;

x_lim=[min(y(:,1))-2 max(y(:,1))+2];
y_lim=[min(y(:,2))-2 max(y(:,2))+2];

X=(min(x_lim):(max(x_lim)-min(x_lim))/100:max(x_lim))';
Y=(min(y_lim):(max(y_lim)-min(y_lim))/100:max(y_lim))';

pdf=zeros(exp.res.BIC.K,101,101);

for kpdf=1:exp.res.BIC.K
    for i=1:101
        pdf(kpdf,:,i)=exp.res.BIC.p(kpdf)*mvnpdf([X,repmat(Y(i),101,1)],exp.res.BIC.mu(kpdf,:),exp.res.BIC.S(:,:,kpdf));
    end
end

color=['r';'g';'b';'c';'m';'y';'k';'g';'r';'b';'c';'m';'y';'k'];
col=color(1:exp.res.BIC.K);

for Kb=2:exp.res.BIC.K-1
    K=exp.res.BIC.K+1-Kb;
    for k=1:K
       coln(k)=col(min(find(exp.res.combi(K).M(k,:)==1)));
    end
    col=coln;
    figure;
  
    hold on;

    xlabel('X_1');
    ylabel('X_2');

    levels=[0.005 0.005];

    M=eye(K,K);

    for KM=K:(exp.res.BIC.K-1)
        M=M*exp.res.combi(KM).M;
    end
    
    pdf_mix=zeros(K,101,101);
    
    for i=1:101
        for j=1:101
            pdf_mix(:,i,j)=M*pdf(:,i,j);
        end
    end

    for k=1:K
        yl=y(find(exp.res.combi(K).labels(:,k)==1),:);
        plot(yl(:,1),yl(:,2),[col(k) '.']);

        Z=zeros(101,101);
        for i=1:101
            Z(:,i)=pdf_mix(k,:,i)';
        end
        cZ=max(Z);
        c_max=find(max(cZ)==cZ);
        r_max=find(Z(:,c_max)==max(Z(:,c_max)));
        
        contour(X,Y,Z',levels,col(k));
    end

    title(['Combined Solution. K=' int2str(K) '. ENT=' int2str(exp.res.combi(K).ent)]);
    box on;

    hold off;

end

%% Plot the ICL solution (optional)

try % In case the user did not provide the ICL solution
    
    figure;

    map_contour(exp.data.obs,exp.res.ICL.mu,exp.res.ICL.S,exp.res.ICL.p,exp.res.ICL.K);

    title(['ICL solution. K=' int2str(exp.res.ICL.K) '. ENT=' int2str(ENT(exp.data.obs,exp.res.ICL.mu,exp.res.ICL.S,exp.res.ICL.p,exp.res.ICL.K,exp.cond.n))]);
    box on;

catch
end