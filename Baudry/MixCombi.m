%% Computes the combined solutions, starting from the BIC
%% solution and ending with the solution with a single class, according to
%% the article.

% Variables necessarily provided by the user:

% exp.data.obs (data)
% exp.res.BIC.mu 
% exp.res.BIC.S 
% exp.res.BIC.p 
% exp.res.BIC.K

% Optional variables: 
% exp.res.ICL.mu
% exp.res.ICL.S
% exp.res.ICL.p
% exp.res.ICL.K

exp.cond.n=max(size(exp.data.obs));

addpath('func');

% To evaluate the cpu time:
tic

%Solution (of BIC) before combining:
exp.res.combi(exp.res.BIC.K).M=ones(exp.res.BIC.K,exp.res.BIC.K);
exp.res.combi(exp.res.BIC.K).labels=MAP_combi(exp.data.obs,exp.res.BIC.mu,exp.res.BIC.S,exp.res.BIC.p,eye(exp.res.BIC.K,exp.res.BIC.K),exp.cond.n);%(=exp.res.BIC.labels)
exp.res.combi(exp.res.BIC.K).tau=posteriorf(exp.data.obs,exp.res.BIC.mu,exp.res.BIC.S,exp.res.BIC.p,exp.res.BIC.K,exp.cond.n); %(=exp.res.BIC.tau)
exp.res.BIC.ent=ENT(exp.data.obs,exp.res.BIC.mu,exp.res.BIC.S,exp.res.BIC.p,exp.res.BIC.K,exp.cond.n);

%Combining:

for K=exp.res.BIC.K-(1:exp.res.BIC.K-1)

    % Try to combine each pair of components (j,k) and choose the one which
    % resulting entropy is minimal

    ent_current=inf;
    
    for k=1:K
        for j=(k+1):(K+1)

            M1=zeros(j-1,K-j+2);
            M1(k,1)=1;

            M=[eye(j-1,j-1) M1;...
                zeros(K-1-j+2,j-1) zeros(K-1-j+2,1) eye(K-j+1,K-j+1)];

            tau=(M*exp.res.combi(K+1).tau')';

            ent=-sum(sum(tau.*log_perso(tau)));

            if ent>ent_current
                continue
            end

            ent_current=ent;

            exp.res.combi(K).M=M;
            exp.res.combi(K).tau=tau;
            exp.res.combi(K).ent=ent;

            M_tot=eye(K,K);

            for KM=K:(exp.res.BIC.K-1)
                M_tot=M_tot*exp.res.combi(KM).M;
            end

            exp.res.combi(K).labels=MAP_combi(exp.data.obs,exp.res.BIC.mu,exp.res.BIC.S,exp.res.BIC.p,M_tot,exp.cond.n);
        end
    end
end

%Put in cell res the chosen models and display it

clear('res');
res(1,:)={'Criterion' 'K' 'ENT' 'Lcc'};
try
    res(2,:)={'ICL' int2str(exp.res.ICL.K) int2str(ENT(exp.data.obs,exp.res.ICL.mu,exp.res.ICL.S,exp.res.ICL.p,exp.res.ICL.K,exp.cond.n)) int2str(Lc(exp.data.obs,exp.res.ICL.mu,exp.res.ICL.S,exp.res.ICL.p,exp.res.ICL.K,exp.cond.n))};
catch % In case the user did not provide the ICL informations
    res(2,:)={'ICL' '???' '???' '???'};
end    
for i=3:exp.res.BIC.K
    res(i,:)={['Combined K=' int2str(i+2-3)] int2str(i+2-3) int2str(exp.res.combi(i+2-3).ent) int2str(L(exp.data.obs,exp.res.BIC.mu,exp.res.BIC.S,exp.res.BIC.p,exp.res.BIC.K,exp.cond.n)-exp.res.combi(i+2-3).ent)};
end
res(exp.res.BIC.K+1,:)={'BIC' int2str(exp.res.BIC.K) int2str(exp.res.BIC.ent) int2str(Lc(exp.data.obs,exp.res.BIC.mu,exp.res.BIC.S,exp.res.BIC.p,exp.res.BIC.K,exp.cond.n))};

disp(res)
 