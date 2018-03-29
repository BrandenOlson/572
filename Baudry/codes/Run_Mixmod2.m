%% Computes the MLE's solutions for K between exp.cond.Kmin and
%% exp.cond.Kmax, the corresponding BIC and ICL values, and returns their
%% solutions such that they can be read by the MixCombi script. 

%% The following variables are necessary

% exp.cond.Kmax (maximum number of gaussian components to be considered)
% exp.cond.Kmin (minimum number of gaussian components to be considered)
% exp.cond.models (Models to be condidered)
% exp.data.obs (Data)

%% EM with Mixmod...

% Mixmod 2.1.1 paths

addpath('c:\Program Files\Mixmod-2.1.1\')
addpath('c:\Program Files\Mixmod-2.1.1\UTIL\MATLAB')

addpath('func');

% First mixmod launch to get parameters and solutions for each K and model
% Initialisation and strategy parameters for mixmod

init= struct('name','RANDOM','param',[],'partition',{{}});
algo = struct('name','EM','stopRule','NBITERATION','stopRuleValue',10);
strategy = struct('initialization',init,'algorithm',algo);
mod=struct('name','','subDimensionFree',[],'subDimensionEqual',[]);
criteria={'BIC'};
mod.name=exp.cond.models{1};

% To construct the 'mix' structure 

output=mixmod(exp.data.obs,1,'criterion',criteria,'model',{mod},'strategy',strategy);

mix=repmat(output,exp.cond.Kmax-exp.cond.Kmin+1,size(exp.cond.models,1));

% Mixmod for each (k,model), to get dimensions and parameters in each case

exp.cond.n=max(size(exp.data.obs));
exp.cond.d=min(size(exp.data.obs));

for k=exp.cond.Kmin:exp.cond.Kmax
    for nmodel=1:size(exp.cond.models,1)
        output=[];
        while isequal(output,[])
            output=mixmod(exp.data.obs,k,'criterion',criteria,'model',{mod},'strategy',strategy);
        end
        mix(k,nmodel).dim=output.modelOutput(1).likelihood.nbFreeParam;
        mix(k,nmodel).param=output.modelOutput(1).param;
     end
end

% And then, a better run of mixmod to get the true choices of ICL and
% BIC

% Solution for BIC and ICL among all possible K and models

BIC_value=repmat(inf,exp.cond.Kmax-exp.cond.Kmin+1,size(exp.cond.models,1));
ICL_value=repmat(inf,exp.cond.Kmax-exp.cond.Kmin+1,size(exp.cond.models,1));

i_actuel=1;
K_actuel=exp.cond.Kmin;

while i_actuel<2
    for n_models=1:size(exp.cond.models,1)
        mod=struct('name','','subDimensionFree',[],'subDimensionEqual',[]);
        mod.name=exp.cond.models{n_models};
        while K_actuel<exp.cond.Kmax+1
                % Initialisation and strategy parameters for mixmod (with a better and more
                % expensive strategy)
                % Rem : The Small-EM strategy is implemented by hand to match the one described in 
                % Biernacki, C., Celeux, G., and Govaert, G. (2003) : "Choosing starting values for 
                % the em algorithm for getting the highest likelihood in
                % multivariate Gaussian mixture models.", 
                % Computational Statistics & Data Analysis, 41(3-4):567 --
                % 575.
                init= struct('name','RANDOM','param',[],'partition',{{}});
                algo = struct('name','EM','stopRule','NBITERATION','stopRuleValue',50);
                strategy = struct('initialization',init,'algorithm',algo);
                criteria={'BIC';'ICL'};
                
            % A few random initializations followed by a short run of EM
            for i=1:10
                output=[];
                for j=1:10 
                    output=mixmod(exp.data.obs,K_actuel,'criterion',criteria,'model',{mod},'strategy',repmat(strategy,1,1));
                    if isequal(output,[])
                        continue; 
                    else
                        break;
                    end
                end
                if output.modelOutput(1).criterion.value<BIC_value(K_actuel-exp.cond.Kmin+1,n_models)
                    out=output;
                    BIC_value(K_actuel-exp.cond.Kmin+1,n_models)=output.modelOutput(1).criterion.value;
                    [exp.res.mix.mu{K_actuel,n_models},exp.res.mix.S{K_actuel,n_models},exp.res.mix.p{K_actuel,n_models}]=extract_param(output.modelOutput(1).param,K_actuel,exp.cond.d);
                    % ICL_value(K_actuel-exp.cond.Kmin+1,n_models)=output.modelOutput(2).criterion.value;
                    % The computation of ICL is made a little differently 
                    % than mixmod does.
                    ICL_value(K_actuel-exp.cond.Kmin+1,n_models)=BIC_value(K_actuel-exp.cond.Kmin+1,n_models)+2*ENT(exp.data.obs,exp.res.mix.mu{K_actuel,n_models},exp.res.mix.S{K_actuel,n_models},exp.res.mix.p{K_actuel,n_models},K_actuel,exp.cond.n);
                end
            end

            init= struct('name','USER','param',[],'partition',{{}});
            init.param=out.modelOutput(1).param;
            algo = struct('name','EM','stopRule','NBITERATION','stopRuleValue',2000);
            strategy = struct('initialization',init,'algorithm',algo);
            output=[];
          % And a long run of EM from the best of the preceding solutions  
          for j=1:10
                output=mixmod(exp.data.obs,K_actuel,'criterion',criteria,'model',{mod},'strategy',strategy);
                if isequal(output,[])
                    continue;
                else
                    if output.modelOutput(1).criterion.value<BIC_value(K_actuel-exp.cond.Kmin+1,n_models) % To work out the NaN's
                        BIC_value(K_actuel-exp.cond.Kmin+1,n_models)=output.modelOutput(1).criterion.value;
                        [exp.res.mix.mu{K_actuel,n_models},exp.res.mix.S{K_actuel,n_models},exp.res.mix.p{K_actuel,n_models}]=extract_param(output.modelOutput(1).param,K_actuel,exp.cond.d);
                        % ICL_value(K_actuel-exp.cond.Kmin+1,n_models)=output.modelOutput(2).criterion.value;
                        % The computation of ICL is made a little differently 
                        % than mixmod does.
                        ICL_value(K_actuel-exp.cond.Kmin+1,n_models)=BIC_value(K_actuel-exp.cond.Kmin+1,n_models)+2*ENT(exp.data.obs,exp.res.mix.mu{K_actuel,n_models},exp.res.mix.S{K_actuel,n_models},exp.res.mix.p{K_actuel,n_models},K_actuel,exp.cond.n);
                    end
                    break;
                end
          end
          
        K_actuel=K_actuel+1;
        end
        K_actuel=exp.cond.Kmin;
    end
i_actuel=i_actuel+1;
end

%% Build the solutions for the MixCombi script

n_model_BIC=find(min(min(BIC_value))==min(BIC_value));
exp.res.BIC.K=exp.cond.Kmin+find(min(BIC_value(:,n_model_BIC))==BIC_value(:,n_model_BIC))-1;
exp.res.BIC.mu=exp.res.mix.mu{exp.res.BIC.K,n_model_BIC};
exp.res.BIC.S=exp.res.mix.S{exp.res.BIC.K,n_model_BIC};
exp.res.BIC.p=exp.res.mix.p{exp.res.BIC.K,n_model_BIC};

n_model_ICL=find(min(min(ICL_value))==min(ICL_value));
exp.res.ICL.K=exp.cond.Kmin+find(min(ICL_value(:,n_model_ICL))==ICL_value(:,n_model_ICL))-1;
exp.res.ICL.mu=exp.res.mix.mu{exp.res.ICL.K,n_model_ICL};
exp.res.ICL.S=exp.res.mix.S{exp.res.ICL.K,n_model_ICL};
exp.res.ICL.p=exp.res.mix.p{exp.res.ICL.K,n_model_ICL};
