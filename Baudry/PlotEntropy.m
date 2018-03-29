addpath('func');

%% Entropy plots

% Build the entropy vector
E=[exp.res];
E=[E.combi];
E=[E.ent];
E=[E,exp.res.BIC.ent];

% Plot the simple entropy plot and the one-change-point regression
figure;
pcws_reg(1:exp.res.BIC.K,E);
xlabel('K');
ylabel('Entropy');

% Plot the simple "difference in entropy" plot
figure;
plot(2:exp.res.BIC.K,E(2:end)-E(1:(end-1)),'ok');
xlabel('K');
ylabel('Difference in entropy');

% Compute the number of observations merged at each step
clear Nbr;
for K=1:(exp.res.BIC.K-1)
    i=find(sum(exp.res.combi(K).M,2)==2);
    Nbr(K)=sum(exp.res.combi(K+1).labels*exp.res.combi(K).M(i,:)');
end

% Plot the rescaled entropy plot and the one-change-point regression
figure;
pcws_reg(cumsum([0,Nbr]),E);
xlabel('Cumul. count of merged obs.');
ylabel('Entropy');

% Plot the "normalized difference in entropy" plot
figure;
plot(2:exp.res.BIC.K,(E(2:end)-E(1:(end-1)))./Nbr,'ok');
xlabel('K');
ylabel('Normalized difference in entropy');
