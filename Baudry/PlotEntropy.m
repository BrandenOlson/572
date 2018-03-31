addpath('func');

%% Entropy plots

% Build the entropy vector
E=[exp.res];
E=[E.combi];
E=[E.ent];
E=[E,exp.res.BIC.ent];

% Plot the simple entropy plot and the one-change-point regression
hold on;
figure('visible', 'off');
pcws_reg(1:exp.res.BIC.K,E);
xlabel('K');
ylabel('Entropy');
figure_name = strcat(filename, "_entropy.pdf");
print("-dpdf", figure_name);
movefile(figure_name, strcat("../BaudryFigures/", filename));

% Plot the simple "difference in entropy" plot
hold on;
figure;
plot(2:exp.res.BIC.K,E(2:end)-E(1:(end-1)),'ok');
xlabel('K');
ylabel('Difference in entropy');

figure_name = strcat(filename, "_diff_entropy.pdf");
print("-dpdf", figure_name);
movefile(figure_name, strcat("../BaudryFigures/", filename));

% Compute the number of observations merged at each step
clear Nbr;
for K=1:(exp.res.BIC.K-1)
    i=find(sum(exp.res.combi(K).M,2)==2);
    Nbr(K)=sum(exp.res.combi(K+1).labels*exp.res.combi(K).M(i,:)');
end

% Plot the rescaled entropy plot and the one-change-point regression
hold on;
figure;
pcws_reg(cumsum([0,Nbr]),E);
xlabel('Cumul. count of merged obs.');
ylabel('Entropy');
figure_name = strcat(filename, "_rescaled_entropy.pdf");
print("-dpdf", figure_name);
movefile(figure_name, strcat("../BaudryFigures/", filename));

% Plot the "normalized difference in entropy" plot
hold on;
figure;
plot(2:exp.res.BIC.K,(E(2:end)-E(1:(end-1)))./Nbr,'ok');
xlabel('K');
ylabel('Normalized difference in entropy');

figure_name = strcat(filename, "_normalized_diff_entropy.pdf");
print("-dpdf", figure_name);
movefile(figure_name, strcat("../BaudryFigures/", filename));
