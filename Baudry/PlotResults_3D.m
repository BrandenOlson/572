addpath('func');

%% 3D  Plot the observed data

hold on;
figure;

plot3(exp.data.obs(:,1),exp.data.obs(:,2),exp.data.obs(:,3),'*b');
grid on;

title(['Observed Data. n=' int2str(exp.cond.n)]);

figure_name = strcat(filename, "_observed.pdf");
print("-dpdf", figure_name);
movefile(figure_name, strcat("../BaudryFigures/", filename));

%% 3D Plot the BIC solution

hold on;
figure;

map_contour(exp.data.obs,exp.res.BIC.mu,exp.res.BIC.S,exp.res.BIC.p,exp.res.BIC.K);

title(['BIC solution. K=' int2str(exp.res.BIC.K) '. ENT=' int2str(exp.res.BIC.ent)]);
grid on;
view(3);
figure_name = strcat(filename, "_BIC.pdf");
print("-dpdf", figure_name);
movefile(figure_name, strcat("../BaudryFigures/", filename));

%% 3D Plot the ICL solution (optional)

try % In case the user did not provide the ICL solution
    hold on;
    figure;

    map_contour(exp.data.obs,exp.res.ICL.mu,exp.res.ICL.S,exp.res.ICL.p,exp.res.ICL.K);

    title(['ICL solution. K=' int2str(exp.res.ICL.K) '. ENT=' int2str(ENT(exp.data.obs,exp.res.ICL.mu,exp.res.ICL.S,exp.res.ICL.p,exp.res.ICL.K,exp.cond.n))]);
    grid on;
    view(3);
    figure_name = strcat(filename, "_ICL.pdf");
    print("-dpdf", figure_name);
    movefile(figure_name, strcat("../BaudryFigures/", filename));
catch
end

%% 3D Plot the combined solutions

x_lim=[min(exp.data.obs(:,1))-2 max(exp.data.obs(:,1))+2];
y_lim=[min(exp.data.obs(:,2))-2 max(exp.data.obs(:,2))+2];

X=(min(x_lim):(max(x_lim)-min(x_lim))/100:max(x_lim))';
Y=(min(y_lim):(max(y_lim)-min(y_lim))/100:max(y_lim))';

color=['r';'g';'b';'c';'m';'y';'k';'g';'r';'b';'c';'m';'y';'k'];
col=color(1:exp.res.BIC.K);

for Kb=2:exp.res.BIC.K-1
    K=exp.res.BIC.K+1-Kb;

    for k=1:K
       coln(k)=col(min(find(exp.res.combi(K).M(k,:)==1)));
    end
    
    figure;
    
    hold on;

    xlabel('X_1');
    ylabel('X_2');
    zlabel('X_3');

%    levels=[0.005 0.005 0.005];

    M=eye(K,K);

    for KM=K:(exp.res.BIC.K-1)
        
        M=M*exp.res.combi(KM).M;

    end
    
    for k=1:K
        
        yl=exp.data.obs(find(exp.res.combi(K).labels(:,k)==1),:);
        plot3(yl(:,1),yl(:,2),yl(:,3),[col(k) '.']);

    end

    title(['Combined solution. K=' int2str(K) '. ENT=' int2str(exp.res.combi(K).ent)]);
    grid on;
    view(3);

    figure_name = strcat(filename, "_K_", int2str(K), ".pdf");
    print("-dpdf", figure_name);
    movefile(figure_name, strcat("../BaudryFigures/", filename));


%    file=['C:\Maths\Article avec Celeux et Raftery\Graphics\Gottardo_2D\Gottardo_2D_C_combi(' int2str(K) ').eps']
    
%    print('-depsc',file)

    hold off;

    % % % waitforbuttonpress;
    
end

