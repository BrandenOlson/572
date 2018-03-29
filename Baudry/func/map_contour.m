function []=map_contour(y,mu,S,p,K)

%figure;

d=size(y,2);

if d==1

    n=size(y,1);

    color=['r';'g';'b';'c';'m';'y';'k';'g';'r';'b';'c';'m';'y';'k'];

    if K==1

        hold on;

        k=1;

        xlabel('X_1');
        ylabel('X_2');

        x=min(y(:,1))-2:0.01:max(y(:,1))+2;
        yd=p(k)*normpdf(x,mu(k,:),S(:,:,k));

        m=max(yd);

        plot(y(:,1),repmat(m/10,n,1),[color(k) '.']);

        xlim([min(y(:,1))-2 max(y(:,1))+2]);
        ylim([0 m*1.1]);

        plot(x,yd,[color(k) '-']);

        hold off;

    else

        for k=1:K
            posterior(:,k)=p(k)*mvnpdf(y,mu(k,:),S(:,:,k));
        end

        Z=(repmat(max(posterior,[],2),1,K)==posterior);

        hold on;

        xlabel('X_1');
        ylabel('X_2');

        m=0;

        x=min(y(:,1))-2:0.01:max(y(:,1))+2;

        for k=1:K
            yd=p(k)*normpdf(x,mu(k,:),S(:,:,k));
            plot(x,yd,[color(k) '-']);
            m=max(max(yd*1.1),m);
        end

        xlim([min(y(:,1))-2 max(y(:,1))+2]);
        ylim([0 m*1.1]);

        for k=1:K
            yl=y(find(Z(:,k)==1),:);
            plot(yl(:,1),repmat(m/10,size(yl,1),1),[color(k) '.']);
        end

        hold off;

    end

elseif d==2
    
    for k=1:K
        posterior(:,k)=p(k)*mvnpdf(y,mu(k,:),S(:,:,k));
    end

    Z=(repmat(max(posterior,[],2),1,K)==posterior);

    color=['r';'g';'b';'c';'m';'y';'k';'r';'g';'b';'c';'m';'y';'k';'r';'g';'b';'c';'m';'y';'k';'r';'g';'b';'c';'m';'y';'k'];
%    color=['k';'k';'k';'k';'k';'k';'k';'k';'k';'k';'k';'k';'k';'k'];
    
    hold on;

    for k=1:K
        yl=y(find(Z(:,k)==1),:);
        plot(yl(:,1),yl(:,2),[color(k) '.']);
    end

    xlabel('X_1');
    ylabel('X_2');

    levels=[0.01 0.01];

    x_lim=[min(y(:,1))-2 max(y(:,1))+2];
    y_lim=[min(y(:,2))-2 max(y(:,2))+2];

         for k=1:K
             mvncontour(x_lim,y_lim,mu(k,:),S(:,:,k),levels,p(k),color(k));
         end

    hold off;

elseif d==3
    
    for k=1:K
        posterior(:,k)=p(k)*mvnpdf(y,mu(k,:),S(:,:,k));
    end

    Z=(repmat(max(posterior,[],2),1,K)==posterior);

    color=['r';'g';'b';'c';'m';'y';'k';'g';'r';'g';'b';'c';'m';'y';'k';'r';'g';'b';'c';'m';'y';'k';'r';'g';'b';'c';'m';'y';'k'];

    hold on;

    for k=1:K
        yl=y(find(Z(:,k)==1),:);
        plot3(yl(:,1),yl(:,2),yl(:,3),[color(k) '.']);
    end

    xlabel('X_1');
    ylabel('X_2');
    zlabel('X_3');

%     levels=[0.005 0.005];

%     x_lim=[min(y(:,1))-2 max(y(:,1))+2];
%     y_lim=[min(y(:,2))-2 max(y(:,2))+2];
%     z_lim=[min(y(:,3))-2 max(y(:,3))+2];
    
    %     for k=1:K
    %         mvncontour(x_lim,y_lim,mu(k,:),S(:,:,k),levels,p(k),color(k));
    %     end

    hold off;

elseif d==4
    
    for k=1:K
        posterior(:,k)=p(k)*mvnpdf(y,mu(k,:),S(:,:,k));
    end

    Z=(repmat(max(posterior,[],2),1,K)==posterior);

    color=['r';'g';'b';'c';'m';'y';'k';'r';'g';'b';'c';'m';'y';'k';'r';'g';'b';'c';'m';'y';'k';'r';'g';'b';'c';'m';'y';'k'];
%    color=['k';'k';'k';'k';'k';'k';'k';'k';'k';'k';'k';'k';'k';'k'];
    
    hold on;

    for k=1:K
        yl=y(find(Z(:,k)==1),1:2);
        plot(yl(:,1),yl(:,2),[color(k) '.']);
    end

    xlabel('X_1');
    ylabel('X_2');

    levels=[0.01 0.01];

    x_lim=[min(y(:,1))-2 max(y(:,1))+2];
    y_lim=[min(y(:,2))-2 max(y(:,2))+2];

    hold off;


end