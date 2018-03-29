function [cf,a1f,b1f,a2f,b2f,ss]=pcws_reg(x,y)

% pcws_reg computes the piecewise linear regression ---with two
% pieces--- to (x,y), for any possible change point, chooses the one leading
% to the smallest least-square error, and returns and plots the
% corresponding parameters.

a1=inf(1,size(x,2)-1);
a2=inf(1,size(x,2)-1);
b1=inf(1,size(x,2)-1);
b2=inf(1,size(x,2)-1);
ss=inf(1,size(x,2)-1);

for c=2:(size(x,2)-1)
    x1=x(1:c);
    y1=y(1:c);
    x2=x(c:end);
    y2=y(c:end);
    
    a1(c)=(sum(x1.*y1)-sum(x1)*mean(y1))/(sum(x1.^2)-sum(x1)^2/c);
    b1(c)=-a1(c)*mean(x1)+mean(y1);
    
    a2(c)=(sum(x2.*y2)-mean(x2)*sum(y2))/(sum(x2.^2)-sum(x2)^2/size(x2,2));
    b2(c)=-a2(c)*mean(x2)+mean(y2);

    ss(c)=sum((a1(c)*x1+b1(c)-y1).^2) + sum((a2(c)*x2+b2(c)-y2).^2);
end

cf=find(ss==min(ss(2:end)));
a1f=a1(cf);
a2f=a2(cf);
b1f=b1(cf);
b2f=b2(cf);

plot(x,y,'ko');
hold on;
plot(x(1:cf),a1f*x(1:cf)+b1f,'r--');
plot(x(cf:end),a2f*x(cf:end)+b2f,'r--');
hold off;
