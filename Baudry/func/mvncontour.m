function []=mvncontour(x_lim,y_lim,mu,S,levels,p,LineSpec)
%MVNCONTOUR   Trace les isodensités d'une gaussienne multivariée.
%   MVNCONTOUR(x_lim,y_lim,mu,S) trace un niveau d'isodensité de la gaussienne de
%   moyenne mu et de matrice de covariance S dans un graphe aux limites
%   données par [min(x_lim):max(x_lim)] et [min(y_lim):max(y_lim)].
%
%   MVNCONTOUR(x_lim_y_lim,mu,S,p,levels) fait de même mais pour les niveaux
%   donnés par levels (si levels=N, mvncontour trace N courbes de niveau; si
%   levels=V, mvncontour trace length(V) courbes de niveau, aux niveaux donnés
%   par les valeurs du vecteur V; pour tracer le niveau v, utiliser levels=[v
%   v]).
%
%   MVNCONTOUR(x_lim,y_lim,mu,S,levels,p) fait de même mais en multipliant la
%   densité de la gaussienne considérée par p. Utile si on veut tracer
%   l'isodensité d'un mélange...
%
%   MVNCONTOUR(x_lim,y_lim,mu,S,levels,p,LineSpec) fait de même en
%   utilisant les spécifications de forme et couleur des lignes et
%   marqueurs...

if nargin==4 
    levels=1;
    p=1;
    LineSpec='';
elseif nargin==5
    p=1;
    LineSpec='';
elseif nargin==6
    LineSpec='';
end

X=(min(x_lim):(max(x_lim)-min(x_lim))/100:max(x_lim))';
Y=(min(y_lim):(max(y_lim)-min(y_lim))/100:max(y_lim))';

Z=zeros(101,101);

for i=1:101
    Z(:,i)=p*mvnpdf([X,repmat(Y(i),101,1)],mu,S);
end

contour(X,Y,Z',levels,LineSpec);
%colormap(cool);
