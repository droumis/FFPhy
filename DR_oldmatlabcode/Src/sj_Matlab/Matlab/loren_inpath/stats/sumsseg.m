function q=sumss(x,y,b2,b)
%function q=sumss(x,y,b2,b)
% calculates sum of squares for regression model ym
% e.g. for ym=p0 + p1.x1 + p2.x2 +p3.x3
% x=[1 x1 x2 x3 ] for n rows n= # data
% if be are plotting joint CR for p1 and p2 then b2=[p1,p2]
% and b=rest of param Least Square estimates b=[p0 p3]
% ym= b(1) + b2(1).*x(:,1) +b2(2).*x(:,2) + b(2).*x(:,3)
%
phi=b2;t=x;
% b (LS estimates of remaining params) not used here- only 2 params
ym=phi(1).*(1-exp(-phi(2)*t));
q=sum((ym-y).^2);
return
