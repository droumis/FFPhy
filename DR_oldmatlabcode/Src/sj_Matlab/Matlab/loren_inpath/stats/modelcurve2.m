function y=modelcurve2(p,t)
y=exp(-p(1).*t(:,1)).*sin(2*pi*p(2).*t(:,2));
return