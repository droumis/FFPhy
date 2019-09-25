% modelcr1 from es413/512 assig	4
% modified for curvefit (x,xdata)
% used by jntcreg for joint conf region
%function y=modelcr1(t,phi)
function y=modelcr1(phi,t)
y=phi(1).*(1-exp(-phi(2)*t));
return
