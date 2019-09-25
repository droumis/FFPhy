function []=normqq(y)
%NORMQQ	NORMQQ(Y) produces a normal probability plot of the elements of Y.

% GKS  27 July 95

y = sort(y(:));
nd = normq( (1:length(y))./(length(y)+1) );
plot(nd,nd,nd,y,'x');
xlabel('Normal Deviates');
ylabel('Sample Deviates');
title('Normal Probability Plot');
