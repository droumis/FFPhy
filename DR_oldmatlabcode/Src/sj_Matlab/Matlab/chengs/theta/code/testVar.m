function [p, F]= testVar(x, y)
%function [p, F]= testVar(x, y)
%
% Test whether the two samples have identical variance
% Null Hypothesis: var1= var2
% Alternate Hypothesis: var1 ~= var2
%
% Both x and y have to be vectors.


v1= var(x);
v2= var(y);

if v1>v2
    F= v1/v2;
    nn= length(x);
    nd= length(y);
else
    F= v2/v1;
    nn= length(y);
    nd= length(x);
end

p= 2*(1-fcdf(F, nn-1, nd-1));
