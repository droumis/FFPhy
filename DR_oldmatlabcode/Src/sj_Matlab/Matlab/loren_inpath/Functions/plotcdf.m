% [h] = plotcdf(data, linestyle)
%       Plots a cdf of data with the line style specified 
function [h] = plotcdf(data, linestyle);

[f x] = ecdf(data);
if (nargin > 1)
    h = stairs(x,f, linestyle);
else
    h = stairs(x,f);
end
