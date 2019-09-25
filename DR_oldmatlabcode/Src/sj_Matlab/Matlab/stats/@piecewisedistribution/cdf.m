function p=cdf(pd,x)
%CDF Cumulative distribution function for piecewise distribution.
%    P=CDF(OBJ,X) returns an array P of values of the cumulative distribution
%    function (CDF) for the piecewise distribution object OBJ, evaluated
%    at the values in the array X.
%
%    See also PARETOTAILS, PIECEWISEDISTRIBUTION/ICDF.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:57:18 $

% Determine the segment that each point occupies
s = segment(pd,x);

% Invoke the appropriate cdf for each segment
p = nan(size(x),class(x));
for j=1:max(s(:))
    t = (s==j);
    if any(t(:))
        p(t) = pd.distribution(j).cdf(x(t));
    end
end
