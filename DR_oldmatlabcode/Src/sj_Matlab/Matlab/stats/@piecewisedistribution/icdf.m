function x=icdf(pd,p)
%ICDF Inverse cumulative distribution function for piecewise distribution.
%    X=ICDF(OBJ,P) returns an array X of values of the inverse cumulative
%    distribution function (ICDF) for the piecewise distribution object OBJ,
%    evaluated at the values in the array P.
%
%    See also PARETOTAILS, PIECEWISEDISTRIBUTION/CDF.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:57:21 $

% Determine the segment that each point occupies
s = segment(pd,[],p);

% Invoke the appropriate icdf for each segment
x = NaN(size(p),class(p));
for j=1:max(s(:))
    t = (s==j);
    if any(t(:))
        x(t) = pd.distribution(j).icdf(p(t));
    end
end
