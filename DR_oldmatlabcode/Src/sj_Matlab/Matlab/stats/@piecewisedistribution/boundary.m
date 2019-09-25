function [p,q] = boundary(obj,j)
%BOUNDARY Boundary between segments of piecewise distribution.
%    [P,Q]=BOUNDARY(OBJ) returns the boundary points between the segments
%    of the piecewise distribution defined by OBJ.  P is a vector of the
%    probability values at each boundary.  Q is a vector of the quantile
%    values at each boundary.
%
%    [P,Q]=BOUNDARY(OBJ,J) returns P and Q for the Jth boundary.
%
%    See also PARETOTAILS, PIECEWISEDISTRIBUTION/NSEGMENTS.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:57:17 $

if nargin==1
    p = obj.P;
    q = obj.Q;
else
    if ~isvector(j)
        j = j(:);
    end
    if any(~ismember(j,1:numel(obj.P)))
        error('stats:piecewisedistribution:boundary:BadIndex',...
              'Input J must be a valid boundary number.');
    end
    p = obj.P(j);
    q = obj.Q(j);
end
