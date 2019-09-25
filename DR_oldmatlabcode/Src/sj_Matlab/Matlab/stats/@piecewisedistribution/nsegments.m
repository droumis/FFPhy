function n = nsegments(obj)
%NSEGMENTS Number of segments defined in piecewise distribution.
%    N=NSEGMENTS(OBJ) returns the number of segments in the piecewise
%    distribution defined by OBJ.
%
%    See also PARETOTAILS, PIECEWISEDISTRIBUTION/BOUNDARY, PIECEWISEDISTRIBUTION/SEGMENT.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:57:22 $

n = numel(obj.P)+1;
        
