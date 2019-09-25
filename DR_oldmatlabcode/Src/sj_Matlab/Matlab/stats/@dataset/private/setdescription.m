function a = setdescription(a,newdescr)
%SETUNITS Set dataset array Description property.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:32:02 $

if nargin < 2
    error('stats:dataset:setdescription:TooFewInputs', ...
          'Requires at least two inputs.');
end

if nargin == 2
    if isempty(newdescr)
        a.props.Description = {};
        return
    end
    a.props.Description = newdescr;
end
