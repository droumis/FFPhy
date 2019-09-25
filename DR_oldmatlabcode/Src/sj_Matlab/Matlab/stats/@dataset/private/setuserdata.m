function a = setuserdata(a,newdata)
%SETUSERDATA Set dataset array UserData property.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:32:07 $

if nargin < 2
    error('stats:dataset:setuserdata:TooFewInputs', ...
          'Requires at least two inputs.');
end

a.props.UserData = newdata;
