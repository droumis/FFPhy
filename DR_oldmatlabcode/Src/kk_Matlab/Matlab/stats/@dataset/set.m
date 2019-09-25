function b = set(a,varargin)
%SET Set a dataset array property value.
%   B = SET(A,'PropertyName',VALUE) sets the property 'PropertyName' of the
%   dataset A to the value VALUE.
%  
%   B = SET(A,'Property1',Value1,'Property2',Value2,...) sets multiple dataset
%   property values with a single statement.
%
%   SET(A,'Property') displays possible values for the specified property in A.
%
%   SET(A) displays all properties of A and their possible values. 
%
%   See also DATASET/GET, DATASET/SUMMARY.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:31:45 $

if nargin < 3
    propertyNames = [ fieldnames(a.props); {'ObsNames'; 'VarNames'} ];
    
    propertyVals   = cell2struct(cell(size(propertyNames)), propertyNames, 1);
    propertyDescrs = cell2struct(cell(size(propertyNames)), propertyNames, 1);
    
    propertyVals.Description   = {};
    propertyDescrs.Description = 'A dataset''s ''Description'' property does not have a fixed set of values.';
    propertyVals.Units         = {};
    propertyDescrs.Units       = 'Cell array of strings, or an empty cell array.';
    propertyVals.DimNames      = {};
    propertyDescrs.DimNames    = 'Cell array of strings, or an empty cell array.';
    propertyVals.UserData      = {};
    propertyDescrs.UserData    = 'A dataset''s ''UserData'' property does not have a fixed set of values.';
    propertyVals.ObsNames      = {};
    propertyDescrs.ObsNames    = 'Cell array of non-empty, unique strings, or an empty cell array.';
    propertyVals.VarNames      = {};
    propertyDescrs.VarNames    = 'Cell array of non-empty, unique strings.';
    
    if nargin == 2
        name = matchpropertyname(a,varargin{1});
        if nargout == 1
            b = propertyVals.(name);
        else
            disp(propertyDescrs.(name));
        end        
    else
        if nargout == 1
            b = propertyVals;
        else
            disp(propertyDescrs);
        end
    end
    
elseif mod(nargin,2) == 1 %name/value pairs
    b = a;
    for i = 1:(nargin-1)/2
        name = varargin{2*i-1};
        p = varargin{2*i};
        b = setproperty(b,name,p);
    end
else
    error('stats:dataset:set:WrongNumberArgs', ...
          'Wrong number of arguments.');
end

