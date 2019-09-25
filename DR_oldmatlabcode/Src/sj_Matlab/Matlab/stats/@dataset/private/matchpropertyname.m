function name = matchpropertyname(a,name)
%MATCHPROPERTYNAME Validate a dataset array property name.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:32:01 $

% This matches names against this list of property names, including
% 'ObsNames'and 'VarNames', even though they are not under the 'props' field.
propertyNames = [ fieldnames(a.props); {'ObsNames'; 'VarNames'} ];

if ~(ischar(name) && isvector(name) && (size(name,1)==1))
    error('stats:dataset:matchpropertyname:InvalidPropertyName', ...
          'Invalid property name.');
end

j = find(strncmp(name,propertyNames,length(name)));
if isempty(j)
    error('stats:dataset:matchpropertyname:UnknownProperty', ...
          'Unknown dataset property: %s.', name);
elseif ~isscalar(j)
    error('stats:dataset:matchpropertyname:AmbiguousProperty', ...
          'Ambiguous dataset property name: %s.', name);
end

name = propertyNames{j};
