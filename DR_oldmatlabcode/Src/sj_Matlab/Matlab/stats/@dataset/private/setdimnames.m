function a = setdimnames(a,newnames)
%SETDIMNAMES Set dataset array DimNames property.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.3 $  $Date: 2006/12/15 19:32:03 $

if nargin < 2
    error('stats:dataset:setdimnames:TooFewInputs', ...
          'Requires at least two inputs.');
end

if nargin == 2
    if isempty(newnames)
        a.props.DimNames = {}; % do this for cosmetics
        return
    end
    if numel(newnames) ~= a.ndims
        error('stats:dataset:setdimnames:WrongLength', ...
              'NEWNAMES must have one element for each dimension in A.');
    elseif ~iscell(newnames)
        error('stats:dataset:setdimnames:InvalidDimnames', ...
              'NEWNAMES must be a cell array of nonempty strings.');
    elseif ~all(cellfun(@isstring,newnames))
        error('stats:dataset:setdimnames:InvalidDimnames', ...
              'NEWNAMES must contain nonempty strings.');
    elseif checkduplicatenames(newnames);
        error('stats:dataset:setdimnames:DuplicateDimnames', ...
              'Duplicate dimension names.');
    end
    a.props.DimNames = newnames(:);
    
end

function tf = isstring(s) % require a nonempty row of chars
tf = ischar(s) && isvector(s) && (size(s,1) == 1) && ~isempty(s);
