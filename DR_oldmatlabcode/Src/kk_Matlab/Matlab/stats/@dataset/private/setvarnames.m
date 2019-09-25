function a = setvarnames(a,newnames,vars)
%SETVARNAMES Set dataset array variable names.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.3 $  $Date: 2006/12/15 19:32:08 $

if nargin < 2
    error('stats:dataset:setvarnames:TooFewInputs', ...
          'Requires at least two inputs.');
end

if nargin == 2
    if isstring(newnames)
        if a.nvars ~= 1
            error('stats:dataset:setvarnames:IncorrectNumberOfVarnames', ...
                  'NEWNAMES must contain one name for each variable in A.');
        end
        newnames = cellstr(newnames);
    elseif iscell(newnames)
        if numel(newnames) ~= a.nvars
            error('stats:dataset:setvarnames:IncorrectNumberOfVarnames', ...
                  'NEWNAMES must have one name for each variable in A.');
        elseif ~all(cellfun(@isstring,newnames)) % require a nonempty row of chars
            error('stats:dataset:setvarnames:InvalidVarnames', ...
                  'NEWNAMES must be a nonempty string or a cell array of nonempty strings.');
        elseif checkduplicatenames(newnames);
            error('stats:dataset:setvarnames:DuplicateVarnames', ...
                  'Duplicate variable names.');
        end
    elseif ~iscell(newnames)
        error('stats:dataset:setvarnames:InvalidVarnames', ...
              'NEWNAMES must be a nonempty string or a cell array of nonempty strings.');
    end
    a.varnames = strtrim(newnames(:))'; % this conveniently converts {} to a 1x0
    
else % if nargin == 3
    varIndices = getvarindices(a,vars);
    if isstring(newnames)
        if ~isscalar(varIndices)
            error('stats:dataset:setvarnames:IncorrectNumberOfVarnames', ...
                  'NEWNAMES must contain one name for each variable name being replaced.');
        end
        newnames = cellstr(newnames);
    elseif iscell(newnames)
        if length(newnames) ~= length(varIndices)
            error('stats:dataset:setvarnames:IncorrectNumberOfVarnames', ...
                  'NEWNAMES must contain one name for each variable name being replaced.');
        elseif ~all(cellfun(@isstring,newnames)) % require a nonempty row of chars
            error('stats:dataset:setvarnames:InvalidVarnames', ...
                  'NEWNAMES must be a nonempty string or a cell array of nonempty strings.');
        end
    else
        error('stats:dataset:setvarnames:InvalidVarnames', ...
              'NEWNAMES must be a nonempty string or a cell array of nonempty strings.');
    end
    if checkduplicatenames(newnames,a.varnames,varIndices);
        error('stats:dataset:setvarnames:DuplicateVarnames', ...
              'Duplicate variable names.');
    end
    a.varnames(varIndices) = strtrim(newnames);
end

checkreservednames(a.varnames);
[a.varnames,wereModified] = genvalidnames(a.varnames);
if wereModified
    warning('stats:dataset:setvarnames:ModifiedVarnames', ...
            'Variable names were modified to make them valid MATLAB identifiers.');
end


function tf = isstring(s) % require a nonempty row of chars
tf = ischar(s) && isvector(s) && (size(s,1) == 1) && ~isempty(s);
