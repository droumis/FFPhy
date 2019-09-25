function a = horzcat(varargin)
%HORZCAT Horizontal concatenation for dataset arrays.
%   DS = HORZCAT(DS1, DS2, ...) horizontally concatenates the dataset arrays
%   DS1, DS2, ... .  You may concatenate dataset arrays that have duplicate
%   variable names, however, the variables must contain identical data, and
%   HORZCAT includes only one copy of the variable in the output dataset.
%
%   Observation names for all dataset arrays that have them must be identical
%   except for order. HORZCAT concatenates by matching observation names when
%   present, or by position for datasets that do not have observation names.
%
%   See also DATASET/CAT, DATASET/VERTCAT.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.3 $  $Date: 2006/12/15 19:31:39 $

b = varargin{1};
if ~isa(b,'dataset')
    error('stats:dataset:horzcat:InvalidInput', ...
          'All input arguments must be datasets.');
end
a = b;
if ~isempty(a.obsnames)
    [a_sortobsnames,a_obsord] = sort(a.obsnames);
end
for i = 2:nargin
    b = varargin{i};
    if ~isa(b,'dataset')
        error('stats:dataset:horzcat:InvalidInput', ...
              'All input arguments must be datasets.');
    elseif a.nobs ~= b.nobs
        error('stats:dataset:horzcat:SizeMismatch', ...
              'All datasets in the bracketed expression must have the same number of observations.');
    end
    [dups,dupBVars] = checkduplicatenames(b.varnames,a.varnames);
    if dups
        for i = find(dupBVars)
            if ~isequal(a.data{getvarindices(a,b.varnames(i))},b.data{i})
                error('stats:dataset:horzcat:DuplicateVarnames', ...
                      'Duplicate variable names with distinct data.');
            end
        end
        uniqBVars = find(~dupBVars);
    else
        uniqBVars = 1:b.nvars;
    end
    
    if ~isempty(a.obsnames) && ~isempty(b.obsnames)
        [b_sortobsnames,b_obsord] = sort(b.obsnames);
        if ~all(strcmp(a_sortobsnames,b_sortobsnames))
            error('stats:dataset:horzcat:UnequalObsnames', ...
                  'All datasets in the bracketed expression must have the same observation names.');
        end
        b_reord(a_obsord) = b_obsord;
        a.data = horzcat(a.data, cell(1,length(uniqBVars)));
        for i = 1:length(uniqBVars)
            bVar = b.data{uniqBVars(i)};
            sizeOut = size(bVar);
            a.data{a.nvars+i} = reshape(bVar(b_reord,:),sizeOut);
        end
    else
        if isempty(a.obsnames) && ~isempty(b.obsnames)
            a.obsnames = b.obsnames;
            [a_sortobsnames,a_obsord] = sort(a.obsnames);
        end
        a.data = horzcat(a.data, b.data(uniqBVars));
    end
    % We've already weeded out duplicate var names, either by erroring or by
    % ignoring the second copy, so there's no need to uniqueify var names.
    a.nvars = a.nvars + length(uniqBVars);
    a.varnames = horzcat(a.varnames, b.varnames(uniqBVars));
end
