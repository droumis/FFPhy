function a = vertcat(varargin)
%VERTCAT Vertical concatenation for dataset arrays.
%   DS = VERTCAT(DS1, DS2, ...) vertically concatenates the dataset arrays
%   DS1, DS2, ... .  Observation names, when present, must be unique across
%   datasets.  VERTCAT fills in default observation names for the output when
%   some of the inputs have names and some do not.
%
%   Variable names for all dataset arrays must be identical except for order.
%   VERTCATCAT concatenates by matching variable names.
%
%   See also DATASET/CAT, DATASET/HORZCAT.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.3 $  $Date: 2006/12/15 19:31:53 $

b = varargin{1};
if ~isa(b,'dataset')
    error('stats:dataset:vertcat:InvalidInput', ...
          'All input arguments must be datasets.');
end
a = b;
[a_varnames,a_varord] = sort(b.varnames);
for i = 2:nargin
    b = varargin{i};
    if ~isa(b,'dataset')
        error('stats:dataset:vertcat:InvalidInput', ...
              'All input arguments must be datasets.');
    elseif a.nvars ~= b.nvars
        error('stats:dataset:vertcat:SizeMismatch', ...
              'All datasets in the bracketed expression must have the same number of variables.');
    end
    
    [b_varnames,b_varord] = sort(b.varnames);
    if ~all(strcmp(a_varnames,b_varnames))
        error('stats:dataset:vertcat:UnequalVarNames', ...
              'All datasets in the bracketed expression must have the same variable names.');
    end
    
    if ~isempty(a.obsnames) && ~isempty(b.obsnames)
        if checkduplicatenames(b.obsnames,a.obsnames)
            error('stats:dataset:vertcat:DuplicateObsnames', ...
                  'Duplicate observation names.');
        end
        a.obsnames = vertcat(a.obsnames, b.obsnames);
    elseif ~isempty(b.obsnames) % && isempty(a.obsnames)
        a.obsnames = vertcat(strcat({'Obs'},num2str((1:a.nobs)','%d')), b.obsnames);
        a.obsnames = genuniquenames(a.obsnames,a.nobs+1);
    elseif ~isempty(a.obsnames) % && isempty(b.obsnames)
        a.obsnames = vertcat(a.obsnames, strcat({'Obs'},num2str(a.nobs+(1:b.nobs)','%d')));
        a.obsnames = genuniquenames(a.obsnames,a.nobs+1);
    end

    b_reord(a_varord) = b_varord;
    for i = 1:a.nvars
        try
            a.data{i} = vertcat(a.data{i}, b.data{b_reord(i)});
        catch, rethrow(lasterror); end
    end
    a.nobs = a.nobs + b.nobs;
end

