function [c,loc] = join(a,b,varargin)
%JOIN Merge observations from two dataset arrays.
%   C = JOIN(A, B) creates a dataset array C by merging observations from the
%   two dataset arrays A and B.  JOIN performs this merge by first finding the
%   key variables, i.e., a pair of dataset variables, one in A and one in B,
%   that share the same name.  B's key must contain unique values, and must
%   contain all the values that are present in A's key.  JOIN then uses these
%   two key variables to define a many-to-one correspondence between
%   observations in A and those in B.  JOIN uses this correspondence to
%   replicate B's observations and combine them with A's observations to
%   create C.
%
%   C contains one observation for each observation in A.  Variables in C
%   include all of the variables from A, as well as one variable corresponding
%   to each variable in B (except for B's key).
%
%   C = JOIN(A, B, KEY) performs the join using the variable specified by KEY
%   as the key variable in both A and B.  KEY is a positive integer, a variable
%   name, a cell array containing a variable name, or a logical vector with
%   one true entry.
%
%   C = JOIN(A, B, 'PARAM1',val1, 'PARAM2',val2, ...) allows you to specify
%   optional parameter name/value pairs to control how the dataset variables
%   in A and B are used in the join.  Parameters are:
%
%      'Key'       - specifies the variable to use as a key in both A and B.
%      'LeftKey'   - specifies the variable to use as a key in A.
%      'RightKey'  - specifies the variable to use as a key in B.
%
%   You may provide either the 'Key' parameter, or both the 'LeftKey' and
%   'RightKey' parameters.  The value for these parameters is a positive
%   integer, a variable name, a cell array containing a variable name, or a
%   logical vector with one true entry.
%
%      'LeftVars'  - specifies the variables from A to include in C.  By
%                    default, JOIN includes all variables from A.
%      'RightVars' - specifies the variables from B to include in C.  By
%                    default, JOIN includes all variables from B except the
%                    key variable.
%
%   The value for these parameters is a positive integer, a vector of positive
%   integers, a variable name, a cell array containing one or more variable
%   names, or a logical vector.
%
%   [C,LOC] = JOIN(...) returns an index vector LOC, where the observations in
%   C are constructed by horizontally concatenating A(:,LEFTVARS) and
%   B(LOC,RIGHTVARS).
%
%   See also DATASET/HORZCAT, DATASET/SORTROWS, DATASET/UNIQUE.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.6.2 $  $Date: 2006/12/15 19:31:41 $

if nargin < 2
    error('stats:dataset:join:TooFewInputs', ...
          'Requires at least two inputs.');
elseif ~isa(a,'dataset') || ~isa(b,'dataset')
    error('stats:dataset:horzcat:InvalidInput', ...
          'A and B must be dataset arrays.');
end

if nargin < 4
    if nargin == 2 % join(a,b), use the variable that they have in common
        key = find(ismember(a.varnames,b.varnames));
        if isempty(key)
            error('stats:dataset:join:CantInferKey', ...
                  'Cannot find a common dataset variable to use as a key variable.');
        end
        key = a.varnames(key);
    elseif nargin == 3 % join(a,b,key), use the same variables on the left and the right
        key = varargin{1};
    end
    leftkey = key;
    rightkey = key;
    leftvars = [];
    rightvars = [];
    
else % join(a,b,'key',keyvar,...) or join(a,b,'leftkey',leftkeyvar,'rightkey',rightkeyvar,...)
    pnames = {'key' 'leftkey' 'rightkey' 'leftvars' 'rightvars'};
    dflts =  {   []        []         []         []          []};
    [eid,errmsg,key,leftkey,rightkey,leftvars,rightvars] ...
                       = dfswitchyard('statgetargs', pnames, dflts, varargin{:});
    if ~isempty(eid)
        error(sprintf('stats:dataset:join:%s',eid),errmsg);
    end

    if isempty(key)
        if isempty(leftkey) && isempty(rightkey)
            key = find(ismember(a.varnames,b.varnames));
            if isempty(key)
                error('stats:dataset:join:CantInferKey', ...
                      'Cannot find a common dataset variable to use as a key variable.');
            end
            key = a.varnames(key);
            leftkey = key;
            rightkey = key;
        elseif isempty(leftkey) || isempty(rightkey)
            error('stats:dataset:join:MissingKeyVar', ...
                  'Must provide either the ''Key'' parameter, or both the ''LeftKey'' and ''RightKey'' parameters.');
        end
    else
        if ~isempty(leftkey) || ~isempty(rightkey)
            error('stats:dataset:join:ConflictingInputs', ...
                  'Must provide either the ''Key'' parameter, or both the ''LeftKey'' and ''RightKey'' parameters.');
        end
        leftkey = key;
        rightkey = key;
    end
end

% Make sure the keys exist, do not allow multiple keys.
leftkey = getvarindices(a,leftkey);
rightkey = getvarindices(b,rightkey);
if ~isscalar(leftkey) || ~isscalar(rightkey)
    error('stats:dataset:join:MultipleKeyVars', ...
          'Multiple key variables are not allowed for a dataset join.');
end

% Use all vars from A and B by default, or use the specified vars.
if isempty(leftvars)
    leftvars = 1:a.nvars;
else
    leftvars = getvarindices(a,leftvars);
    if length(unique(leftvars)) < length(leftvars)
        error('stats:dataset:join:DuplicateVars', ...
              'Cannot include a variable twice.');
    end
end
if isempty(rightvars)
    % Leave out B's key var
    rightvars = [1:(rightkey-1) (rightkey+1):b.nvars];
else
    rightvars = getvarindices(b,rightvars);
    if length(unique(rightvars)) < length(rightvars)
        error('stats:dataset:join:DuplicateVars', ...
              'Cannot include a variable twice.');
    end
end
if any(ismember(b.varnames(rightvars),a.varnames(leftvars)))
    error('stats:dataset:join:DuplicateVars', ...
          'Duplicate variable names.');
end

% Get the key var values, and check that they are scalar-valued or
% char-valued.
leftkey = a.data{leftkey};
rightkey = b.data{rightkey};
if (ndims(leftkey) > 2) || (ndims(rightkey) > 2)
    error('stats:dataset:join:NDKeyVar', ...
          'An N-D variable cannot be used as the key for a join.');
end

% Check that B's key contains no duplicates.
if size(rightkey,2) > 1
    if isnumeric(rightkey) || islogical(rightkey) || ischar(rightkey)
        nunique = length(unique(rightkey,'rows'));
    else
        error('stats:dataset:join:MulticolumnKeyVar', ...
              'A %s variable cannot be used as the key for a join if it has multiple columns.', ...
              class(rightkey));
    end
else
    nunique = length(unique(rightkey));
end
if nunique < size(rightkey,1)
    error('stats:dataset:join:DuplicateRightKeyVarValues', ...
          'The key variable for B must have unique values.');
end

% Do the ismember by rows if either key has multiple columns, except if either
% one is a cell array -- the cell method does not accept a 'rows' flag.
byrows = ( (size(leftkey,2) > 1) || (size(rightkey,2) > 1) ) ...
                                 && ~(iscell(leftkey) || iscell(rightkey));

% Use the key vars to find indices from A into B, and make sure every
% observation in A has a corresponding one in B.
try
    if byrows
        [tf,loc] = ismember(leftkey,rightkey,'rows');
    else
        [tf,loc] = ismember(leftkey,rightkey);
    end
catch
    error('stats:dataset:join:KeyIsmemberMethodFailed', ...
          ['Unable to find observations in B corresponding to those in A using the\n', ...
           'key variables.  The ISMEMBER method generated the following error:\n\n%s'], lasterr);
    end
if ~isequal(size(tf),[a.nobs,1])
    error('stats:dataset:join:KeyIsmemberMethodFailed', ...
          ['Unable to find observations in B corresponding to those in A using the\n', ...
           'key variables.  The value returned by the ISMEMBER method had the wrong\n' ...
           'number of rows.']);
elseif any(~tf)
    error('stats:dataset:join:LeftKeyValueNotFound', ...
          'The key variable for B must have contain all values in the key variable for A.');
end

% Create a new dataset by combining the specified variables from A with those
% from B, the latter broadcasted out to A's length using the key variable
% indices.
c = a;
c.nvars = length(leftvars) + length(rightvars);
c.varnames = [a.varnames(leftvars) b.varnames(rightvars)];
c.data = [a.data(leftvars) cell(1,length(rightvars))];
for j = 1:length(rightvars)
    var_j = b.data{rightvars(j)};
    szOut = size(var_j); szOut(1) = a.nobs;
    c.data{length(leftvars)+j} = reshape(var_j(loc,:),szOut);
end
if ~isempty(a.props.Units) || ~isempty(b.props.Units)
    units = repmat({''},1,c.nvars);
    if ~isempty(a.props.Units)
        units(1:length(leftvars)) = a.props.Units(leftvars);
    end
    if ~isempty(b.props.Units)
        units((length(leftvars)+1):c.nvars) = b.props.Units(rightvars);
    end
    c.props.Units = units;
end
