function D = dummyvar(group)
%DUMMYVAR Dummy variable coding.
%   X=DUMMYVAR(GROUP) makes a matrix X containing a dummy column for each
%   unique value of the grouping variables in GROUP.  GROUP can be a
%   categorical variable, a cell array of multiple categorical variables,
%   or a matrix of grouping variable values.  If GROUP is a matrix, the
%   values of the elements in any column of GROUP go from one to the number
%   of members in the group defined by that column.
%
%   Example: Suppose we are studying the effects of two machines and three
%   operators on a process.  The first column of GROUP would have the
%   values one or two depending on which machine was used.  The second
%   column of GROUP would have the values one, two, or three depending on
%   which operator ran the machine.
%
%       machine = [1 1 1 1 2 2 2 2]';
%       oper    = [1 2 3 1 2 3 1 2]';
%       x = dummyvar([machine oper])

%   Copyright 1993-2006 The MathWorks, Inc. 
%   $Revision: 2.8.2.5 $  $Date: 2006/11/11 22:55:07 $

if isa(group,'categorical')
    % group by a categorical (nominal or ordinal) variable
    [m,n] = size(group);
    if n~=1
        error('stats:dummyvar:BadGroup',...
              'Categorical grouping variable must have one column.');
    end
    maxg = length(getlabels(group));
    group = double(group);
elseif iscell(group)
    % collection of grouping variables in a cell array
    n = numel(group);
    for j=1:n
        gj = group{j};
        gj = grp2idx(gj);
        if j==1
            m = size(gj,1);
            G = zeros(m,n);
        else
            if size(gj,1)~=m
                error('stats:dummyvar:InputSizeMismatch',...
                      'All grouping variables must the same number of rows.');
            end
        end
        G(:,j) = double(gj);
    end
    group = G;
    maxg = max(group);
else
    % vector or matrix of grouping variables
    [m,n] = size(group);
    if m == 1
        m = n;
        n = 1;
        group = group(:);
    end
    maxg = max(group);
end

if any(any(group - round(group))) ~= 0
   error('stats:dummyvar:BadGroup',...
         'Each element of GROUP must be a positive integer.')
end
if ~isfloat(group)
   group = double(group);
end

colstart = [0 cumsum(maxg)];
colstart(n+1) = [];
colstart = reshape(colstart(ones(m,1),:),m*n,1);

colD = sum(maxg);
D = zeros(m,colD,class(group));

row = (1:m)';
row = reshape(row(:,ones(n,1)),m*n,1);

idx = m*(colstart + group(:) - 1) + row;
D(idx) = 1;
