function m = trimmean(x,percent,dim)
%TRIMMEAN Trimmed mean.
%   M = TRIMMEAN(X,PERCENT) calculates the trimmed mean of the values in X.
%   For a vector input, M is the mean of X, excluding the highest and
%   lowest (PERCENT/2)% of the data.  For a matrix input, M is a row vector
%   containing the trimmed mean of each column of X.  For N-D arrays,
%   TRIMMEAN operates along the first non-singleton dimension.  PERCENT is
%   a scalar between 0 and 100.
%
%   TRIMMEAN(X,PERCENT,DIM) takes the trimmed mean along dimension DIM of X.
%
%   The trimmed mean is a robust estimate of the sample location.
%
%   TRIMMEAN treats NaNs as missing values, and removes them.
%
%   See also MEAN, NANMEAN, IQR.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 2.10.2.3 $  $Date: 2004/07/28 04:39:40 $

if nargin < 2
    error('stats:trimmean:TooFewInputs', 'TRIMMEAN requires two input arguments.');
elseif percent >= 100 || percent < 0
    error('stats:trimmean:InvalidPercent', 'PERCENT must be between 0 and 100.');
end

if nargin < 3 || isempty(dim)
    % The output size for [] is a special case, handle it here.
    if isequal(x,[]), m = NaN; return; end;

    % Figure out which dimension nanmean will work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Keep track of columns that were all missing data, or length zero.
allmissing = all(isnan(x),dim);

% Need to tile the output of prctile to trim X.
tile = ones(1,max(ndims(x),dim));
tile(dim) = size(x,dim);

% Find the upper and lower percentiles of X, and trim out values that are
% more extreme.
zlo = repmat(prctile(x, (percent / 2),dim), tile);
zhi = repmat(prctile(x, 100 - percent / 2,dim), tile);
x(x < zlo | zhi < x) = NaN;

% Compute the mean of X, excluding the trimmed values.
m = nanmean(x,dim);

% Warn if everything was trimmed, but not if all missing to begin with.
alltrimmed = (all(isnan(x),dim) & ~allmissing);
if any(alltrimmed(:))
    if all(alltrimmed(:))
        warning('stats:trimmean:NoDataRemaining','No data remain after trimming.');
    else
        warning('stats:trimmean:NoDataRemaining','No data remain in some columns after trimming.');
    end
end
