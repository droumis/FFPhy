function r = mnrnd(n,p,m)
%MNRND Random vectors from the multinomial distribution.
%   R = MNRND(N,PROB) returns random vectors chosen from the multinomial
%   distribution with sample sizes N and probabilities PROB. PROB is an M-by-K
%   matrix or a 1-by-K vector of multinomial probabilities, where K is the
%   number of multinomial bins or categories.  Each row of PROB must sum to
%   one.  N is an M-by-1 vector of positive integers or a positive scalar
%   integer.  R is an M-by-K matrix, and MNRND generates each row of Y using
%   the corresponding rows of the inputs, or replicates them if needed.
%
%   R = MNRND(N,PROB,M) returns an M-by-K matrix of random vectors.
%
%   See also MNPDF, MNRFIT, MNRVAL, RANDSAMPLE

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2006/06/20 20:51:15 $

if nargin < 2
    error('stats:mnrnd:TooFewInputs', ...
          'Requires two input arguments.');
elseif nargin < 3
    [m,k] = size(p);
elseif ~isscalar(m)
    error('stats:mnrnd:NonscalarM', ...
          'M must be a scalar.');
else
    [mm,k] = size(p);
    if ~(mm == 1 || mm == m)
        error('stats:mnrnd:InputSizeMismatch', ...
              'P must be a row vector or have M rows.');
    end
end
if k < 1
    error('stats:mnrnd:NoCategories', ...
          'P must have at least one column.');
end

[mm,kk] = size(n);
if kk ~= 1
    error('stats:mnrnd:InputSizeMismatch', ...
          'N must be a scalar, or a column vector with as many rows as P.');
elseif m == 1 && ~isscalar(n)
    m = mm; % p will replicate out to match n
end

outClass = superiorfloat(n,p);

edges = [zeros(size(p,1),1) cumsum(p,2)];
pOK = all(0 <= p & p <= 1, 2) & (abs(edges(:,end)-1) <= eps);
nOK = (0 <= n & round(n) == n);

% If all cases have the same size and probs, histc can do them all at once
if isscalar(n) && isvector(p)
    if pOK && nOK
        r = histc(rand(m,n),edges,2);
        if strcmp(outClass,'single'), r = single(r); end
    else
        r = NaN(m,k,outClass);
    end

% Otherwise, treat each case individually
else
    r = NaN(m,k+1,outClass);
    if isvector(p) % && ~isscalar(n)
        if pOK
            for i = 1:m
                if nOK(i)
                    r(i,:) = histc(rand(1,n(i)),edges);
                end
            end
        end
    elseif isscalar(n) % && ~isvector(p)
        if nOK
            for i = 1:m
                if pOK(i)
                    r(i,:) = histc(rand(1,n),edges(i,:));
                end
            end
        end
    else
        for i = 1:m
            if pOK(i) && nOK(i)
                r(i,:) = histc(rand(1,n(i)),edges(i,:));
            end
        end
    end
end
r(:,end) = []; % remove histc's degenerate uppermost bin
