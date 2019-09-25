function [p, h, stats] = ranksum(x,y,alpha)
%RANKSUM Wilcoxon rank sum test that two populations are identical.
%   P = RANKSUM(X,Y,ALPHA) returns the significance for testing the
%   null hypothesis that the populations generating two independent
%   samples, X and Y, are identical. X and Y are vectors but can have
%   different lengths.  The alternative is that the median of the X
%   population is shifted from the median of the Y population by a
%   non-zero amount.
%
%   ALPHA is the desired level of significance and must be a scalar
%   between zero and one.  Its default value is 0.05.
%
%   P is the probability of observing a result equally or more 
%   extreme than the one using the data (X and Y) if the null  
%   hypothesis is true. If P is near zero, this casts doubt on
%   this hypothesis.
%
%   [P, H] = RANKSUM(X,Y,ALPHA) also returns H, the result of the
%   hypothesis test.  H is 0 if the medians of X and Y are not
%   significantly different, and 1 if they are significantly
%   different.
%
%   [P, H, STATS] = RANKSUM(X,Y,ALPHA) also returns a STATS structure
%   with one or two fields.  The field 'ranksum' contains the value of
%   the rank sum statistic.  If the sample size is large, then P is
%   calculated using a normal approximation and the field 'zval'
%   contains the value of the normal (Z) statistic.

%   B.A. Jones 12-28-96
%   Copyright 1993-2002 The MathWorks, Inc. 
% $Revision: 1.10 $

if nargin < 3
   alpha = 0.05;
end

if (length(alpha)>1)
   error('RANKSUM requires a scalar ALPHA value.');
end
if ((alpha <= 0) | (alpha >= 1))
   error('RANKSUM requires 0 < ALPHA < 1.');
end

[nx, colx] = size(x);
[ny, coly] = size(y);

if min(nx, colx) ~= 1 | min(ny,coly) ~= 1,
   error('RANKSUM requires vector rather than matrix data.');
end 
if nx == 1
   nx = colx;
   x = x';
end
if ny == 1,
   ny = coly;
   y = y';
end

if nx <= ny
   smsample = x;
   lgsample = y;
   ns = nx;
   nl = ny;
else
   smsample = y;
   lgsample = x;
   ns = ny;
   nl = nx;
end

% Compute the rank sum statistic based on the smaller sample
[ranks, tieadj] = tiedrank([smsample; lgsample]);
xrank = ranks(1:ns);
w = sum(xrank);

wmean = ns*(nx + ny + 1)/2;
if ns < 10 & (nx+ny) < 20     % Use the sampling distribution of W.
   allpos = nchoosek(ranks,ns);
   sumranks = sum(allpos,2);
   np = length(sumranks);
   if w < wmean
      p = 2*length(find(sumranks <= w))./np;
   else 
      p = 2*length(find(sumranks >= w))./np;
   end
   p = min(p, 1);        % p>1 means w is in the middle and double-counted
else    % Use the normal distribution approximation of W.
   tiescor = 2 * tieadj / ((nx+ny) * (nx+ny-1));
   wvar  = nx*ny*((nx + ny + 1) - tiescor)/12;
   wc = w - wmean;
	z = (wc - 0.5 * sign(wc))/sqrt(wvar);
	p = normcdf(z,0,1);
	p = 2*min(p,1-p);
   if (nargout > 2)
      stats.zval = z;
   end
end

if nargout > 1,
   h = (p<=alpha);
   if (nargout > 2)
      stats.ranksum = w;
   end
end
