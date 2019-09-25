function [h,p,ci,stats] = ttest(x,m,alpha,tail)
%TTEST  Hypothesis test: Compares the sample average to a constant.
%   [H,P,CI,STATS] = TTEST(X,M,ALPHA,TAIL) performs a T-test to determine
%   if a sample from a normal distribution (in X) could have mean M.
%   M = 0, ALPHA = 0.05 and TAIL = 0 by default.
%
%   The Null hypothesis is: "mean is equal to M".
%   For TAIL=0,  alternative: "mean is not M".
%   For TAIL=1,  alternative: "mean is greater than M"
%   For TAIL=-1, alternative: "mean is less than M"
%   TAIL = 0 by default.
%
%   ALPHA is desired significance level. 
%   P is the p-value, or the probability of observing the given result
%     by chance given that the null hypothesis is true. Small values
%     of P cast doubt on the validity of the null hypothesis.
%   CI is a confidence interval for the true mean.  Its confidence
%     level is 1-ALPHA.
%   STATS is a structure with two elements named 'tstat' (the value
%     of the test statistic) and 'df' (its degrees of freedom).
%
%   H=0 => "Do not reject null hypothesis at significance level of alpha."
%   H=1 => "Reject null hypothesis at significance level of alpha."

%   References:
%      [1] E. Kreyszig, "Introductory Mathematical Statistics",
%      John Wiley, 1970, page 206. 

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 2.14 $  $Date: 2002/05/08 18:43:58 $

if nargin < 1, 
    error('Requires at least one input argument.'); 
end

[m1 n1] = size(x);
if (m1 ~= 1 & n1 ~= 1) 
    error('First argument has to be a vector.');
end
x = x(~isnan(x));

if nargin < 2
    m = 0;
end
 
if nargin < 4, 
    tail = 0; 
end 

if nargin < 3, 
    alpha = 0.05; 
end 

if (prod(size(alpha))>1), error('ALPHA must be a scalar.'); end
if (alpha<=0 | alpha>=1), error('ALPHA must be between 0 and 1'); end

samplesize  = length(x);
xmean = mean(x);
ser = std(x) ./ sqrt(samplesize);
tval = (xmean - m) / ser;
if (nargout > 3), stats = struct('tstat', tval, 'df', samplesize-1); end;
if isnan(tval)
   p = NaN;
else
   p = tcdf(tval,samplesize - 1);
end

% the p-value just found is for the  tail = -1 test

if (tail == 0)
    p = 2 * min(p, 1-p);
    if (nargout>2)
        crit = tinv((1 - alpha / 2),samplesize - 1) .* ser;
        ci = [(xmean - crit) (xmean + crit)];
    end
else
    if tail == 1
        p = 1 - p;
        if (nargout>2)
            crit = tinv(1 - alpha,samplesize - 1) .* ser;
            ci = [(xmean - crit), Inf];
        end
    else
        if (nargout>2)
            crit = tinv(1 - alpha,samplesize - 1) .* ser;
            ci = [-Inf, (xmean + crit)];
        end
    end
end


% Determine if the actual significance exceeds the desired significance
h = 0;
if p <= alpha, 
    h = 1; 
end 

if isnan(p), 
    h = NaN; 
end
