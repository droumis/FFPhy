function stats = bf_rank_test(obs,labels,varargin)
%BF_RANK_TEST Two-sample rank test in a Behrens-Fisher situation.
%
%   BF_RANK_TEST performs a nonparametric test on the hypothesis that the two
%   samples are drawn from distributions with the same center but possibly
%   unequal variances. Note that this test is different from the Wilcoxon rank
%   sum test, which tests the null hypothesis that the samples are drawn from
%   the *same* distribution; the Wilcoxon test, unlike the present test, is
%   sensitive to heteroscedasticity. For small sample sizes, a permutation test
%   is used to compute the p-value; otherwise, the null distribution of the test
%   statistic is approximated with a t-distribution.
%
%   STATS = BF_RANK_TEST(X1,X2) returns a STATS struct with the following fields
%   (which roughly match the variable names used by Neubert & Brunner)
%        n: vector of group sizes [n1, n2]
%       TN: test statistic; positive-valued when X2 > X1 
%        f: degrees of freedom for t distribution which approximates the null
%           distribution of TN (Satterthwaite-Smith-Welch approximation). If 
%           METHOD is not specified, this approximation is used by default to 
%           compute the p-value when n1 >= 10 and n2 >= 10.
%     nsim: number of random permutations drawn for computing p-value. If
%           METHOD is not specified, the permutation test is used by default 
%           when n1 < 10 or n2 < 10. If NSIM is not provided, the default 
%           number of permutations is 1e6.
%        p: p-value for a 2-sided comparison
%   method: 'approximation' or 'permutation'
%
%   BF_RANK_TEST(X1,X2,'approximation') uses the Satterthwaite-Smith-Welch
%   approximation to compute the p-value. This argument overrides the default
%   behavior for small sample sizes.
%
%   BF_RANK_TEST(X1,X2,'permutation') uses a permutation test to compute the
%   p-value. The test statistic TN that is computed from the original
%   observations is compared to the permutational distribution of TN*
%   recalculated for all permutations of the labels. By default, 1e6 random
%   permutations are computed.
%
%   BF_RANK_TEST(X1,X2,'permutation',NSIM) does the permutation test with NSIM
%   random permutations.
%
%   The test is described in Neubert K., Brunner E. (2007) A studentized
%   permutation test for the non-parametric Behrens-Fisher problem.
%   _Computational Statistics & Data Analysis_ 51: 5192-5204.
%
%Depends on:
%   TIEDRANK (MATLAB Statistics toolbox)
%   TCDF (MATLAB Statistics toolbox)
%
%Written by smk, 4 January 2009
%

if (exist('tiedrank') ~= 2)
  error(['BF_RANK_TEST depends on m-file TIEDRANK ' ...
      '(in MATLAB Statistics Toolbox)']);
end
if (exist('tcdf') ~= 2)
  error(['BF_RANK_TEST depends on m-file TCDF ' ...
      '(in MATLAB Statistics Toolbox)']);
end

% default number of permutation draws if nsim is not specified
DEFAULT_NSIM = 1e6;

if ~isnumeric(x1) || ~isnumeric(x2) || ~isvector(x1) || ~isvector(x2)
  error('data arguments must be vectors')
end
% number of subjects
n1 = numel(x1); 
n2 = numel(x2); 
stats.n = [n1 n2];
% pick default method according to sample size
if (n1 <= 10) || (n2 <= 10)
  stats.method = 'permutation';
  stats.nsim = DEFAULT_NSIM; 
else
  stats.method = 'approximation';
end

if (length(varargin) == 1)
  if ischar(varargin{1}) && ...
      any(strcmp(varargin{1},{'permutation','approximation'}))
    stats.method = varargin{1};
    stats.nsim = DEFAULT_NSIM;
  else
    error('invalid METHOD argument specified');
  end
elseif (length(varargin) == 2)
  if ~strcmp(stats.method,'permutation')
    error('NSIM is a valid argument only with the permutation method');
  end
  if isnumeric(varargin{2}) && isscalar(varargin{2}) && ...
      isreal(varargin{2}) && (varargin{2} > 0) && ...
      (round(varargin{2}) == varargin{2})
    stats.nsim = varargin{2};
  else
    error('NSIM must be a positive integer');
  end
else
  error('too many arguments');
end

% concatenate data into a long vector
x = [x1(:); x2(:)]; 
N = n1+n2;
% vector of 1s and 2s as integer labels
groupidx = [repmat(1,size(x1(:))); repmat(2,size(x2(:)))];
if strcmp(stats.method,'approximation')
  [stats.TN, var1, var2] = studentizedstat(x,groupidx);
  % under the null hypothesis, TN obeys approximately a 
  % central t distribution with degrees of freedom given by
  stats.f = ( var1/(N-n1) + var2/(N-n2) )^2 / ...
      ( (var1/(N-n1))^2 / (n1-1) + (var2/(N-n2))^2 / (n2-1) );
  % compute a two-sided p-value, exploiting the convenient 
  % property that the t distribution is even
  stats.p = 2*(1 - tcdf(abs(stats.TN),stats.f));
elseif strcmp(stats.method,'permutation')
  stats.TN = studentizedstat(x,groupidx);
  for i = 1:stats.nsim
    % permute the observations
    Tperm(i) = studentizedstat(x(randperm(N)),groupidx);
  end
  % two-sided test
  stats.p = nnz((Tperm <= -abs(stats.TN)) | (Tperm >= abs(stats.TN))) / ...
      stats.nsim;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction for computing studentized test statistic
function [TN, var1, var2] = studentizedstat(x,groupidx)
n1 = nnz(groupidx==1);
n2 = nnz(groupidx==2);
N = n1+n2;
% first, rank *all* observations. note that the MATLAB
% tiedrank function handles ties correctly
Rik = tiedrank(x);
% mean of ranks for the two samples
R1_ = mean(Rik(groupidx==1));
R2_ = mean(Rik(groupidx==2));
% ranks among observations *within* each sample
R1k = tiedrank(x(groupidx==1));
R2k = tiedrank(x(groupidx==2));
% estimated variances
var1 = 1/(n1-1) * sum(( ...
    Rik(groupidx==1) - R1k - repmat(R1_,[n1 1]) + (n1+1)/2 ).^2);
var2 = 1/(n2-1) * sum(( ...
    Rik(groupidx==2) - R2k - repmat(R2_,[n2 1]) + (n2+1)/2 ).^2);
varN = N/(n1*n2) * (n1*var1 + n2*var2);
% test statistic
TN = (R2_ - R1_)/sqrt(varN) * sqrt(n1*n2/N); 

