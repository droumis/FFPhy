function stats = rm_rank_test(obs,labels)
% RM_RANK_TEST Repeated-measures rank test.
%
%   RM_RANK_TEST performs a nonparametric test for a repeated-measures,
%   between-groups design, in which N subjects are partitioned into A groups and
%   observed over T timepoints.
%
%   OBS is a NxT array. Note that every subject must be at all K timepoints;
%   missing observations are not permitted. 
%
%   LABELS is a cell array or character array whose N rows are string labels for
%   grouping subjects; i.e. numel(unique(labels)) == A. If the labels argument
%   is omitted, or if A==1, then the function tests for a simple time effect
%   under the assumption that all subjects belong to the same group.
%
%   STATS = RM_RANK_TEST(OBS,LABELS) returns a struct with the following fields
%   (which correspond roughly to the variable names used in Brunner et al.)
%   groups: cell array of group labels
%        n: vector of group sizes
%       FT: F test statistic for average time effect
%       fT: 1st  df for null distribution of FT
%       pT: p-value for average time effect
%       FA: F test statistic for main effect of group
%       fA: 1st df for null distribution of FA
%       f0: 2nd df for null distribution of FA
%       pA: p-value for main effect of group
%      FAT: F test statistic for groupxtime interaction
%      fAT: 1st df for null distribution of FAT
%      pAT: p-value for groupxtime interaction
%
%   The test is described in section 8.3 (pages 130-152) Brunner E., Domhof S.,
%   Langer F. (2002) Nonparametric Analysis of Longitudinal Data in Factorial
%   Experiments. Wiley-Interscience
%
%Depends on:
%   TIEDRANK (MATLAB Statistics toolbox)
%   FCDF (MATLAB Statistics toolbox)
%   GRP2IDX (MATLAB Statistics toolbox)
%
%Written by smk, 18 October 2008
%

if (exist('tiedrank') ~= 2)
  error(['BF_RANK_TEST depends on m-file TIEDRANK ' ...
      '(in MATLAB Statistics Toolbox)']);
end
if (exist('fcdf') ~= 2)
  error(['BF_RANK_TEST depends on m-file FCDF ' ...
      '(in MATLAB Statistics Toolbox)']);
end
if (exist('grp2idx') ~= 2)
  error(['RM_RANK_TEST depends on m-file GRP2IDX ' ...
      '(in MATLAB Statistics Toolbox)']);
end

% if labels are not provided, create a labels array whose elements are all the
% empty string ''
if (nargin == 1)
  labels = cellstr(char(ones([size(obs,1),1])));
end
if ~isnumeric(obs) 
  error('observations must be numeric');
end
if any(isnan(obs(:)))
  error('NaN values are not permitted');
end
if ischar(labels) && ~isempty(labels)
  % convert character array to cell array
  labels = cellstr(labels);
elseif iscellstr(labels)
  % pass
else
  error('labels argument must be a character array or cell array of string');
end

% total number of observations
N = numel(obs);
% number of subjects and number of timepoints
[n t] = size(obs);
if size(labels,1) ~= n
  error('obs and labels must have same number of rows (subjects)');
end
% count groups
[groupidx, stats.labels] = grp2idx(labels);
% number of groups
a = max(groupidx);
% count number of subjects in each group
ni = diff(find([1; diff([sort(groupidx); a+1])]));
if n ~= sum(ni)
  error('bug: subjects are not counted correctly');
end
stats.n = ni;

% centering matrices
Pt = eye(t) - ones(t)/t;
Pa = eye(a) - ones(a)/a;
% first, rank *all* observations. note that the MATLAB tiedrank function handles
% ties correctly
Riks = reshape(tiedrank(obs(:)),size(obs));
% within-group means of ranks over all subjects in group i at each timepoint
for i = 1:a
  Ri_s(i,:) = mean(Riks(groupidx==i,:),1);
end
% unweighted mean of the group rank means; note that each group contributes
% equally regardless of its number of subjects
R__s = mean(Ri_s,1);
% covariance between timepoints for each group
for i = 1:a
  temp = Riks(groupidx==i,:);
  res = Riks(groupidx==i,:) - repmat(Ri_s(i,:),[ni(i) 1]);
  % size(V) = [t t a]
  Vi(:,:,i) = n/(N^2*ni(i)*(ni(i)-1)) * res' * res;
end
% covariance matrix for the average time effect over all groups
St = 1/a^2 * sum(Vi,3);
% ANOVA-type test statistic for average time effect
stats.FT = n/(N^2 * trace(Pt*St)) * ...
  sum((R__s - mean(R__s)).^2);
% the null distribution of FT can be approximated by a central F(fT,Inf)
% distribution with the following (non-integer) 1st degree of freedom:
stats.fT = (trace(Pt*St))^2 / trace(Pt*St*Pt*St);
stats.pT = 1 - fcdf(stats.FT,stats.fT,Inf);

% test group effects only if there is more than one group
if a > 1
  % within-subject mean of ranks over all timepoints
  Rik_ = mean(Riks,2);
  % within-group means of ranks over all subjects and 
  % over all timepoints
  Ri__ = mean(Ri_s,2);
  % within-group variance of ranks over all timepoints and over all subjects in
  % group i
  for i = 1:a
    vari(i,1) = var(Rik_(groupidx==i),0);
  end
  % Vn is the kronecker sum of Vi(:,:,i)
  Vn = zeros(t*a);
  for i = 1:a
    blockidx = (t*i-t+1):(t*i);
    Vn(blockidx,blockidx) = Vi(:,:,i);
  end

  % ANOVA-type statistic for the main effect of group. Note that this formula
  % includes the *unweighted* mean of the group means Ri__; this is not a
  % mistake! Each group contributes equally regardless of its number of subjects
  stats.FA = a/((a-1)*sum(vari./ni)) * sum((Ri__ - mean(Ri__)).^2);
  % the null distribution of FA can be approximated by a central F(fA,f0)
  % distribution with the following (non-integer) degrees of freedom:
  stats.fA = (a-1)^2 / (1 + a*(a-2)*(sum((vari./ni).^2)/(sum(vari./ni))^2));
  stats.f0 = (sum(vari./ni))^2 / (sum((vari./ni).^2 ./ (ni-1)));
  stats.pA = 1 - fcdf(stats.FA,stats.fA,stats.f0);

  % ANOVA-type test statistic for group x time interaction
  TAT = kron(Pa,Pt);
  stats.FAT = n/(N^2 * trace(TAT*Vn)) * ...
      sum(sum((Ri_s - repmat(Ri__,[1 t]) - repmat(R__s,[a 1]) ...
      + repmat(mean(R__s),[a t])).^2,2),1);
  % the null distribution of FAT can be approximated by a central F(fAT,Inf)
  % distribution with the following (non-integer) 1st degree of freedom:
  stats.fAT = (trace(TAT*Vn))^2 / trace(TAT*Vn*TAT*Vn);
  stats.pAT = 1 - fcdf(stats.FAT,stats.fAT,Inf);
end


