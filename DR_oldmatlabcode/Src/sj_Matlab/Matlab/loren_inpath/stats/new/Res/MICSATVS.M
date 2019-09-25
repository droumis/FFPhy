% MICSATVS: Calculates variance statistics (W,I,S,T) from allele frequencies of a
%           single locus for microsatellite data, based on a nested 
%           random-effects ANOVA of allele sizes.  Ignores zeros (=missing data).
%
%     Usage: [W,I,S,T] = micsatvs(alleles,grps)
%
%           alleles = [n x 2] matrix of allele frequencies for a single locus
%           grps =    [n x 1] vector of group-membership identifiers
%           --------------------------------------------------------------------
%           W =       within-individual, among-allele variance component
%           I =       within-population, among-individual variance component
%           S =       among-population variance component
%           T =       total variance in allele sizes
%

% Justification: Slatkin, M. 1995. A measure of population subdivision based
%   on microsatellite allele frequencies.  Genetics 139:457-462.
% Computations:  Sokal & Rohlf (1981: 294-299): 2-level random-effects nested 
%   anova with unequal sample sizes.

% RE Strauss, 3/9/97

function  [W,I,S,T] = micsatvs(alleles,grps)
  if (size(alleles,2)~=2)
    error('  MICSATVS: single locus only.');
  end;

  indx = (alleles(:,1)==0 | alleles(:,2)==0); % Delete obs with missing data
  if (sum(indx) > 0)
    alleles = alleles(~indx,:);
    grps = grps(~indx);
  end;

  grpid = uniquef(grps);              % Group identifiers
  ngrps = length(grpid);              % Number of groups
  nobs = size(alleles,1);             % Number of observations and loci

  x = zeros(2*nobs,1);
  g = zeros(2*nobs,1);
  s = zeros(2*nobs,1);

  for i = 1:nobs
    j = i*2-1;
    x(j:(j+1)) = alleles(i,:)';
    g(j:(j+1)) = [grps(i) grps(i)]';
    s(j:(j+1)) = [i i]';
  end;

  [F,pr,df,ss,ms,varcomp] = anovanst(x,g,s);   % Nested random-effects anova
  W = varcomp(3);
  I = varcomp(2);
  S = varcomp(1);
  T = sum(varcomp);

  return;
