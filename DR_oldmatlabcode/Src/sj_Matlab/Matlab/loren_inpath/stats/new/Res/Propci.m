% PROPCI: Confidence intervals for proportions and differences between 
%         proportions.
%
%     Usage: [lower,upper] = propci(numer,denom,{paired},{alpha})
%
%         numer =         vector (length k) of numerators of proportions.
%         denom =         corresponding vector of denominators.
%         paired =        optional boolean vector, for k=2, indicating that the
%                           two proportions are paired (eg, pretest & posttest) 
%                           [default = 0].
%         alpha =         confidence level [default = 95].
%         ---------------------------------------------------------------------
%         lower, upper =  limits of confidence interval: scalars for a single 
%                           proportion or two paired proportions, or [k x k]
%                           matrices of pairwise limits among k proportions.
%

% CI for single proportion:  EB Wilson. 1927. JASA 22:209-212.
% CI for difference between paired or independent proportions:  RG Newcombe. 
%     1998. Stats in Medicine 17:2635-2650.  Method 10.

% RE Strauss, 5/17/99

function [lower,upper] = propci(numer,denom,paired,alpha)
  if (nargin < 3) paired = []; end;
  if (nargin < 4) alpha = []; end;

  if (isempty(paired))
    paired = 0;
  end;
  if (isempty(alpha))
    alpha = 95;
  end;

  if (min(size(numer))>1 | min(size(denom))>1)
    error('PROPCI: numerators and denominators must be vectors');
  end;

  k = length(numer);
  if (length(denom)~=k)
    error('PROPCI: numerator and denominator vectors not compatible');
  end;


  return;
