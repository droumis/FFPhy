% BASEDIFF: Given a matrix of nucleotide sequence data, or a transition matrix 
%           of numbers of base differences between two taxa, or a vector of base 
%           differences, returns the proportional differences between sequences 
%           and its standard error.
%
%     Usage: [diff,se] = basediff(T)
%
%         T =     [n x m] matrix (m~=10) of base sequences of length m, for n 
%                   taxa, with bases represented by four character values or 
%                   integers (designations immaterial) and missing bases by NaNs.
%                     OR
%                 [4 x 4] symmetric matrix of numbers of base differences between 
%                   two taxa, with numbers of base identities on the diagonal.
%                     OR
%                 [n x 10] matrix of counts, of which each row represents a 
%                   different pairwise contrast between two taxa (in no 
%                   particular sequence).  The first 6 columns are counts of 
%                   differences between taxa, and the last 4 are numbers of bases 
%                   that are identical between taxa.  
%                   Sequence of base designations is immaterial.
%         -----------------------------------------------------------------------
%         diff =  scalar or [n x 1] vector or [n x n] symmetric matrix of 
%                   proportional difference (sum of all differences divided by 
%                   total number of bases).
%         se =    corresponding scalar or matrix of standard errors.
%

% RE Strauss, 10/20/00
%   1/23/01 - added input in form of sequence data.
%   5/29/02 - allow for character input; provide scalar output for 2 taxa.

function [diff,se] = basediff(T)
  retsym = 0;

  T = double(T);
  
  if (isvector(T))
    t = T(:)';
  else
    [r,c] = size(T);
    if (c==10)
      t = T;
    elseif (issqsym(T))
      t = [trilow(T); diag(T)]';
    else
      t = basediffp(T);
      retsym = 1;
    end;
  end;

  if (~isintegr(t))
    error('  BASEDIFF: input matrix must consist of counts (integers)');
  end;

  n = size(t,1);
  diff = zeros(n,1);
  se = zeros(n,1);

  for i = 1:n
    tt = t(i,:)';
    N = sum(tt);
    diff(i) = sum(tt(1:6))/N;
    se(i) = sqrt(diff(i)*(1-diff(i))/N);
  end;

  if (retsym)
    diff = trisqmat(diff);
    se = trisqmat(se);
  end;
  
  if (size(diff)==[2 2])
    diff = diff(1,2);
    se = se(1,2);
  end;

  return;
