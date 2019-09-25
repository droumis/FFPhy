% PARTCORM: Given a full [r x r] correlation matrix, produces a corresponding 
%           [r-1 x r-1] partial correlation matrix from which the covariate has 
%           been removed.
%
%     Usage: pc = partcorm(corrmat,{v})
%
%           corrmat = full correlation matrix.
%           v =       address (subscript) of the covariate [default = r].
%           -------------------------------------------------------------
%           pc =      reduced partial-correlation matrix.
%

function pc = partcorm(corrmat,v)
  if (nargin < 2) v = []; end;

  [r,c] = size(corrmat);

  if (~iscorr(corrmat))
    error('  PARTCORM: input matrix must be a correlation matrix.');
  end;

  if (isempty(v))
    v = r;
  end;

  if (v~=r)                             % Move var to be partialled, to last position
    rx = corrmat(v,:);
    cx = corrmat(:,v);
    corrmat(v,:) = corrmat(r,:);
    corrmat(:,v) = corrmat(:,r);
    corrmat(r,:) = rx;
    corrmat(:,r) = cx;
  end;
  r = r-1;

  [r12,i,j] = trilow(corrmat(1:r,1:r));
  r13 = zeros(size(r12));
  r23 = zeros(size(r12));

  for k = 1:length(i)
    r13(k) = corrmat(i(k),r+1);
    r23(k) = corrmat(j(k),r+1);
  end;

  pc = partcorr(r12,r13,r23);
  pc = trisqmat(pc,1);

  return;
