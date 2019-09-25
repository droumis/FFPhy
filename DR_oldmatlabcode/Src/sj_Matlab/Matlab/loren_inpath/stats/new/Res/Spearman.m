% SPEARMAN: Spearman's rank correlation between two variables, corrected for 
%           tied ranks.
%
%     Usage: rs = spearman(x,y)
%
%           x,y = data vectors.
%           --------------------------------
%           rs =  Spearman rank correlation.
%

% Zar (1984), p 320.

% RE Strauss, 12/26/99

function rs = spearman(x,y)
  x = x(:);
  y = y(:);

  n = length(x);
  if (length(y) ~= n)
    error('  SPEARMAN: data vectors not compatible');
  end;

  x = ranks(x);
  y = ranks(y);

  [ux,tx] = uniquef(x);
  [uy,ty] = uniquef(y);

  sumTx = sum(tx.^3 - tx)./12;
  sumTy = sum(ty.^3 - ty)./12;

  sumd2 = sum((x-y).^2);
  n3n = (n.^3-n)./6;

  num =   n3n-sumd2-sumTx-sumTy;
  denom = sqrt((n3n-2*sumTx).*(n3n-2*sumTy));
  if (denom > 0)
    rs = num ./ denom;
  else
    rs = 0;
  end;

  return;

