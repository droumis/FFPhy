% RANDINDX: Adjusted Rand index for the comparison of different data
%           partitions.  Uses Krauth's (1989) modifications to account
%           for chance agreements and to standardize the range.
%
%           Properties: range -1 <= A <= +1;
%                       A = 0 in the case of independence and simultaneous
%                           identity of marginal distributions;
%
%     Usage: A = randindx(part1,part2)
%
%           part1 = group-identification vector (length n) indicating
%                     the first partition.
%           part2 = matching vector indicating the second partition.
%           ------------------------------------------------------------
%           A =     adjusted Rand index value
%

% Krauth, J. 1989. A new modification of the Rand index for comparing
%   partitions.  In: O. Opitz (ed.), Conceptual and Numerical Analysis
%   of Data, pp. 135-146. Springer-Verlag.

% RE Strauss, 1/3/99

function A = randindx(part1,part2)
  [r1,c1] = size(part1);
  [r2,c2] = size(part2);

  err = 0;
  if (min([r1,c1])>1 | min([r2,c2])>1)
    err = 1;
  end;
  len1 = max([r1,c1]);
  len2 = max([r2,c2]);
  if (len1~=len2)
    err = 1;
  end;
  if (err)
    error('RANDINDX: input matrices must be matching vectors');
  end;

  ntot = len1;
  u1 = uniquef(part1);
  u2 = uniquef(part2);
  k1 = length(u1);
  k2 = length(u2);

  n = zeros(k1,k2);
  for i = 1:ntot
    p1 = find(u1==part1(i));
    p2 = find(u2==part2(i));
    n(p1,p2) = n(p1,p2)+1;
  end;

  sumrowcomb = 0;
  sumcolcomb = 0;
  sumcellcomb = 0;
  combtot = comb(ntot,2);

  for i = 1:k1
    sumrowcomb = sumrowcomb + comb(sum(n(i,:)),2);
    for j = 1:k2
      if (i==1)
        sumcolcomb = sumcolcomb + comb(sum(n(:,j)),2);
      end;
      sumcellcomb = sumcellcomb + comb(n(i,j),2);
    end;
  end;

  m00 = combtot - sumrowcomb - sumcolcomb + sumcellcomb;
  m01 = sumcolcomb - sumcellcomb;
  m10 = sumrowcomb - sumcellcomb;
  m11 = sumcellcomb;
  m0dot = combtot - sumrowcomb;
  m1dot = sumrowcomb;
  mdot0 = combtot - sumcolcomb;
  mdot1 = sumcolcomb;
  m = combtot;

  num = (1/m)*(m00+m11);                              % A
  num = num - (1/(m*m))*(m0dot*m0dot + m1dot*m1dot);  % E(A)
  num = num - (1/(m*m))*(m01-m10)*(m01-m10);          % Adj for diff clust sizes

  den = 1 - (1/(m*m))*(m0dot*m0dot+m1dot*m1dot);      % Max value

  A = num/den;

  return;
