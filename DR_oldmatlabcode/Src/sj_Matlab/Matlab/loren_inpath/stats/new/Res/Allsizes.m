% Allsizes: Returns list of all possible decompositions of a total sample size 
%           of N into k group-sample sizes, omitting permutations.
%
%     Usage: sizes = allsizes(N,k)
%
%           N = total sample size.
%           k = number of groups.
%
%           sizes = [nsizes x k] matrix of all possible decompositions.
%

% RE Strauss, 8/29/99

function sizes = allsizes(N,k)
  s = [ones(1,k-1) N-k+1];
  sizes = s;

  smin = min(s);
  smax = max(s);

  i = 1;
  while (~isempty(i) & (smin < smax-1))
    i = max(find(s(1:k-1) < s(2:k)-1));

    if (~isempty(i))
      s(i) = s(i)+1;
      s(i+1) = s(i+1)-1;
      sizes = [sizes; s];
    else
      i = max(find(s < smax-1));
      j = min(find(s == smax));
      if (~isempty(i))
        s(i) = s(i)+1;
        s(j) = s(j)-1;
        sizes = [sizes; s];
      end;
    end;

    smin = min(s);
    smax = max(s);
  end;

  return;