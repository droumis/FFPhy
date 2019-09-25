% COMBVALS: Generates the combinations of n integers taken r at a time.  The 
%           number of such combinations is given by function nc=combin().  
%
%     Usage: c = combvals(n,r)
%
%           n = number of integers (1:n) to be combined.
%           r = number to be taken at a time (0 < r <= n).
%           -------------------------------------------------------
%           c = [nc x r] matrix of combinations.
%

% Based on ACM Algorithm 94, J. Kurtzberg, Comm. ACM, June 1962.

% RE Strauss, 12/18/98

function c = combvals(n,r)
  n = round(n);
  r = round(r);

  if (n<1 | r<1 | r>n)
    error('COMBVALS: n or r out of range');
  end;

  ncomb = comb(n,r);
  c = zeros(ncomb,r);
  j = zeros(1,r);

  for i = 1:ncomb
    b = 1;
    endflag = 0;
    while(~endflag)
      if (j(b)>=b)
        a = j(b)-b-1;
        for l = 1:b
          j(l) = l+a;
        end;
        endflag = 1;;     
      else
        if (b==r)
          for b = 1:r
            j(b) = n-r-1+b;
          end;
          endflag = 1;
        end;
        b = b+1;
      end;
    end;
    c(i,:) = n-j(r:-1:1);
  end;

  return;
