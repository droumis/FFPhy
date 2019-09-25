% WILCOXON: Wilcoxon's signed-ranks test for a single sample and a null median.
%
%     Usage: [pr,Tpos,Tneg] = wilcoxon(x,m0,tail)
%
%           x =          data vector.
%           m0 =         null median value.
%           tail =       direction of test (-1 = left-tailed, 0 = 2-tailed,
%                          +1 = right-tailed [default = 0].
%           ---------------------------------------------------------------
%           pr =         significance level.
%           N =          sample size for test.
%           Tpos, Tneg = sums of ranks of positive and negative deviations.
%

% RE Strauss, 5/8/99

function [pr,N,Tpos,Tneg] = wilcoxon(x,m0,tail)
  if (nargin < 3) tail = []; end;

  if (isempty(tail))
    tail = 0;
  end;

  if (min(size(x)) > 1)
    error('WILCOXON: data matrix must be a vector');
  end;
  if (max(size(m0)) > 1)
    error('WILCOXON: null median value must be a scalar');
  end;

  n = [4:30]';                                      % Critical values table
  alpha = [0.005, 0.01, 0.025, 0.05, 0.10, 0.50];   % Gardiner 1997, B7
  critval = [  0   0   0   0   1   5.0
               0   0   0   1   3   7.5
               0   0   1   3   4  10.5
               0   1   3   4   6  14.0
               1   2   4   6   9  18.0
               2   4   6   9  11  22.5
               4   6   9  11  15  27.5
               6   8  11  14  18  33.0
               8  10  14  18  22  39.0
              10  13  18  22  27  45.5
              13  16  22  26  32  52.5
              16  20  26  31  37  60.0
              20  24  30  36  43  68.0
              24  28  35  42  49  76.5
              28  33  41  48  56  85.5
              33  38  47  54  63  98.0
              38  44  53  61  70 105.0
              43  50  59  68  78 115.5
              49  56  66  76  87 126.5
              55  63  74  84  95 138.0
              62  70  82  92 105 150.0
              69  77  90 101 114 162.5
              76  85  99 111 125 175.5
              84  93 108 120 135 189.0
              92 102 117 131 146 203.0
             101 111 127 141 158 217.5
             110 121 138 152 170 232.5 ];

  Tpos = 0;
  Tneg = 0;

  d = x - m0;
  i = find(d == 0)
  if (~isempty(i))
    d(i) = [];
  end;

  N = length(d);
  if (N < min(n))
    error('WILCOXON: insufficient sample size');
  end;
  if (N > max(n))
    error('WILCOXON: need, but do not have, large-sample approximation');
  end;

  r = ranks(abs(d));
  
  i = find(d<0);
  if (~isempty(i))
    Tneg = sum(r(i));
  end;

  i = find(d>0);
  if (~isempty(i))
    Tpos = sum(r(i));
  end;

  T = min([Tpos Tneg]);

  id = find(n==N);
  critval = critval(id,:);

  if (T < critval(1))
    pr = alpha(1);
  elseif (T > critval(length(critval)))
    pr = alpha(length(critval));
  else
    pr = interp1(critval,alpha,T,'spline');
  end;

  if (tail==0)
    pr = 2*pr;
  end;

  return;

%$$ z_{+} = {{T_{+} - N(N+1)/4 - .5} \over {\sqrt{N(N+1)(2N+1)/24}}} $$

%$$ z_{-} = {{T_{-} - N(N+1)/4 - .5} \over {\sqrt{N(N+1)(2N+1)/24}}} $$

