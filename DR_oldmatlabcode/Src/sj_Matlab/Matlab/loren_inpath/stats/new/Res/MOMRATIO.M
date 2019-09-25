% MOMRATIO: For given vector X, calculates the RANDSTRB sse over the
%           Pearson moment-ratio landscape
%
%       Usage: [b1,b2,sse] = momratio(X)
%

% RE Strauss, 5/5/95

function [b1,b2,sse] = momratio(X)

  g1_values = [-2:0.2:2];
  g2_values = [1:0.5:8];

  sse = zeros(length(g2_values)*length(g1_values),1);
  b1 =  zeros(length(g2_values)*length(g1_values),1);
  b2 =  zeros(length(g2_values)*length(g1_values),1);
  k = 0;

%  sse = zeros(length(b2),length(b1));

  i = 0;
  for g1 = g1_values
    i = i+1;
    j = 0;
    for g2 = g2_values
      j = j+1;
      [x,m,s,skew,kurt,l1,l2,e] = randstrb(X,25,5,g1,g2,1);
[g1 g2 skewness(x) kurtosis(x) e]
%      sse(j,i) = e;

      k = k+1;
      b1(k) = g1;
      b2(k) = g2;
      sse(k) = e;
    end;
  end;

%  figure(1), contour(b1,b2,sse);
%  axis([min(b1) max(b1) min(b2) max(b2)]);

%  figure(2), surf(b1,b2,sse);
%  axis([min(b1) max(b1) min(b2) max(b2) 0 max(max(sse))]);
%  view(15,30);

  return;

