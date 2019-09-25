% MANNW: Mann-Whitney two-sample rank-sum test
%
%     [T,m,n] = mannw(x,grps)
%

% RE Strauss, 4/26/99
%   11/29/99 - changed calling sequence.

function [T,m,n] = mannw(x,grps)
  r = ranks(x);
  g = uniquef(grps);
  
  i = find(grps==g(1));
  m = length(i);
  n = length(x)-m;
  W = sum(r(i));
  T = W - m*(m+1)/2;


  return;
