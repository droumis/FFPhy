% CONTINEX: Fisher's exact test for a 2x2 contingency table
%
%     Usage: p = continex(obs)
%
%           obs = [2 x 2] matrix of counts.
%           --------------------------------------------------------------------
%           p =   probability of being as or more extreme in cell-count
%                   distributions; returns NaN if the cell counts are too large.
%

% Ostle, B & RW Mensing. 1975. Statistics in Research, 3rd ed. Iowa State
%   University Press.  p 135.

% RE Strauss, 9/29/95

function p = continex(obs)
  p = obs(:,1)./sum(obs')';
  if (p(1)>p(2))
    obs = obs([2,1],:);
  end;

  a = obs(1,1);
  b = obs(1,2);
  c = obs(2,1);
  d = obs(2,2);

  n = sum(sum(obs));
  prob = zeros((a+1),1);

  num = sum(log(1:(a+b))) + sum(log(1:(c+d))) + sum(log(1:(a+c))) + ...
        sum(log(1:(b+d))) - sum(log(1:n));

  for i = 1:(a+1)
    denom = sum(log(1:(a-i+1))) + sum(log(1:(b+i-1))) + sum(log(1:(c+i-1))) + ...
            sum(log(1:(d-i+1)));
    prob(i) = exp(num-denom);
  end;

  p = sum(prob);

  return;
