% RandState: Given the transitional probabilities of changing from the current
%   state to the next state, predicts the next state.
%
%     Usage:  nextstate = randstate(probabilities)
%
%         probabilities = [1 x n] matrix of transitional probabilities
%                           for n states.
%         ------------------------------------------------------------
%         nextstate = scalar indicating the next state, from 1 to n.
%

function nextstate = randstate(probabilities)
  if (sum(probabilities)~=1)
    error('  RandState: probabilities must sum to 1.');
  end;
  
  p = cumsum(probabilities);
  r = rand;

  nextstate = min(find(p-r >= 0));

  return;
  