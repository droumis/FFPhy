% STEPLENGTH: Given the p-dimensional coordinates of a path, finds the distances 
%             between consecutive points.
%
%     Usage: steplen = steplength(crds)
%
%         crds = [n x p] coordinates of a path through a p-dimensional space.
%         -------------------------------------------------------------------
%         steplen = column vector (length n-1) of step lengths.
%

% RE Strauss, 5/26/00

function steplen = steplength(crds)
  [n,p] = size(crds);

  d = [crds(2:n,:) - crds(1:(n-1),:)]';
  steplen = sqrt(sum(d.^2))';

  return;
