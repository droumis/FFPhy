% SCORE: Calculates unstandardized factor scores for a subset of factors.
%
%     Usage: scr = score(X,evect,{nf})
%
%         X =     [n x p] data matrix.
%         evect = [p x p] matrix of eigenvectors.
%         nf =    optional number of leading factors for which scores needed
%                   [default = p].
%         -------------------------------------------------------------------
%         scr =   [p x nf] matrix of factor scores (columns).
%

% RE Strauss, 10/14/93
%   10/4/99 - minor changes to increase efficiency.

function scr = score(X,evect,nf)
  if (nargin < 3) nf = []; end;

  [n,p] = size(X);
  if (isempty(nf))                % If nf not passed,
    nf = p;                       %   set = p
  end;

  nf = max([nf,1]);               % Check range of nf
  nf = min([nf,p]);

  scr = X * evect(:,1:nf);

  return;
