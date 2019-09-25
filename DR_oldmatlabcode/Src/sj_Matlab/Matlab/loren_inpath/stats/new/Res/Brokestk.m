% BROKESTK: Predicts the eigenvalues from a principal component analysis of 
%           random data based on the null broken-stick model of Frontier (1976).  
%           If a vector of eigenvalues is supplied, the function returns the 
%           number of 'significant' eigenvalues.  See Jackson (1993).
%
%     [pred_evals,nvect] = brokestk(nvars,{totvar},{evals})
%
%           nvars =   number of variables in analysis.
%           totvar =  optional total variance, if the PCA is based on a 
%                       covariance matrix.
%           evals =   optional vector of observed eigenvalues.
%           -------------------------------------------------------------------
%           pred_evals = col vector of predicted eigenvalues from null model.
%           nvect =   number of 'significant' eigenvalues, if evals is passed.
%

% Frontier, S.  1976.  Etude de la decroissance des valeurs propres dans une 
%   analyze en composantes principales: comparison avec le modele de baton 
%   brise.  J. Exp. Mar. Biol. Ecol. 25:67-75.
% Jackson, D.A.  1993.  Stopping rules in principal components analysis: a 
%   comparison of heuristial and statistical approaches.  Ecology 74:2204-2214.

% RE Strauss, 12/7/99

function [pred_evals,nvect] = brokestk(nvars,totvar,evals)
  if (nargin < 2) totvar = []; end;
  if (nargin < 3) evals = []; end;

  if (~isempty(evals))
    if (length(evals)~=nvars)
      error('  BROKESTK: requires full set of observed eigenvalues');
    end;
  end;

  pred_evals = zeros(nvars,1);
  nvect = [];

  for i = 1:nvars                         % Predicted eigenvalues from broken-stick model
    pred_evals(i) = sum(1./(i:nvars));
  end;

  if (~isempty(totvar))                   % Adjust by total variance
    pred_evals = totvar .* pred_evals ./ nvars;
  end;

  if (~isempty(evals))                    % Number of significant eigenvalues
    b = (evals > pred_evals);
    i = min(find(~b));
    nvect = length(1:(i-1));
  end;

  return;
