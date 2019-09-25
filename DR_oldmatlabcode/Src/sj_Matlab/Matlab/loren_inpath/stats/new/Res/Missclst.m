% MISSCLST: Given a data matrix, number of modes, and proportion of missing 
%           data, randomly introduces missing values into matrix based on 
%           normal probabilities about the modes.
%
%     Usage: Xmiss = missclst(X,nval,[kvar,kobs,kboth],{intensity})
%
%         X =         [n x p] data matrix.
%         nval =      number of missing data values.
%         k =         3-element vector of number of cluster modes for variables, 
%                       observations, or both (only one of the three value 
%                       can be nonzero).
%         intensity = normalized (0-1) measure of the inverse variance of 
%                       the normal distributions [default = 0.5].
%         -----------------------------------------------------------------
%         Xmiss =     [n x p] data matrix with missing values.
%

function Xmiss = missclst(X,nval,k,intensity)
  if (nargin < 4) intensity = []; end;

  if (isempty(intensity))
    intensity = 0.5;
  end;



  return;
