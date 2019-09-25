% CORRCLST: Finds the optimal permutation for a square symmetric matrix, as
%           measured by Moran's spatial autocorrelation index I, either by
%           complete enumeration or by simulated annealing.
%
%       Usage: [permute,obj,rc,stats,distrib] = corrclst(r,{nbest},{maxperm})
%
%           r =       [n x n] correlation matrix
%           nbest =   saves the n-best permutations (default=1)
%           maxperm = max matrix order for exact solution by permutations
%                       (default = 7)
%           permute = [nbest x n] best row/col permutations
%           obj =     [nbest x 1] vector of objective-function values
%           rc =      permuted correlation matrix
%           stats =   vector of obj-fn values (for enumeration) or matrix of
%                       performance statistics (for simulated annealing)
%           distrib = vector of sampling distribution of objective-function 
%                       values (default = min(1000 or n!))
%
%       Calls function CORRPERM for enumeration and CORRANEL for simulated
%       annealing.
%

% RE Strauss, 2/25/95

function [permute,obj,rc,stats,distrib] = corrclst(r,nbest,maxperm)
  default_maxperm = 7;
  default_nbest = 1;

  if (nargin < 3)
    maxperm = default_maxperm;
  end;
  if (nargin < 2)
    nbest = default_nbest;
  end;

  n = size(r,1);                            % Size of correlation matrix
  if (n <= maxperm)                         % If order less than max,
    disp('  using enumeration...');         %   evaluate all permutations
    [permute,obj,rc,stats,distrib] = corrperm(r,nbest);
    disp(' ');
    return;
  else
    disp('  using simulated annealing...'); % Else use simulated annealing
    [permute,obj,rc,stats,distrib] = corranel(r,nbest);
    disp(' ');
    return;
  end;

  return;
