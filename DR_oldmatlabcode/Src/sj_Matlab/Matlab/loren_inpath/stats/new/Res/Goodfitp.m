% GOODFITP: Performs a randomized goodness-of-fit test based on the mean
%           squared deviation between observed and expected values
%           (Edgington 1987; Manly 1991:170-171).  Equivalent to a 1-way 
%           contingency analysis.
%             See also GOODFIT, which uses a chi-square statistic that is 
%           relative to the expected values.
%
%     Usage: [prob,t] = goodfitp(O,E,{iter},{doplot})
%
%           O =       column vector of observed values.
%           E =       column vector of corresponding expected values.
%           iter =    number of permutations [default=1000].
%           doplot =  optional boolean flag indicating, if true, that plot(s) 
%                       of distribution and critical value is to be produced 
%                       [default = 0].
%           ------------------------------------------------------------
%           prob =    scaler probability for the agreement of observed with
%                       expected.
%           t =       observed test-statistic value.
%

% RE Strauss, 4/11/96
%   7/1/99 -  renamed from goodfit().
%   9/22/99 - changed usage of null input arguments; added error message for 
%               non-vectors.
%   2/27/01 - return mean-square rather than sum-of-squares; added plot option.

function [prob,t] = goodfitp(O,E,iter,doplot)
  if (nargin < 3) iter = []; end;
  if (nargin < 4) doplot = []; end;

  if (isempty(iter))                    % Set argument to default value
    iter = 1000;                        %   if not passed
  end;
  if (isempty(doplot))
    doplot = 0;
  end;

  if (min(size(O))>1 | min(size(E))>1)
    error('GOODFITP: matrices of observed and expected values must be vectors');
  end;

  O = O(:);
  E = E(:);

  No = length(O);
  Ne = length(E);

  if (No ~= Ne)                         % Check for compatibility
    disp('GOODFITP: vector sizes incompatible');
    return;
  end;

  % Test statistic value for observed data

  t = (O-E)'*(O-E)/No;                  % Mean-squared error

  % Randomized probability

  Bt = t*ones(iter,1);                  % Allocate result vector
  for bit = 2:iter                      % Randomization iterations
    BO = O([randperm(No)]);             % Random permutation of obs values
    Bt(bit) = (BO-E)'*(BO-E)/No;        % Randomized test-statistic value
  end;

  Bt = sort(Bt);                        % Sort randomized stat values
  prob = randprob(t,Bt,doplot);         % Get right-tailed probability
prob

  t = t./length(O);                     % Convert sum-of-squares to mean-square

  return;
