% MATCOND:  Examines the condition of the correlation matrix of a data matrix 
%           as a function of reducing (rarifying) either variables or 
%           observations.  Uses the Linpack reciprocal condition factor returned 
%           by the Matlab function rcond(): if the correlation matrix is well 
%           conditioned, the condition factor is near 1; if it is badly 
%           conditioned, the condition factor is near 0.  The condition factor 
%           is identical for correlation or covariance matrices.
%
%     Usage: [v,cond,condci,bestcomb,ncomb] = ...
%                             matcond(X,{do_obs},{maxcomb},{noplot},{CI_level})
%
%           X =         [n x p] data matrix.
%           do_obs =     optional boolean flag indicating, if true, that observations 
%                         are to be reduced [default = 0 = variables reduced].
%           maxcomb =   optional maximum number of combinations to be examined at 
%                         each level of variables (or observations) [default = 1000].  
%                         If the exact number needed is greater than this, a random 
%                         sample of combinations of this size is used.
%           noplot =    optional boolean flag indicating, if true, that a plot 
%                         is not to be produced [default = 0].
%           CI_level =  optional level of confidence intervals for condition 
%                         factors [default = 0.95].
%           ------------------------------------------------------------------------
%           v =         column vector [length q] of values of rarified variables 
%                         (or observations).
%           cond =      [q x 1] vector of corresponding condition factors.
%           condci =    [q x 2] matrix of confidence intervals of condition factors. 
%           bestcomb =  [q x p] (or [q x n]) matrix of particular combinations of 
%                         variables (or observations) giving greatest condition 
%                         factors, from the set of min(all,maxcomb) combinations 
%                         examined at each rarification level.  Gauranteed to be 
%                         globally optimal only if all combinations are examined.
%           ncomb =     corresponding vector of numbers of combinations for each 
%                         value of rarified variables (or observations).
%

% RE Strauss, 10/4/00

function [v,cond,condci,bestcomb,ncomb] = matcond(X,do_obs,maxcomb,noplot,CI_level)
  if (nargin < 2) do_obs = []; end;
  if (nargin < 3) maxcomb = []; end;
  if (nargin < 4) noplot = []; end;
  if (nargin < 5) CI_level = []; end;

  if (isempty(do_obs))
    do_obs = 0;
  end;
  if (isempty(maxcomb))
    maxcomb = 1000;
  end;
  if (isempty(noplot))
    noplot = 0;
  end;
  if (isempty(CI_level))
    CI_level = 0.95;
  end;
  if (CI_level > 1)
    CI_level = CI_level/100;
  end;

  [n,p] = size(X);

  if (~do_obs)                           % Rarify variables
    v = [3:p]';                           % Initialize output and working matrices
    cond = zeros(p-2,1);
    condci = zeros(p-2,2);
    bestcomb = zeros(p-2,p);
    rc = zeros(maxcomb,1);

    for pc = 3:p                          % Vary number of variables
      ncombs = comb(p,pc);                  % Number of combinations
      if (ncombs <= maxcomb)                % Examine all combinations
        combs = randcomb(p,pc);
      else                                  % Or examine random combinations
        combs = randcomb(p,pc,maxcomb,1);   %   without replacement
        ncombs = maxcomb;
      end;

      maxrc = 0;
      for i = 1:ncombs                      % For all combinations
        x = X(:,combs(i,:));                  % Subsample variables
        rc(i) = rcond(corrcoef(x));           % Stash condition factor
        if (rc(i) > maxrc)
          maxrc = rc(i);
          bestcomb(pc-2,1:pc) = combs(i,:);
        end;
      end;
      
      cond(pc-2) = mean(rc(1:ncombs));
      if (ncombs > 3)
        ci = bootci(rc(1:ncombs),[],[],CI_level);
        condci(pc-2,:) = ci';
      else
        condci(pc-2,:) = rc([1,ncombs])';
      end;
    end;

    ncomb = [comb(p,3:p)]';                 % Numbers of combinations

    if (~noplot)
      figure;
      plot(v,cond,'k',v,condci(:,1),'k--',v,condci(:,2),'k--');
      putbnd([v;v;v],[cond;condci(:)]);
      putxlab('Number of variables');
      putylab('Matrix condition');
    end;
  end;

  if (do_obs)                            % Rarify observations
    v = [3:n]';                           % Initialize output and working matrices
    cond = zeros(n-2,1);
    condci = zeros(n-2,2);
    bestcomb = zeros(n-2,n);
    rc = zeros(maxcomb,1);

    for nc = 3:n                          % Vary number of observations
      ncombs = comb(n,nc);                  % Number of combinations
      if (ncombs <= maxcomb)                % Examine all combinations
        combs = randcomb(n,nc);
      else                                  % Or examine random combinations
        combs = randcomb(n,nc,maxcomb);
        ncombs = maxcomb;
      end;

      maxrc = 0;
      for i = 1:ncombs                      % For all combinations
        x = X(combs(i,:),:);                  % Subsample observations
        rc(i) = rcond(corrcoef(x));           % Stash condition factor
        if (rc(i) > maxrc)
          maxrc = rc(i);
          bestcomb(nc-2,1:nc) = combs(i,:);
        end;
      end;
      
      cond(nc-2) = mean(rc(1:ncombs));
      if (ncombs > 3)
        ci = bootci(rc(1:ncombs),[],[],CI_level);
        condci(nc-2,:) = ci';
      else
        condci(nc-2,:) = rc([1,ncombs])';
      end;
    end;

    ncomb = [comb(n,3:n)]';                 % Numbers of combinations

    if (~noplot)
      figure;
      plot(v,cond,'k',v,condci(:,1),'k--',v,condci(:,2),'k--');
      putbnd([v;v;v],[cond;condci(:)]);
      putxlab('Number of observations');
      putylab('Matrix condition');
    end;
  end;

  return;
