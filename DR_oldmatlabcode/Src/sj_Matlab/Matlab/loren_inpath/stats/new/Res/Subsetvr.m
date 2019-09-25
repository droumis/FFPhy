% SUBSETVR: Finds the 'best' subset of r<p variables that maximizes the 
%           reciprocal condition number (rcond) of the matrix.  Tests all 
%           possible combinations of r variables if possible; if not, 
%           uses a stepwise procedure.
%
%     Usage: [vars,condnum,allcomb] = subsetvr(X,{r},{maxcomb})
%         
%           C =       [p x p] square symmetric matrix.
%           r =       desired number of variables.  If omitted, searches for best 
%                       subsets for r = (p-1) to 2.
%           maxcomb = maximum threshold for all possible combinations to try 
%                       [default = 5000].
%           ---------------------------------------------------------------------
%           vars =    [1 x r] list of var indices for single solution (if r is 
%                       given), or [p-2 x r] matrix of lists of vars for stepwise 
%                       solutions.
%           condnum = scalar (if r is given) or [p-2 x 1] vector of condition 
%                       numbers for solutions.  If the matrix is well 
%                       conditoned, condnum is near 1; if badly conditioned, 
%                       condnum is near 0.
%           allcomb = scalar (if r is given) or [p-2 x 1] vector of boolean flags 
%                       indicating whether solution was obtained by stepwise 
%                       procedure (=1) or by examining all possible combinations 
%                       of variables (=0).
%

% RE Strauss, 4/24/99; modified from 'steprank'.
%   9/3/99 - changed plot colors for Matlab v5.

function [vars,condnum,allcomb] = subsetvr(C,r,maxcomb)
  if (nargin < 2) r = []; end;
  if (nargin < 3) maxcomb = []; end;

  [p,p1] = size(C);

  err = 0;
  if (p<3 | p~=p1)
    err = 1;
  elseif (any(any(abs(C-C')>eps)))
    err = 1;
  end;
  if (err)
    error('SUBSETVR: input matrix must be square symmetric and of minimum order 3');
  end;

  if (isempty(r))
    r = (p-1):-1:2;
  end;
  if (isempty(maxcomb))
    maxcomb = 5000;
  end;

  vars = zeros(length(r),max(r)); % Allocate output matrices
  condnum = zeros(length(r),1);
  allcomb = zeros(length(r),1);

  for ir = 1:length(r)              % Cycle thru solution lengths
    nc = comb(p,r(ir));             % Number of possible combinations
    if (nc <= maxcomb)              % Examine all possible combinations of vars
      allcomb(ir) = 1;
      comblist = combvals(p,r(ir));   % List of all possible combinations
      best_rcond = 0;
      best_vars = [];

      for ic = 1:nc                   % Find combination giving max condition
        v = comblist(ic,:);
        bc = rcond(C(v,v));
        if (bc > best_rcond)
          best_rcond = bc;
          best_vars = v;
        end;
      end;
      
      vars(ir,1:length(v)) = best_vars; % Stash best solution
      condnum(ir) = best_rcond;

    else                            % Use stepwise procedure
      best_rcond = 0;
      v = [1:p];

      for i = 1:(p-1)                   % Find two best vars
        for j = (i+1):p
          bc = rcond(C([i j],[i j]));
          if (bc > best_rcond)
            best_rcond = bc;
            isave = i;
            jsave = j;
          end;
        end;
      end;

      vlist = [isave jsave];              % Initialize lists of variables
      v(vlist) = [];

      if (r(ir) > 2)                      % Find vars 3 ... r
        for irr = 3:r(ir)
          best_rcond = 0;
          for i = 1:length(v)
            j = [vlist v(i)];
            bc = rcond(C(j,j));
            if (bc > best_rcond)
              best_rcond = bc;
              isave = i;
            end;
          end;
          vlist = [vlist v(isave)];
          v(isave) = [];
        end;
      end;

      vars(ir,1:length(vlist)) = sort(vlist); % Stash best solution
      condnum(ir) = best_rcond;

    end;  % if (nc <= maxcomb)
  end;  % for ir = 1:length(r)

  if (length(r)>1)
    nvars = sum(vars'>0)';
    plot(nvars,log(condnum)./nvars,'k');
    putbnd(nvars,log(condnum)./nvars);
    putxlab('Number of variables');
    putylab('Relative matrix condition');
  end;

  return;