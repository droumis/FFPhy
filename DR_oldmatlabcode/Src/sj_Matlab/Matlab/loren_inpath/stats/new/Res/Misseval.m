% MISSEVAL: Evaluate performance of missing-data estimation methods (EM and PCA) 
%           by randomly introducing missing values into a complete data matrix, 
%           and assessing the difference between the original and predicted 
%           values.  Criterion of performance is the mean squared error (MSE) 
%           from the regression of predicted missing values on original values, 
%           accumulated across iterations.
%
%     Usage: [goodresults,allbouts] = misseval(X,type,{pmiss},{iter},{nobs},{nchars})
%
%         X =           complete [n x p] data matrix.
%         type =        estimation method: 0 for EM, 1 for PCA.
%         pmiss =       optional vector of proportions (0-1) or numbers (>1) of 
%                         missing data to be randomly introduced into data matrix
%                         [default = 0.01:0.01:0.50].
%         iter =        vector (length 3) containing optional numbers of 
%                         randomization iterations for each value of 'pmiss': 
%                         [total iterations, iterations per bout, maximum bouts]
%                         where number of bouts = (total iterations)/(iterations per 
%                         bout) but max number of bouts can be greater
%                         (default = [1000,200,10]).
%         nobs =        optional number of randomized observations to be included 
%                         in matrix having missing data [default = all].
%         nchars =      optional number of randomized variables to be included 
%                         in matrix having missing data [default = all].
%         --------------------------------------------------------------------------
%         goodresults = matrix containing results from bouts in which all 
%                         iterations were successful.  Also written to 
%                         'MissevalGoodResults.txt' (see below).
%         allbouts =    3-column matrix containing numbers of successful iterations
%                         per bout: [cur_pmiss, bout_number, percent_success]
%                         (see below).  Also written to 'MissevalAllBouts.txt'.
%
%     Output stored in 'goodresults',and written to 'MissevalGoodResults':
%         cur_pmiss =   current value of pmiss.
%         bout_number = current bout within current value of pmiss.
%         mse =         residual variance for regression of predicted missing 
%                         values versus original.
%         deltavar =    signed proportional change in total variance from original 
%                         to predicted missing values: [mean, stderr] 
%                         across iterations.
%         nonuni =      [1 x 4] vector of mean-squared differences between observed and 
%                       expected marginal totals.
%                         1) relative difference for matrix;
%                         2) relative difference for rows;
%                         3) relative difference for columns;
%                         4) contingency MAD statistic;
%

% RE Strauss, 3/17/99
%   10/2/00 -  added time limit to call to missem().
%   10/3/00 -  allow proportion or number of missing values.
%   10/6/00 -  upgraded standardization and documentation;
%                added doplot option.
%   10/10/00 - allowed for 'pmiss' to be a vector rather than scalar.
%   2/2/01 -   removed 'do_deltavar' flag from input arguments.
%   3/6/01 -   move tofile to within 1:np loop.
%   2/28/02 -  general overhaul; automated many functions, provided for 
%                iterations within bouts.
%   4/12/02 -  change default value for 'iter'.
%   10/30/02 - changed misspca() to misspc().

function [goodresults,allbouts] = misseval(X,type,pmiss,iter,nobs,nchars)
  if (nargin < 3) pmiss = []; end;
  if (nargin < 4) iter = []; end;
  if (nargin < 5) nobs = []; end;
  if (nargin < 6) nchars = []; end;
  
  maxtime = 30*60;                      % 30 minute time limit per iteration
  
  if (sum(~isfinite(X(:))>0))
    error('  MISSEVAL: input data matrix contains missing data');
  end;

  [r,c] = size(X);

  if (isempty(pmiss))
    pmiss = [0.01:0.01:0.50];
  end;
  if (isempty(iter))
    iter = [1000,200,10];
  end;
  if (isempty(nobs))
    nobs = r;
  end;
  if (isempty(nchars))
    nchars = c;
  end;

  if (length(iter)~=3)
    error('  MISSEVAL: iterations vector must contain 3 elements.');
  end;

  npmiss = length(pmiss);
  nmiss = pmiss;
  for i = 1:npmiss
    if (pmiss(i)<1)
      nmiss(i) = round(pmiss(i)*nobs*nchars);
    end;
  end;

  if (max(max(X)) > 10)                 % Log-transform if needed
    X = log(X);
  end;

  max_total_iterations = iter(1);       % Isolate iteration information
  iter_per_bout = iter(2);
  max_bouts = iter(3);
  nbouts = round(max_total_iterations / iter_per_bout);

  result_cols = 10;
  result_rows = npmiss*max_bouts;
  results = zeros(result_rows,result_cols);
  cur_res = 0;
  
  oldX = X;

  for np = 1:npmiss                     % Cycle through proportions of missing values
    cum_oldX = zeros(nmiss(np)*iter_per_bout,1);
    cum_newX = zeros(nmiss(np)*iter_per_bout,1);
    dvar = zeros(iter_per_bout,1);

    timeout = 0;
    breakflag = 0;
    failflag = 0;
    tot_iter = 0;
    in = 0;
    
    all_bouts_failed = 1;
    
    for bout_number = 1:max_bouts           % Cycle thru bouts
      percent_success = 0;
      nonunival = zeros(iter_per_bout,4);
      
      for it = 1:iter_per_bout              % Randomization iterations
        X = oldX;
  
        if (nobs < r)
          robs = randperm(r);
          X = X(robs(1:nobs),:);
        end;
        if (nchars < c)
          rchars = randperm(c);
          X = X(:,rchars(1:nchars));
        end;

        prevX = X;
        [X,failflag] = randmiss(X,nmiss(np));  % Introduce random missing values
        
        if (~failflag)
          nonunival(it,:) = missrandmeas(X);    % Measure nonuniformity of missing values

          [rn,cn] = find(~finite(X));           % Save addresses of missing values

          if (type==0)
            newX = missem(X,maxtime);
            [newX,M,C,propmiss,failflag,stoptime] = missem(X,maxtime);
          else
            [newX,propmiss,failflag] = misspc(X);
          end;
        end;

        if (failflag)
          save missevalmat;
          break;
        end;
        percent_success = percent_success+(100/iter_per_bout);

        diff = newX - prevX;
        nonzero = diff(diff~=0);

        oldvar = sum(diag(cov(prevX)));
        newvar = sum(diag(cov(newX)));
        dvar(it) = (newvar-oldvar)/oldvar;
 
        for i = 1:nmiss(np)
          in = in+1;
          cum_oldX(in) = oldX(rn(i),cn(i));
          cum_newX(in) = newX(rn(i),cn(i));
        end;
      end; % for it = 1:iter
      
      if (~failflag)                            % Results for this bout
        tot_iter = tot_iter + iter_per_bout;
        cum_nonzero = cum_newX - cum_oldX;
  
        [b,stats] = linregr(cum_oldX(:),cum_newX(:));
        mse = stats(2);
        m = mean(dvar);
        s = std(dvar);
        deltavar = [m, s];
        nonuni = mean(nonunival);
        all_bouts_failed = 0;
      else
        mse = NaN;
        deltavar = [NaN, NaN];
        nonuni = [NaN, NaN, NaN, NaN];
      end;
      
      pmiss_bout_percsuccess_totiter = ...
        [pmiss(np) bout_number percent_success tot_iter]      

      cur_res = cur_res+1;
      res = [pmiss(np) bout_number percent_success mse deltavar nonuni];
      results(cur_res,:) = res;
    
      save missevalmat;
      
      if (tot_iter >= max_total_iterations)
        break;
      end;
    end; % for bout_number = 1:max_bouts
    
    if (all_bouts_failed)
      break;
    end;
  end;  % for np = 1:npmiss

  allbouts = results(1:cur_res,1:3);
  tofile(allbouts,'MissevalAllBouts.txt',-5);
  
  i = find(isfinite(results(1:cur_res,4)));
  goodresults = results(i,[1,2,4:end]);
  tofile(goodresults,'MissevalGoodResults.txt',-5);

  return;




