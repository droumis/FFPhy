% RegrCrossValidate: Cross-validates the multiple regression model by jackknifing
%         1, 2, ... m observations (where k < n-2) and tracking the increase in 
%         variance of the residuals of predicted values, as well as the 
%         degredation of the MSE and the adjusted R^2 as the jackknifed predicted
%         residuals are merged with those from the regression.  
%         For jackknifing 2...m observations, the function reports the accumulation 
%         over k randomized iterations or all possible subsets, whichever is fewer.
%
%     Usage: [predstd,njack,b] = regrcrossvalidate(X,y,{maxiter},{maxdel},{noplot})
%
%         X =       [n x p] matrix of independent (predictor) variables.
%         y =       [n x 1] vector of dependent (response) variable.
%         maxiter = optional k, the maximum number of randomized iterations at each 
%                     jackknifing step; if the number of all possible subsets is 
%                     less than this value, uses the subsets instead [default = 100].
%         maxdel =  optional m, the maximum number of observations deleted by 
%                     jackknifing [default = half of the observations = ceil(n/2)].
%         noplot =  optional boolean flag indicating, if true, that plot of predvar
%                     versus njack is not to be produced [default = 0].
%         ---------------------------------------------------------------------------
%         predstd = [m+1 x 1] vector of standard deviations of 0...m predicted 
%                     (jackknifed) residuals, accumulated over k iterations.  The
%                     zeroeth value in the matrix is based on the residuals of all
%                     observations used in fitting the regression model.
%         njack =   [m+1 x 1] vector of number of jackknifed points at each step
%                     (=[0:m]');
%         b =       coefficients from regression of predstd(2:end) on njack(2:end);
%         

% RE Strauss, 11.26.02

function [predstd,njack,b] = regrcrossvalidate(X,y,maxiter,maxdel,noplot)
  if (nargin < 3) maxiter = []; end;
  if (nargin < 4) maxdel = []; end;
  if (nargin < 5) noplot = []; end;
  
  get_b = 0;
  if (nargout>2)
    get_b = 1;
  end;
  
  if (isvector(X))
    X = X(:);
  end;
  [n,p] = size(X);
  y = y(:);
  if (length(y)~=n)
    error('  RegrCrossValidate: input matrices X and y are incompatible.');
  end;
  
  if (isempty(maxiter))
    maxiter = 100;
  end;
  if (isempty(maxdel))
    maxdel = ceil(n/2);
  end;
  if (isempty(noplot))
    noplot = 0;
  end;
  
  if (maxdel > n-2)
    error('  RegrCrossValidate: maxdel too large.');
  end;
  
  predstd = zeros(maxdel+1,1);
  njack = [0:maxdel]';
  
  [b,stats,pred,resid] = linregr(X,y);  % Original regression solution
  predstd(1) = std(resid);
  
  for i = 1:maxdel                      % Cycle thur number of obs to be jackknifed
    sumx = 0;
    sumx2 = 0;
    N = 0;
    
    nc = comb(n,i);                       % Number of possible subsets

    if (nc <= maxiter)                    % If all combinations are feasible,
      c = combvals(n,i);                    % Get combinations    
      for ir = 1:size(c,1);                 % Iterate thru combinations
        ic = c(ir,:);                         % Current combination
        inc = 1:n;
        inc(ic) = [];                         % Complement of combination

        Xp = X(ic,:);                        % Partition into calibration & test sets
        yp = y(ic);
        XX = X(inc,:);
        yy = y(inc);
        [b,stats,pred,resid] = linregr(XX,yy,[],Xp,yp); % Predicted residuals
        
        sumx = sumx + sum(resid);
        sumx2 = sumx2 + sum(resid.*resid);
        N = N + length(resid);
      end;
      
    else                                  % Else use random permutations
      for ir = 1:maxiter
        r = randperm(n);
        X = X(r,:);
        y = y(r);
        
        Xp = X(1:i,:);                      % Partition into calibration & test sets
        yp = y(1:i);
        XX = X(i+1:end,:);
        yy = y(i+1:end);
        [b,stats,pred,resid] = linregr(XX,yy,[],Xp,yp); % Predicted residuals
        
        sumx = sumx + sum(resid);
        sumx2 = sumx2 + sum(resid.*resid);
        N = N + length(resid);
      end;
    end;
    predstd(i+1) = sqrt((sumx2-(sumx*sumx/N))/N);
  end;
  
  if (get_b)
    b = linregr(njack(2:end),predstd(2:end));
  end;
  
  if (~noplot)
    plot(njack,predstd,'k');
    putbnd(njack,predstd);
    putregrline(njack(2:end),predstd(2:end))
    if (maxdel<10)
      puttick(njack);
    end;
    putxlab('Number of jackknifed observations');
    putylab('Std( jackknifed residuals )');
  end;

  return;
  