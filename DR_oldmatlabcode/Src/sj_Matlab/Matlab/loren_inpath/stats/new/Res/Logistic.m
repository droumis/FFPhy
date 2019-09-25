% LOGISTIC: Least-squares fitting of logistic regression model for continuous 
%           Y = Ya/(1+exp(-k*(t-h))) = Ya/(1+exp((2.2/c)*(h-t))).
%
%             Ya = asymptotic function maximum
%             k =  rate constant
%             c =  time spread between 0.1*Ya and 0.9*Ya
%             h =  half-saturation time
%
%     Syntax:  [param,mse,ci,Ypred] = logistic(t,Y,{tpred},{iter},{doplot})
%
%         t =         vector of abscissa values.
%         Y =         corresponding vector of ordinate (size) values.
%         tpred =     optional vector of t values for which predicted Y 
%                       values are desired.  Pass as [] if not used.
%         iter = number of bootstrap iterations for confidence intervals 
%                       on parameters [default = 0].
%         doplot =    optional boolean flag indicating that plot of data and 
%                       fitted function is to be produced [default = 0].
%         --------------------------------------------------------------------
%         param =     row vector of parameter estimates: [k,h,Y0,Ya].
%         mse =       mean-squared-residual.
%         ci =        95% bootstrapped confidence intervals of parameters, if 
%                       iter > 1.
%         Ypred =     vector of predicted Y values corresponding to tpred (if 
%                        given) or to Y (if not).
%

% RE Strauss, 4/7/99, modified from previous function.
%   9/3/99 -  convert input row vectors to col vectors; 
%             changed plot colors for Matlab v5.
%   1/4/00 -  changed fminu() to fmins().

function [param,mse,ci,Ypred] = logistic(t,Y,tpred,iter,doplot)
  if (nargin < 3) tpred = []; end;
  if (nargin < 4) iter = []; end;
  if (nargin < 5) doplot = []; end;

  get_mse = 0;
  get_ci = 0;
  get_Ypred = 0;

  if (nargout >= 2)
    get_mse = 1;
  end;
  if (nargout >= 3)
    get_ci = 1;
  end;
  if (nargout >= 4)
    get_Ypred = 1;
  end;

  if (isempty(tpred))                         % Default input arguments
    tpred = t;
  end;
  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(doplot))
    doplot = 0;
  end;

  if (size(t,1)==1)                           % Convert row vectors to col vectors
    t = t';
  end;
  if (size(Y,1)==1)
    Y = Y';
  end;
  if (size(t) ~= size(Y))
    error('LOGISTIC: input data vectors not compatible');
  end;

  k = 1;                                      % Initial parameter estimates
  h = mean(t);                    
  Ya = max(Y);

  param = fmins('logistf',[k h Ya],[],[],t,Y);    % Fit model

  if (get_mse)                                % Mean squared error
    mse = logistf(param,t,Y);
  end;

  k =  param(1);
  h =  param(2);
  Ya = param(3);

  if (get_Ypred)                              % Predicted values
    Ypred = Ya./(1+exp(-k*(tpred-h)));
  end;

  if (iter>0 & get_ci)                   % Bootstrapped confidence intervals
    procs = [1 0 0 0];
    alpha = 0.05;
    ci = bootstrp('logistbt',procs,iter,alpha,[t,Y])
  end;

  if (doplot)                                 % Plot
    tval = linspace(min(t),max(t));
    Yval = Ya./(1+exp(-k*(tval-h)));
    hval = Ya./(1+exp(0));

    figure;
    plot(t,Y,'ok');
    putbnd(t,Y);
    hold on;
    plot(tval,Yval,'k');
    hold off;

    figure;
    plot(t,Y,'ok');
    putbnd(t,Y);
    hold on;
    plot(tval,Yval,'k');
    plot([h h],[0 hval],':k');
    plot([min(t) h],[hval hval],':k');
    plot([min(t) max(t)],[Ya Ya],':k');
    hold off;
  end;

  return;
