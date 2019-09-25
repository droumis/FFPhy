% Trajectory: Fit mean (1-parameter), linear (2-param), quadratic (3-param), 
%             bilinear (4-param), and transitional bilinear (5-param) regression 
%             models to trajectory data, and determine model of best fit by 
%             three criteria: minimum MSE, maximum R2a, and significant 
%             improvement-of-fit F-test.
%
%     Usage: [param,se,R2a,fit,iof,inflect] = trajectory(x,y,{noplot})
%
%         x =       vector of independent-variable values. 
%         y =       corresponding vector of dependent-variable values.
%         noplot =  optional boolean variable indicating that plots of data and 
%                     fitted regression models are not to be produced 
%                     [default = 0].
%         ---------------------------------------------------------------------
%         param =   [5 x 5] matrix of model parameters, where row i contains the 
%                     i parameters of the ith model.
%         se =      [5 x 2] vector of mean squared errors and sum of squared 
%                     errors for the five models: [mse, sse].
%         R2a =     [5 x 1] vector of coefficients of determination, adjusted 
%                     for the number of parameters in the model.
%         fit =     [5 x 4] matrix of F statistics for testing the fit of 
%                     model i: [prob, F, df].
%         iof =     [5 x 4] matrix of improvement-of-fit F statistics for
%                     testing the fit of model i over i-1, given the increase 
%                     in the number of parameters: [prob, F, df].
%         inflect = [5 x 1] matrix of sizes (x) at which inflection occurs.
%

% Need to bootstrap all models simultaneously for parameters and F-tests.

% RE Strauss, 10/9/00
%   11/18/00 - added F-test for overall fits of models.
%   12/20/00 - return inflection points for nonlinear trajectories.
%   1/2/01 -   use r2adj() to return adjusted R^2 for 3-5 param models.
%   1/3/01 -   use fiof() for improvement-of-fit F-test.

function [param,se,R2a,fit,iof,inflect] = trajectory(x,y,noplot)
  if (nargin < 3) noplot = []; end;

  if (isempty(noplot))
    noplot = 0;
  end;

  if (~isvector(x) | ~isvector(y))
    error(  'TRAJECTORY: input must be vectors.');
  end;
  if (length(x) ~= length(y))
    error(  'TRAJECTORY: input must be vectors of identical length.');
  end;

  x = x(:);                               % Convert to col vectors
  y = y(:);
  n = length(x);

  xmin = min(x);
  xmax = max(x);
  msto = var(y);
  ssto = msto*(n-1);

  param = zeros(5,5);
  sse = zeros(5,1);
  mse = zeros(5,1);
  R2a = zeros(5,1);
  F_fit = NaN*ones(5,1);                % Overall fit F-statistics
  df_fit = NaN*ones(5,2);
  prob_fit = NaN*ones(5,1);
  F_iof = NaN*ones(5,1);                % Improvement of fit F-statistics
  df_iof = NaN*ones(5,2);
  prob_iof = NaN*ones(5,1);
  inflect = NaN*ones(5,1);

  xlin = linspace(xmin,xmax,1000);      % Abscissa values for estimating inflection


  % 1-parameter model (mean)

  ymean = mean(y);
  param(1,1) = ymean;
  e = y - ymean;
  sse(1) = e'*e;
  mse(1) = sse(1)/(n-1);
  R2a(1) = 0;

  if (~noplot)
    scatter(x,y);
    hold on;
    plot([xmin xmax],[ymean ymean],'k');
    hold off;
  end;


  % 2-parameter model (linear)

  p = 2;
  [b,stats,pred] = linregr(x,y,[],[xmin xmax]);

  param(p,1:p) = b';
  sse(p) = stats(2)*(n-p);
  mse(p) = stats(2);
  R2a(p) = max([0, stats(1)]);
  F_fit(p) = stats(3);
  df_fit(p,:) = stats(4:5)';
  prob_fit(p) = stats(6);
  [F_iof(p),prob_iof(p),df_iof(p,:)] = fiof(sse(p-1),sse(p),p-1,p);

  if (~noplot)
    scatter(x,y);
    hold on;
    plot([xmin xmax],pred,'k');
    hold off;
  end;


  % 3-parameter model (quadratic)

  p = 3;
  [b,msep] = quadratic(x,y,~noplot);

  param(p,1:p) = b';
  sse(p) = msep*(n-p);
  mse(p) = msep;
  [R2a(p),F_fit(p),pr_fit(p),df_fit(p,:)] = r2adj(ssto-sse(p),ssto,n,p);
[sse(p-1),sse(p),p-1,p]  
  [F_iof(p),prob_iof(p),df_iof(p,:)] = fiof(sse(p-1),sse(p),p-1,p);

  ylin = b(1) + b(2)*xlin + b(3)*(xlin.^2);
  if (b(3) < 0)
    [xlinmax,yinfl] = max(ylin);
  else
    [xlinmin,yinfl] = min(ylin);
  end;
  inflect(p) = xlin(yinfl);

  % 4-parameter model (bilinear)

  p = 4;
  [b,t,msep] = bilinear(x,y,~noplot);

  param(p,1:p) = b;
  sse(p) = msep*(n-p);
  mse(p) = msep;
  [R2a(p),F_fit(p),pr_fit(p),df_fit(p,:)] = r2adj(ssto-sse(p),ssto,n,p);
  [F_iof(p),prob_iof(p),df_iof(p,:)] = fiof(sse(p-1),sse(p),p-1,p);
  inflect(p) = t;

  % 5-parameter model (transitional bilinear)

  p = 5;
  [b,msep] = bilintrans(x,y,~noplot);

  param(p,1:p) = b;
  sse(p) = min([msep*(n-p), sse(p-1)]);
  mse(p) = msep;
  [R2a(p),F_fit(p),pr_fit(p),df_fit(p,:)] = r2adj(ssto-sse(p),ssto,n,p);
  [F_iof(p),prob_iof(p),df_iof(p,:)] = fiof(sse(p-1),sse(p),p-1,p);
  inflect(p) = b(4);


  % Condense parameters

  se = [mse sse];
  fit = [prob_fit F_fit df_fit];
  iof = [prob_iof F_iof df_iof];

  return;
