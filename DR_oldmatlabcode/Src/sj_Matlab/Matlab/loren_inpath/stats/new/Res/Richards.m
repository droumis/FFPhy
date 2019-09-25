% RICHARDS: Fits the Richards (1959) 4-parameter growth model, using Fletcher's
%           (1979) integrated expression, as modified by Brisbin et al (1986):
%
%             W(t) = [Wa^c + (Wa^c-W0^c) * exp(-2t/T * (m+1))]^(1/c)
%
%           where Wa is asymptotic size, W0 is initial size at t=0, T is the 
%           total growth period, m (>0) is the Richards curve-shape parameter, 
%           and c=1-m.  Wa, W0 or both may be set to fixed values for reduced 
%           models.
%
%     Usage: [param,ci,mse,Wpred,infl] = 
%               richards(t,W,{tpred},{doplot},{paramfix},{boot_iter},{ci_level})
%
%           t =         vector of times.
%           W =         corresponding vector of sizes.
%           tpred =     optional vector of times for which predicted values 
%                         are to be generated [default = t].
%           doplot =    optional boolean flag indicating that plot of data and 
%                         function is to be produced [default = 0].
%           paramfix =  optional vector of parameter values [W0,Wa], 
%                         either of which may be held to a fixed value for the 
%                         fitting of a reduced model. Set parameter values to 
%                         NaN if they are to be fitted, and provide fixed values.
%                         Default = [NaN NaN].
%           boot_iter = optional number of bootstrap iterations for confidence 
%                         intervals on parameters [default = 0].
%           ci_level =  optional confidence value for confidence intervals 
%                         [default = 95].
%           ----------------------------------------------------------------------
%           param =     4-element row vector of parameter values [W0,Wa,T,m].
%           ci =        [2 x 4] matrix of lower and upper alpha-level confidence 
%                         limits on parameters, asymptotic (if boot_iter = 0) or 
%                         bootstrapped.
%           mse =       mean squared error.
%           Wpred =     predicted values of W corresponding to tpred.
%           infl =      3-element row vector of characteristics of inflection 
%                         point [tinfl,Winfl,maxslope]: 
%                           tinfl = time at inflection;
%                           Winfl = size at inflection; 
%                           maxslope: instantaneous slope at inflection.
%

% Richards, F.J. 1959.  A flexible growth function for empirical use.  
%   J.Exp.Botany 10:290-300.
% Fletcher, R.I. 1975.  A general solution for the complete Richards function.
%   Math.Biosci. 27:349-360.
% Brisbin, I.L. Jr., G.C. White, P.B. Bush. 1986.  PCB intake and the growth of 
%   waterfowl: multivariate analyses based on a reparameterized Richards sigmoid 
%   model.  Growth 50:1-11.

% RE Strauss, 3/29/99

function [param,ci,mse,Wpred,infl] = ...
            richards(t,W,tpred,doplot,paramfix,boot_iter,ci_level)

  if (nargin < 3) tpred = []; end;
  if (nargin < 4) doplot = []; end;
  if (nargin < 5) paramfix = []; end;
  if (nargin < 6) boot_iter = []; end;
  if (nargin < 7) ci_level = []; end;

  get_ci = 0;
  get_mse = 0;
  get_Wpred = 0;
  get_infl = 0;

  if (nargout >= 2) get_ci = 1; end;
  if (nargout >= 3) get_mse = 1; end;
  if (nargout >= 4) get_Wpred = 1; end;
  if (nargout >= 5) get_infl = 1; end;

  errflag = 0;                                % Check input vectors
  if (min(size(t))>1 | min(size(W))>1)
    errflag = 1;
  end;
  N = length(t);
  if (length(W) ~= N)
    errflag = 1;
  end;
  if (errflag)
    error('RICHARDS: t and W must be corresponding vectors');
  end;

  if (size(t,2) > 1)                          % Transpose input vectors to cols
    t = t';
  end;
  if (size(W,2) > 1)
    W = W';
  end;

  if (isempty(tpred))                         % Default parameter values
    tpred = t;
  end;
  if (isempty(doplot))
    doplot = 0;
  end;
  if (isempty(boot_iter))
    boot_iter = 0;
  end;
  
  if (isempty(paramfix))
    paramfix = [NaN NaN];
  elseif (length(paramfix)==1)
    paramfix = [paramfix NaN];
  end;

  if (isempty(ci_level))
    ci_level = .95;
  elseif (ci_level > 1)
    ci_level = ci_level/100;
  end;

  if (doplot)
    get_Wpred = 1;
    get_infl = 1;
  end;



  return;
