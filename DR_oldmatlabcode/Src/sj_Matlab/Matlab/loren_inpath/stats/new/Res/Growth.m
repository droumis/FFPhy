% GROWTH: Fits one or more longitudinal (size-at-age) growth models.  
%         See Wan, Zhong & Wang 1998 for description and notation.
%         See function growthfn() for parameters of models.
%
%     Usage: [param,R2adj,pred_size] = growth(age,size,{model},{age_eval},{doplot})
%
%         age =       vector of ages.
%         size =      corresponding vector of sizes.
%         model =     boolean vector indicating the 'm' models to be fitted to the data:
%                       pos 1: von Bertalanffy, = monomolecular (3 parameters);
%                           2: Gompertz   (3 parameters)
%                           3: logistic   (3 parameters)
%                           4: Richards   (4 parameters)
%                           5: Janoschek  (4 parameters)
%                           6: Hill       (4 parameters)
%                           7: Wan        (4 parameters)
%                           8: France     (5 parameters)
%                       [default = all].
%         age_eval =  optional vector of ages at which the growth function is to be 
%                       evaluated [default = age if 'pred_size' is requested, null 
%                       otherwise].
%         doplot =    optional boolean flag indicating that a plot of data and fitted 
%                       growth functions is to be produced [default = 0].
%         -----------------------------------------------------------------------------------
%         param =     [m x p] matrix of estimated parameters for the 'm' models requested.  
%                       The number of columns p, which varies from 3-5, is the maximum number 
%                       of parameters for the models requested; the latter columns are set 
%                       to NaN if not used.
%         R2adj =     [m x 1] vector of adjusted R^2 goodness-of-fit values.
%         pred_size = [m x n] matrix of predicted size values corresponding to 'n' age 
%                       values, for each of 'm' models [default = predicted values for 
%                       the input age vector].
%

% RE Strauss, 10/21/98
%   9/3/99 -  changed plot colors for Matlab v5.
%   1/4/00 -  changed fminu() to fmins().
%   5/4/00 -  changed fmins() to fminsearch();
%             set default for models to 'all'.

function [param,R2adj,pred_size] = growth(age,siz,model,age_eval,doplot)
  if (nargin < 3) model = []; end;
  if (nargin < 4) age_eval = []; end;
  if (nargin < 5) doplot = []; end;

  nmodel = 8;
  nparam = [3 3 3 4 4 4 4 5];

  get_R2adj = 0;
  get_pred_size = 0;
  if (nargout > 1)
    get_R2adj = 1;
  end;
  if (nargout > 2)
    get_pred_size = 1;
  end;

  if (min([size(age) size(siz)])~=1 | length(age)~=length(siz))
    error('  GROWTH: age and size matrices must be corresponding vectors.');
  end;

  if (isempty(age_eval) & get_pred_size)
    age_eval = uniquef(age)';
  end;

  if (isempty(model))
    model = ones(1,nmodel);
  end;

  if (isempty(doplot))
    doplot = 0;
  end;

  if (doplot & ~get_pred_size)
    get_pred_size = 1;
    age_eval = linspace(min(age),max(age));
  end;

  if (length(model) < nmodel)           % Extend model vector if too short
    model = [model zeros(1,nmodel-length(model))];
  end;

  nmodel = sum(model>0);
  maxparam = max(nparam(model==1));

  param = zeros(nmodel,maxparam);       % Allocate output matrices
  if (get_R2adj)
    R2adj = zeros(nmodel,1);
  end;
  if (get_pred_size)
    len_age_eval = length(age_eval);
    pred_size = zeros(nmodel,len_age_eval);
  end;

  nm = 0;
  for m = 1:length(model)               % Cycle thru models
    if (model(m))
current_model = model(m)
      nm = nm+1;
      np = nparam(m);
      p = ones(1,np);
      if (isin(m,[3:6,8]))
        p(3) = max(siz);
      end;

      p = fminsearch('growthfn',p,optimset('Display','off'),m,age,siz); % Optimization
      param(nm,1:np) = abs(p);

      if (get_R2adj)
        sse = growthfn(p,m,age,siz);
        n = length(age);
        ssto = sum((siz-mean(siz)).^2);
        R2adj(nm) = 1 - (((n-1)/(n-np))*(sse/ssto));
      end;
      if (get_pred_size)
        [x,predsize] = growthfn(p,m,age_eval,age_eval);
        if (size(predsize,1)>1)
          predsize = pred_size';
        end;
        pred_size(nm,:) = predsize;
      end;
    end;
  end;

  if (doplot)
    plot(age,siz,'ok');
    hold on;
    for i = 1:size(pred_size,1)
      plot(age_eval,pred_size(i,:),'k');
    end;
    hold off;

    xmin = min([0 min(age) min(age_eval)]);
    xmax = max([max(age) max(age_eval)]);
    ymin = min([min(pred_size) min(siz)]);
    ymax = max([max(pred_size) max(siz)]);

    delta = 0.05;
    deltax = delta*(xmax-xmin);
    deltay = delta*(ymax-ymin);
    axis([xmin-deltax xmax+deltax ymin-deltay ymax+deltay]);

    putxlab('Age');
    putylab('Size');
  end;

  return;
