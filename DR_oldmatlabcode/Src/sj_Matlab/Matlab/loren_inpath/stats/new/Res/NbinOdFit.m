% NBINODFIT: Finds the least-squares fit of a the negative binomial distribution 
%            (using the mean-variance [m,v] parameterization) to a set of observed
%            frequencies.
%
%     Usage: [fxpred,m,v] = nbinodfit(fx,{stat},{doplot})
%
%         fx =      vector of absolute frequencies (counts) corresponding to
%                     X=0:(length(fx)-1).
%         stat =    optional value indicating the optimization statistic to be 
%                     minimized in the fitting:
%                       0 : sum of squared deviations [default].
%                       1 : chi-squared deviations.
%                       2 : sum of absolute deviations.
%         doplot =  optional boolean variable indicating, if true, that a plot
%                     of the squared-deviation optimization surface is to be produced
%                     [default = 0].
%         ---------------------------------------------------------------------------
%         fxpred =  predicted values of f(x) from negative binomial distribution.
%         m =       mean parameter.
%         v =       variance parameter.
%

% Hilborn, R. & M. Mangel. 1997. The ecological detective: confronting models with data.  
%   Princeton Univ. Press.

% RE Strauss, 2/19/03
%   3/4/03 -  added chi-squared and absolute deviations.
%   4/19/03 - change plottype for plotsurface() due to change in that function.

function [fxpred,m,v] = nbinodfit(fx,stat,doplot)
  if (~nargin) help nbinodfit; return; end;
  
  if (nargin < 2) stat = []; end;
  if (nargin < 3) doplot = []; end;
  
  if (isempty(stat))   stat = 0; end;
  if (isempty(doplot)) doplot = 0; end;
  
  fx = fx(:);
  x = 0:length(fx)-1;
  [m_init,v_init] = meanwt(x,fx);                 % Initial parameter estimates
  param = [m_init,v_init];
  
  param = fminsearch('nbinodfitf',param,[],fx,stat);   % Optimization
  m = param(1);
  v = param(2);
  fxpred = nbinodcdf(x,m,v)*sum(fx);
  
  if (doplot)
    mmin = min([0.90*m, m_init]);
    mmax = max([1.10*m, m_init]);
    vmin = min([0.90*v, v_init]);
    vmax = max([1.10*v, v_init]);
    
    mmin = 0.95*mmin;
    mmax = 1.05*mmax;
    vmin = 0.95*vmin;
    vmax = 1.05*vmax;
    
    incr = 50;
    
    mm = linspace(mmin,mmax,incr);
    vv = linspace(vmin,vmax,incr);
    x = zeros(incr*incr,1);
    y = zeros(incr*incr,1);
    z = zeros(incr*incr,1);
    
    i = 0;
    for im = 1:incr
      for iv = 1:incr
        i = i+1;
        x(i) = mm(im);
        y(i) = vv(iv);
        z(i) = nbinodfitf([mm(im),vv(iv)],fx,stat);
      end;
    end;
    
    plottype = 2;
    ngrid = [];
    ptsymbol = 'ko';
    extrapts = [m,v,nbinodfitf([m,v],fx,stat); m_init,v_init,nbinodfitf([m_init,v_init],fx,stat)];
    contours = 1;
    color = 2;
    notsmoothed = 1;
    
    plotsurface(x,y,z,plottype,ngrid,ptsymbol,extrapts,contours,color,notsmoothed);
    hold on;
    plot3([m,m],[v,v],[0,extrapts(1,3)],'k');
    plot3([m_init,m_init],[v_init,v_init],[0,extrapts(2,3)],'k');
    hold off;
    box on;
    putxlab('Mean');
    putylab('Variance');
    switch (stat)
      case 0,
        putzlab('SSE');
      case 1,
        putzlab('Chi-squared deviation');
      case 2,
        putzlab('Sum of absolute deviations');
    end;
  end;
    
  return;
  
