% PearsonFitData:  Fit a continuous function to a histogram using Pearson distribution functions.
%           Returns the pdf as a set of 100 coordinates.    <NOTE: Make optional.>
%
%     Usage: crds = pearsonfitdata(midpoints,freqs,{relative},{noplot})
%                     OR
%            crds = pearsonfitdata(x,{nbins},{relative},{noplot})
%
%         midpoints = vector of midpoints of histogram bars.
%         freqs =     corresponding vector of frequencies.
%         x =         vector of data values
%         nbins =     optional integer indicating number of histogram bins
%                       [default = output of histbins()].
%         relative =  optional boolean flag indicating, if true, that histogram 
%                       and continuous function are to consist of relative 
%                       rather than absolute frequencies [default = 0 = absolute].
%         noplot =    optional boolean flag indicating, if true, that a plot of
%                       the histogram and function are not to be produced
%                       [default = 0 = produce plot].
%         ---------------------------------------------------------------------
%         crds =      [nmbrest x 2] matrix of point coordinates describing the 
%                       continuous function.
%

% Elderton, W.P. and N.L. Johnson. 1969. Systems of frequency curves.
%   Cambridge Univ Press.

% RE Strauss, 2/25/01
%   11/19/02 - Renamed from histfit().

function crds = pearsonfitdata(x,freqs,relative,noplot)
  if (nargin < 2) freqs = []; end;
  if (nargin < 3) relative = []; end;
  if (nargin < 4) noplot = []; end;

  if (isempty(relative))
    relative = 0;
  end;
  if (isempty(noplot))
    noplot = 0;
  end;

  if (isempty(freqs))                   % Bin the data
    [freqs,x] = histbins(x);
  elseif (isscalar(freqs))
    nbins = freqs;
    [freqs,x] = histbins(x,[],nbins);
  elseif (~isvector(x) | ~isvector(freqs))
    error('  PearsonFitData: midpoints and frequencies must be vectors.');
  elseif (length(x)~=length(freqs))
    error('  PearsonFitData: midpoint and frequency vectors not of same length.');
  end;
  x = x(:);
  freqs = freqs(:);
x_fx = [x freqs]

  delta = x(2)-x(1);                        % Extent of function
  xmin = x(1)-delta;
  xmax = x(end)+delta;
  x = [xmin; x; xmax];                      % Extend data
  freqs = [0; freqs; 0];

  N = sum(freqs);
  fx = freqs/N;

  m = meanwt(x,fx);
  nu(1) = 0;
  nu(2) = meanwt((x-m).^2,fx);
  nu(3) = meanwt((x-m).^3,fx);
  nu(4) = meanwt((x-m).^4,fx);
  
  mu(1) = m;
  mu(2) = nu(2) - 1/12;                     % Sheppard's adjustments
  mu(3) = nu(3);
  mu(4) = nu(4) - nu(2)/2 + 7/240;

  b1 = mu(3)^2 / mu(2)^3
  b2 = mu(4) / mu(2)^2
  kappa = (b1*(b2+3)^2)/(4*(4*b2-3*b1)*(2*b2-3*b1-6))   % Pearson criterion

  if (kappa<=-10)                     % Asymptotic
    disp('Not done yet');

  elseif (kappa>-10 & kappa<-0.2)         % Main
    [xx,yy] = pearson1(x,fx,mu,b1,b2);
    disp('Type 1');

  elseif (kappa>=-0.2 & kappa<=0.2)       % Transitional
    [xx,yy,mse1] = pearson1(x,fx,mu,b1,b2);
    [xx,yy,mse2] = pearson2(x,fx,mu,b1,b2,xmin,xmax);
    [xx,yy,mse3] = pearsonn(x,fx,mu,xmin,xmax);
    [xx,yy,mse4] = pearson7(x,fx,mu,b1,b2,xmin,xmax);
    [xx,yy,mse5] = pearson4(x,fx,mu,b1,b2,xmin,xmax);
mse = [mse1 mse2 mse3 mse4 mse5]

    [m,i] = min([mse1 mse2 mse3 mse4 mse5]);
i
    switch (i)
      case 1
        [xx,yy] = pearson1(x,fx,mu,b1,b2);
        disp('Type 1');
      case 2
        [xx,yy] = pearson2(x,fx,mu,b1,b2,xmin,xmax);
        disp('Type 2');
      case 3
        [xx,yy] = pearsonn(x,fx,mu,xmin,xmax);
        disp('Type normal');
      case 4
        [xx,yy] = pearson7(x,fx,mu,b1,b2,xmin,xmax);
        disp('Type 7');
      case 5
        [xx,yy] = pearson4(x,fx,mu,b1,b2,xmin,xmax);
        disp('Type 4');
    end;

  elseif (kappa>0.2 & kappa<0.8)          % Main
    [xx,yy] = pearson4(x,fx,mu,b1,b2,xmin,xmax);
    disp('Type 4');

  elseif (kappa>=0.8 & kappa<1.2)         % Transitional
    disp('Not done yet');

  elseif (kappa>=1.2 & kappa<10)          % Main
    [xx,yy] = pearson6(x,fx,mu,b1,b2,xmin,xmax);
    disp('Type 6');

  elseif (kappa>=10)                      % Asymptotic
    disp('Not done yet');

  end;

xx = real(xx);
yy = real(yy);
  crds = [xx yy];

  if (~noplot)
    figure;
    histgramb(x,fx,[],[],1,[],'w');
    hold on;
    plot(xx,yy,'k');
    hold off;
    putybnd(0,1.05*max([max(fx),max(yy)]));
    putylab('Relative frequency');
  end;

  return;
