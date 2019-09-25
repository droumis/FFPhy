% HISTFIT:  Fit a continuous function to a histogram using the first four 
%           moments about the mean:
%             y = b0 + b1(d) + b2(d^2) + b3(d^3) + b4(d^4)
%           where d = x-xbar, the deviations from the mean.
%
%     Usage: [crds,b] = histfit(midpoints,freqs,{nmbrest},{noplot})
%
%         midpoints = vector of midpoints of histogram bars.
%         freqs =     corresponding vector of frequencies.
%         nmbrest =   optional max number of predicted points (x,fx) to be 
%                       returned (fewer might be returned) [default = 100].
%         noplot =    optional boolean flag indicating, if true, that a plot of
%                       the histogram and function are not to be produced
%                       [default = 0 = produce plot].
%         ---------------------------------------------------------------------
%         crds =      [nmbrest x 2] matrix of point coordinates describing the 
%                       continuous function, for plotting.
%         b =         column vector of regression coefficients.
%

% RE Strauss, 2/25/01

function [crds,b] = histfit(midpoints,freqs,nmbrest,noplot)
  if (nargin < 3) nmbrest = []; end;
  if (nargin < 4) noplot = []; end;

  if (isempty(nmbrest))
    nmbrest = 100;
  end;
  if (isempty(noplot))
    noplot = 0;
  end;

  if (~isvector(midpoints) | ~isvector(freqs))
    error('  HISTFIT: midpoints and frequencies must be vectors.');
  end;

  if (length(midpoints)~=length(freqs))
    error('  HISTFIT: midpoint and frequency vectors not of same length.');
  end;

  midpoints = midpoints(:);
  freqs = freqs(:);

  delta = midpoints(2)-midpoints(1);    % Append zeros on either side
  midpoints = [midpoints(1)-delta; midpoints; midpoints(end)+delta];
  freqs = [0; freqs; 0];

  save_midpoints = midpoints;
  save_freqs = freqs;
  midpoints = [midpoints*ones(1,3)]';
  midpoints = midpoints(:);
  freqs = [freqs*ones(1,3)]';
  freqs = freqs(:);

  m = mean(midpoints);
  d1 = midpoints - m;
  d2 = d1.^2;
  d3 = d1.^3;
  d4 = d1.^4;

  x = linspace(min(midpoints),max(midpoints),nmbrest)';
  xd1 = x - m;
  xd2 = xd1.^2;
  xd3 = xd1.^3;
  xd4 = xd1.^4;
  
  b0 = linregr([d1 d2 d3 d4],freqs);                          % Initial fitting
  options = optimset('Display','off');
  b = fminsearch('histfitf',b0,options,[d1 d2 d3 d4],freqs);  % Constrain to nonneg 
                                                              %   predicted values
  fx = [ones(nmbrest,1) xd1 xd2 xd3 xd4]*b;                   % Predicted

  q = ceil(nmbrest/4);
  i = find(fx(1:q-1)-fx(2:q) > 0);
  x(i) = [];
  fx(i) = [];
  fx = flipud(fx);
  i = find(fx(1:q-1)-fx(2:q) > 0);
  x(i) = [];
  fx(i) = [];
  fx = flipud(fx);

  crds = [x fx];

  if (~noplot)
    figure;
    histgramb(save_midpoints,save_freqs,[],[],1,[],'w');
    hold on;
    plot(x,fx,'k');
    hold off;
  end;

  return;
