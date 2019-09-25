% HISTBINS: Puts data vector into bins for histogram.
%
%     Usage: [n,centers,intrvl,v] = histbins(x,{xfull},{bins},{ignint})
%
%           x =       [n x 1] vector of data values.
%           xfull =   full vector, if x represents one group or variable of a 
%                       set [default xfull = x].
%           bins =    number of bins, if scalar [default = Sturges' value];
%                       bin centers, if vector.
%           ignint =  optional boolean flag indicating, if true, that a vector 
%                       of integers is to be treated as vector of real numbers 
%                       [default = 0].
%           -------------------------------------------------------------------
%           n =       vector of counts per bin.
%           centers = corresponding centers of bins.
%           intrvl =  distance between bin centers.
%           v =       axis ranges: [xmin xmax ymin ymax].
%

% RE Strauss, 11/25/99
%   9/12/00 - allow for single bin.
%  10/17/00 - correct problems with histograms for multiple groups.

function [n,centers,intrvl,v] = histbins(x,xfull,bins,ignint)
  if (nargin < 2) xfull = []; end;
  if (nargin < 3) bins = []; end;
  if (nargin < 4) ignint = []; end;

  if (isempty(xfull))
    xfull = x;
    x = [];
  end;

  if (isempty(ignint))
    ignint = 0;
  end;

  if (isempty(bins))                      % Default bins
    bins = 1 + ceil((10/3)*log10(length(xfull)));  % Sturges' value
  end;
  if (length(bins)<2 & isintegr(xfull) & ~ignint)  % Bins for integers
    perbin = floor(range(xfull)./bins)+1;
    bins = [min(xfull):perbin:max(xfull)];
  end;

  [n,centers] = hist(xfull,bins);         % Bin counts for complete data
  if (~isempty(x))
    [n,centers] = hist(x,centers);        % Bin counts for current data
  end;

  if (length(centers)==1)                 % Single bin
    n = [0 n 0];
    centers = [centers-1 centers centers+1];
  end;

  intrvl = centers(2)-centers(1);         % Interval between bin centers
  xmin = centers(1)-intrvl;
  xmax = centers(end)+intrvl;
  ymin = 0;
  ymax = 1.05 * max(n);

  v = [xmin xmax ymin ymax];

  return;


