% PrctileRange:  Finds the range containing the central P% of the data for a single 
%               variable, where P is the percentage of central points included.  If
%               P = 50, the result is a median range.  Proceeds by finding the range 
%               of increasing numbers of central points and predicting the range of the 
%               central half of the points.  Optionally produces a plot of the functions 
%               of number of central points vs range.
%
%     Usage: pctrange = prctilerange(x,{prctile},{doplot}}
%
%         x =       vector of data values.
%         prctile = optional percentile value [default = 50].
%         doplot =  optional boolean variable indicating, if true, that a plot is to be
%                     produced.
%         -----------------------------------------------------------------------------
%         medrange = predicted median range.
%

% RE Strauss, 3/25/02

function pctrange = prctilerange(x,prctile,doplot)
  if (prctile<2) prctile = []; end;
  if (nargin<3) doplot = []; end;
  
  if (isempty(prctile))
    prctile = 50;
  end;
  if (isempty(doplot))
    doplot = 0;
  end;
  
  if (prctile>1)
    prctile = prctile/100;
  end;
  if (prctile<0 | prctile>1)
    error('  PrctileRange: percentile value out of range.');
  end;
  
  x = sort(x(:));
  npts = length(x);
  
  rng1 = zeros(size(x));
  rng2 = zeros(size(x));
  
  seq1 = [1:npts; npts:-1:1];
  seq1 = seq1(:);
  seq1 = flipud(seq1(1:npts));

  seq2 = [npts:-1:1; 1:npts];
  seq2 = seq2(:);
  seq2 = flipud(seq2(1:npts));
  
  for i = 1:npts
    rng1(i) = range(x(seq1(1:i)));
    rng2(i) = range(x(seq2(1:i)));
  end;
  
  pctpos = prctile*(npts+1);
  pctrange = mean([rng1(floor(pctpos)) rng1(ceil(pctpos)) rng2(floor(pctpos)) rng2(ceil(pctpos))]);
  
  if (doplot)
    plot(1:npts,rng1,'b',1:npts,rng2,'r');
    putxlab('Number of central points');
    putylab('Range');
  end;

  return;
  