% GETAVGCOLOR: Given a color image (RGB or HSV) and a polygon, returns the median color 
%           parameters for the pixels within the polygon.
%
%     Usage: [median_values,univar_stats,cov_stats] = ...
%                                             getavgcolor(I,{crds},{convert},{doplots})
%
%           I =       3D matrix [M x N x 3] of true-color (RGB) or HSV image.
%           crds =    optional [npts x 2] matrix of polygon coordinates. If provided,
%                       finds the subset of pixels within the polygon specified
%                       by the coordinates.
%           convert = optional flag indicating that:
%                       0 = no conversion of values [default];
%                       1 = RGB should be converted to HSV.
%           doplots = optional boolean flag indicating, if true, that histograms
%                       of the distributions of the 3 color parameters are to be
%                       plotted [default = 0].
%           ----------------------------------------------------------------------
%           median_values = median color-parameter values within the polygon.
%           univar_stats =  [1 x 24] vector of univariate statistics from univar(),
%                             8 statistics for each of the 3 color parameters.
%           cov_stats =     [3 x 3 matrix of covariances and correlations among the
%                             3 color parameters:  variances on the diagonal,
%                             covariances in the upper triangular matrix, and
%                             correlations in the lower triangular matrix.
%

% RE Strauss, 06/26/02

function [median_values,univar_stats,cov_stats] = getavgcolor(I,crds,convert,doplots)
  if (nargin < 2) crds = []; end;
  if (nargin < 3) npixels = []; end;
  if (nargin < 4) doplots = []; end;
  
  do_univar = 0;
  do_cov = 0;
  if (nargout > 1)
    do_univar = 1;
  end;
  if (nargout > 2)
    do_cov = 1;
  end;
  
  if (isempty(doplots))
    doplots = 0;
  end;
  
  if (convert)
    if (isrgb(I))
      I = rgb2hsv(I);
    else
      I = hsv2rgb(I);
    end;
  end;

  [r,c,d] = size(I);
  
  Ipts = zeros(r*c,2);
  Icolor = zeros(r*c,3);

  i = 0;
  for ir = 1:r
    for ic = 1:c
      i = i+1;
      Ipts(i,:) = [ir,ic];
      Icolor(i,:) = double(I(ir,ic,:));
    end;
  end;
  
  if (~isempty(crds))                           % Find subset of pixels within polygon
    isin = isinpoly(Ipts,crds);
    ind = find(isin>0);  
  else
    ind = [1:size(Icolor,1)]';
  end;
 
  if (do_univar)
    univar_stats = zeros(1,24);
    median_values = zeros(1,3);
    
    univar_stats(1:8) = univar(Icolor(ind,1));
    univar_stats(9:16) = univar(Icolor(ind,2));
    univar_stats(17:24) = univar(Icolor(ind,3));

    median_values(1) = univar_stats(2);
    median_values(2) = univar_stats(10);
    median_values(3) = univar_stats(18);
  else
    median_values = median(Icolor(ind,:));
  end;
  
  if (do_cov)
    cov_stats = cov(Icolor(ind,:));  
    if (all(diag(cov_stats)>eps))  
      cov_stats(2,1) = cov_stats(1,2)/(sqrt(cov_stats(1,1))*sqrt(cov_stats(2,2)));
      cov_stats(3,1) = cov_stats(1,3)/(sqrt(cov_stats(1,1))*sqrt(cov_stats(3,3)));
      cov_stats(3,2) = cov_stats(2,3)/(sqrt(cov_stats(2,2))*sqrt(cov_stats(3,3)));
    else
      cov_stats(2,1) = NaN;
      cov_stats(3,1) = NaN;
      cov_stats(3,2) = NaN;
    end;
  end;
  
  if (doplots)
    figure;
    histgram(Icolor(ind,1));
    puttitle('Color parameter 1');
    figure;
    histgram(Icolor(ind,2));
    puttitle('Color parameter 2');
    figure;
    histgram(Icolor(ind,3));
    puttitle('Color parameter 3');
  end;
  
  return;
