% TRIALLOM: Given log-transformed data on three sides of a triangular structure, 
%           a log-space size variable, and a specified set of standard sizes, 
%           predicts and triangle sizes at the standard sizes (by separate linear 
%           regressions of the three side variables on size) and reconstructs 
%           the triangles as sets of coordinates.
%             The first side variable is taken for the base, the second variable 
%           the left size, and the third variable the right side of the triangle.
%
%     Usage: [crds,b,stats] = triallom(sides,tsize,stdsize)
%
%           sides =   [n x 3] matrix of values of n observations of a triangular 
%                       structure.
%           tsize =   size vector (length n).
%           stdsize = vector (length k) of standard sizes.
%           ----------------------------------------------------------------------
%           crds =    [k x 6] matrix of coordinates of triangle apices, in which 
%                       each row consists of [0 0, x2 0, x3 y3]. 
%           b =       [2 x 3] matrix of regression coefficients; the first row
%                       gives the intercepts, the second row gives the slopes.
%           stats =   [6 x 3] matrix of regression statistics; each column gives the
%                     intercepts for a single 'side' variable:
%                       row 1: adjusted coefficient of determination (R_a^2)
%                           2: MSE (residual variances)
%                           3: F-statistic value
%                           4: df1 (numerator degrees of freedom, =dfr)
%                           5: df2 (denominator degrees of freedom = dfe)
%                           6: pr >= F (probability under the null hypothesis)
%

% RE Strauss, 7/18/99

function [crds,b,stats] = triallom(sides,tsize,stdsize)
  [n,p] = size(sides);
  if (p~=3)
    error('  TRIALLOM: sides matrix not of size [n x 3]');
  end;

  [r,c] = size(tsize);
  if (r~=n)
    if (r==1 & c==n)
      tsize = tsize';
    else
      error('TRIALLOM: sides and tsize matrices not compatible');
    end;
  end;

  [r,k] = size(stdsize);
  if (min([r,c])~=1)
    error('TRIALLOM: stdsize must be a vector of standard sizes');
  end;
  if (r==1)
    stdsize = stdsize';
  else
    k = r;
  end; 

  [b,stats,pred_sides] = linregr(tsize,sides,[],stdsize);
  pred_sides = exp(pred_sides);       % Inverse transformation of logarithmic data

  crds = zeros(k,6);
  for i = 1:k
    crds(i,3) = pred_sides(i,1);      % Base
    c = [0 0; crds(i,3:4)];
    d = pred_sides(i,2:3);
    P = triangpt(c,d);            % Triangulate third point
    crds(i,5:6) = [max(P)];
  end;

  return;
