% Vectplot: Plot vector diagram of loadings (correlations), centered on origin
%
%     Syntax:  vectplot(crds,{labl},{ci},{floatbnds},{nolabl},{labldist},...
%                         {thresh},{prtarrow},{fontsize})
%
%         crds -      [n x 1] or [n x 2] matrix coordinates (loadings).
%         labl -      optional matrix of labels (by row) corresponding to vectors
%                       [default: integer labels].
%         ci -        optional [n x 2] matrix of lower & upper confidence limits for a 
%                       single column of loadings.  Ignored for 2D vector plots.
%         floatbnds - optional boolean flag indicating, if true, that axis 
%                       bounds are not to be fixed at {-1, +1} [default = 0].
%         nolabl -    optional boolean flag indicating, if true, that labels are not 
%                       to be printed [default = 0].
%         labldist -  optional multiplicative scale factor to increase or decrease 
%                       distances of labels from vectors [default = 1].
%         thresh -    optional threshhold vector length for plotting labels [default = 0].
%         prtarrow -  optional boolean flag indicating, if true, that vectors are 
%                       to be plotted with arrowheads [default = 1].
%         fontsize -  optional font size for labels [default = 10].
%

% RE Strauss, 6/95
%   7/17/98 - added options: nolabl, labldist, thresh, prtarrow.
%   9/10/98 - added plot for 1-dimensional loadings, fontsize.
%   9/16/98 - added ability to plot confidence limits for 1-dimensional loadings.
%   7/13/99 - added ability to plot 1-dimensional allometric coefficients.
%   9/3/99 -  changes to plots for Matlab v5.
%   1/2/00 -  allow for confidence limits out of sequence (high, low).
%   2/26/00 - allow for axis bounds beyond 1; changed use of arrow().
%   5/4/00 -  allow for crds matrix to have >2 cols.
%   5/27/00 - print all labels to the right for 1-dimensional plots.
%   6/27/00 - added floatbnds option; isolated vectplot1 and vectplot2.
%   7/1/02 -  added size-check for crds and labl matrices.
%   9/9/02 -  corrected problem with 7/1/02 correction.

function vectplot(crds,labl,ci,floatbnds,nolabl,labldist,thresh,prtarrow,fontsize)
  if (nargin < 2) labl = []; end;
  if (nargin < 3) ci = []; end;
  if (nargin < 4) floatbnds = []; end;
  if (nargin < 5) nolabl = []; end;
  if (nargin < 6) labldist = []; end;
  if (nargin < 7) thresh = []; end;
  if (nargin < 8) prtarrow = []; end;
  if (nargin < 9) fontsize = []; end;

  [nvars,naxes] = size(crds);
  if (naxes>2)
    crds = crds(:,1:2);
    [nvars,naxes] = size(crds);
  end;

  plotci = 0;
  if (~isempty(ci) & naxes==1)
    [nci,nlim] = size(ci);
    if (nlim~=2)
      error('  VECTPLOT: must have 2 cols of confidence limits for loadings.');
    end;
    if (nci~=nvars)
      disp('  VECTPLOT: loadings and confidence-interval matrices');
      error('            must have same number of rows.');
    end;
    plotci = 1;
  end;
  
  if (isvector(labl))
    labl = labl(:);
  end;
  if (~isempty(labl) & size(labl,1)~=size(crds,1))
    error('  VECTPLOT: loadings and labels must have same number of rows.');
  end;

  if (isempty(floatbnds))
    floatbnds = 0;
  end;
  if (isempty(nolabl))
    nolabl = 0;
  end;
  if (isempty(labldist))
    labldist = 1;
  end;
  if (isempty(thresh))
    thresh = 0;
  end;
  if (isempty(prtarrow))
    prtarrow = 1;
  end;
  if (isempty(fontsize))
    fontsize = 10;
  end;

  if (~isempty(labl) & ~ischar(labl))             % Convert numeric labels to char
    labl = tostr(labl);
  end;
  if (isempty(labl))                              % Default labels
    labl = tostr(1:nvars);
  end;

  switch (naxes)
    case 1,
      vectplot1(crds,labl,ci,plotci,floatbnds,nolabl,labldist,thresh,prtarrow,fontsize);

    case 2,
      vectplot2(crds,labl,floatbnds,nolabl,labldist,thresh,prtarrow,fontsize);

    otherwise
      error('  VECTPLOT: loadings matrix must have 1 or 2 columns (axes)');
  end;

  return;
