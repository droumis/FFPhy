% Vectplot2: Plot 2D vector diagram of loadings (correlations), centered on origin
%
%     Syntax:  vectplot2(crds,labl,floatbnds,nolabl,labldist,...
%                         thresh,prtarrow,fontsize)
%
%         crds -      [n x 1] or [n x 2] matrix coordinates (loadings).
%         labl -      matrix of labels (by row) corresponding to vectors.
%         floatbnds - boolean flag indicating, if true, that axis 
%                       bounds are not to be fixed at {-1, +1}.
%         nolabl -    boolean flag indicating, if true, that labels are not 
%                       to be printed.
%         labldist -  multiplicative scale factor to increase or decrease 
%                       distances of labels from vectors.
%         thresh -    threshhold vector length for plotting labels.
%         prtarrow -  boolean flag indicating, if true, that vectors are 
%                       to be plotted with arrowheads.
%         fontsize -  font size for labels.
%

% RE Strauss, 6/27/00 (isolated from vectplot)
%   1/2/02 - changed arrow() to putarrow().
%   4/12/02 - corrected syntax problem with call to putarrow().
%   6/11/03 -  change headsize in putarrow() to default value.

function vectplot2(crds,labl,floatbnds,nolabl,labldist,thresh,prtarrow,fontsize)
  if (floatbnds)
    bnd = 1.1*max([abs(min(crds)), abs(max(crds))]); % Square plot bounds
  else
    bnd = 1.1*max([abs(min(crds)), abs(max(crds)), 1]); % Square plot bounds
  end;

  plot([0 0],[0 0]);
  axis([-bnd bnd -bnd bnd]);
  axis('square');
  hold on;

  for i = 1:max(size(crds))
    lenvect = sqrt(sum(crds(i,:).^2));
    if (prtarrow)
      putarrow([0 0],crds(i,:),1);
    else
      putarrow([0 0],crds(i,:),1,0);
    end;

    if (crds(i,1)>=0 & crds(i,2)>=0)
      deltax = 0.03;
      deltay = 0.03;
    elseif (crds(i,1)<0 & crds(i,2)>=0)     
      deltax = -0.07;
      deltay = 0.03;
    elseif (crds(i,1)<0 & crds(i,2)<0)   
      deltax = -0.07;
      deltay = -0.04;
    elseif (crds(i,1)>=0 & crds(i,2)<0)
      deltax = 0.03;
      deltay = -0.04;
    end;
    deltax = deltax*bnd;
    deltay = deltay*bnd;

    if (lenvect > thresh)
      if (~nolabl)
        h = text(crds(i,1)+labldist*deltax, crds(i,2)+labldist*deltay, labl(i,:));
        set(h,'fontsize',fontsize);
      end;
    end;
  end;

  hold off;

  return;