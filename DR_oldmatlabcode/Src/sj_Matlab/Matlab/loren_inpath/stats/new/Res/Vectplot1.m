% Vectplot1: Plot 1D vector diagram of loadings (correlations), centered on origin
%
%     Syntax:  vectplot1(crds,{labl},{ci},{floatbnds},{nolabl},{labldist},...
%                         {thresh},{prtarrow},{fontsize})
%
%         crds -      [n x 1] or [n x 2] matrix coordinates (loadings).
%         labl -      optional matrix of labels (by row) corresponding to vectors
%                       [default: integer labels].
%         ci -        optional [n x 2] matrix of lower & upper confidence limits 
%                       sifor a ngle column of loadings.
%         plotci -    boolean flag indicating whether confidence intervals are 
%                       to be plotted.
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
%   6/12/02 - corrected syntax error for single-variable plot.
%   6/11/03 -  change headsize in putarrow() to default value.

function vectplot1(crds,labl,ci,plotci,floatbnds,nolabl,labldist,thresh,prtarrow,fontsize)
  [nvars,naxes] = size(crds);

  if (~plotci)                                   % Find min & max values to be plotted
    vals = [crds(:)];
  else
    vals = [crds(:); ci(:)];
  end;
  minval = min(vals);
  maxval = max(vals);

  alloms = 0;                                   % Check if values are allometries
  if (minval<1 & maxval>1)
    alloms = 1;
  end;

  if (alloms | floatbnds)
    rval = maxval - minval;
    v = [minval-0.05*rval maxval+0.05*rval -(nvars+1) 0];
  else
    v = [-1.1 1.1 -(nvars+1) 0];
  end;
  axis_range = v(2)-v(1);

  plot([0 0],[0 0]);
  axis(v);
  axis('square');
  hold on;
    
  if (alloms)
    plot([1 1],[-(nvars+1) 0],'k:');
    origin = 1;
  else
    plot([0 0],[-(nvars+1) 0],'k:');
    origin = 0;
  end;

  if (~plotci)
    ci = [crds crds];
  end;

  for i = 1:nvars
    lenvect = abs(crds(i)-origin);               

    ci1 = min(ci(i,:));                         % Confidence limits
    ci2 = max(ci(i,:));
    ci1 = ci1(1);
    ci2 = ci2(1);

   if (~plotci)                                 % Plot arrows or lines
      if (prtarrow & lenvect>0.06)
        putarrow([origin -i],[crds(i) -i],1);
      else
        putarrow([origin -i],[crds(i) -i],1,0);
      end;

    else                                        % Plot confidence limits
      deltax = 0.01*axis_range;
      deltay = 0.2;

      plot([ci1 ci2],[-i -i],'k');                          % Horizontal line
      plot([crds(i) crds(i)],[-i-deltay -i+deltay],'k');    % Tick mark

      plot([ci1 ci1],[-i-deltay -i+deltay],'k');            % Lower conf limit
      plot([ci1 ci1+deltax],[-i-deltay -i-deltay],'k');
      plot([ci1 ci1+deltax],[-i+deltay -i+deltay],'k');
      plot([ci2 ci2],[-i-deltay -i+deltay],'k');            % Upper conf limit
      plot([ci2 ci2-deltax],[-i-deltay -i-deltay],'k');
      plot([ci2 ci2-deltax],[-i+deltay -i+deltay],'k');
    end;

    deltax = 0.03*axis_range;
    deltay = 0.03;

    if (lenvect > thresh)                       % Plot labels
      if (~nolabl)
        if (crds(i) > origin)
          h = text(ci2+labldist*deltax, deltay-i, labl(i,:));
          set(h,'fontsize',fontsize);
        else
          h = text(labldist*deltax, deltay-i, labl(i,:));
          set(h,'fontsize',fontsize);
        end;
      end;
    end;
  end;

  puttick([],'off');
  hold off;

  return;
