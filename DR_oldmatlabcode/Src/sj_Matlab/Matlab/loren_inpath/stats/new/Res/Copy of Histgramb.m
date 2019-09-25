% HISTGRAMB: Produces histogram, given centers and heights of bars.
%
%     Usage: xmean = histgramb(centers,heights,{usecent},{uselines},{noarrow},...
%                              {sqrplot},{C},{fontsize})
%
%         centers = sorted vector of midpoints of bins.
%         heights = corresponding vector of bar heights.
%         usecent = optional boolean vector indicating, if true, that the 
%                     specified bar midpoints are to be used as axis tickmarks 
%                     [default = standard axis labeling].
%         uselines = optional boolean vector indicating, if true, that vertical 
%                     lines rather than bars are to be used to portray the 
%                     histogram; useful for large numbers of bars.
%         noarrow = optional boolean variable indicating, if true, that the 
%                     arrow indicating the position of the mean is to be 
%                     suppressed [default = 0].
%         sqrplot = optional boolean variable indicating, if true, that the plot 
%                     is to be square; if false, prints rectangular plot 
%                     [default = 0].
%         C =       optional color [default 'k'] : a single character string 
%                     chosen from the list 'r','g','b','c','m','y','w','k',
%                     or an RGB row-vector triple, [r g b]
%                     ([.5 .5 .5] for medium gray).
%         fontsize = optional font size for axis tick labels [default = 10].
%         ----------------------------------------------------------------------
%         xmean =   weighted mean of data represented by bar centers.
%

% RE Strauss, 2/8/01
%   2/19/00 - allow for vertical lines rather than bars.
%   1/2/02 -  changed arrow() to putarrow().

function xmean = histgramb(centers,heights,usecent,uselines,noarrow,sqrplot,C,fontsize)
  if (nargin < 3) usecent = []; end;
  if (nargin < 4) uselines = []; end;
  if (nargin < 5) noarrow = []; end;
  if (nargin < 6) sqrplot = []; end;
  if (nargin < 7) C = []; end;
  if (nargin < 8) fontsize = []; end;

  if (isempty(noarrow))
    noarrow = 0;
  end;
  if (isempty(C))
    C = 'k';
  end;
  if (isempty(sqrplot))
    sqrplot = 0;
  end;
  if (isempty(fontsize))
    fontsize = 10;
  end;
  if (isempty(usecent))
    usecent = 0;
  end;
  if (isempty(uselines))
    uselines = 0;
  end;

  nbar = length(centers);
  if (length(heights)~=nbar)
    error('  HISTGRAMB: centers and heights vectors not compatible.');
  end;
  centers = centers(:);
  heights = heights(:);

  intrvl = mean(centers(2:nbar)-centers(1:nbar-1)); % Interval between midpoints

  if (uselines)
    hold on;
    v = putbnd(centers,heights);
    for i = 1:length(centers)
      plot([centers([i,i])],[0 heights(i)],C);
    end;
    hold off;
  else
    histpltb(heights,centers,C);            % Produce initial plot
    v = axis;
  end;

  v(1) = centers(1)-intrvl;
  v(2) = centers(nbar)+intrvl;

  xmean = meanwt(centers,heights);        % Weighted mean
  if (noarrow)
    v(4) = 1.05 * max(heights);
  else                                    % Prepare for arrow
    shaftlen = 0.12;                        % Relative length of shaft
    [position,shaftlen,v] = histarw(xmean,heights,centers,v,shaftlen);
  end;
  axis(v);                                % Adjust axis ranges

  if (sqrplot)
    axis('square');
  end;
  if (~noarrow)                           % Draw arrow
    putarrow(position(3:4),position(1:2),sqrplot,0.035);   
  end;
  if (~uselines)
    plot(v(1:2),[0 0],'k');                 % Re-draw plot base
  end;
  set(gca,'FontSize',fontsize);           % Set font size of tick labels

  if (usecent)                            % Use midpoints as x tickmarks
    puttick(centers);                       
  end;

  putylab('Frequency');

  return;

