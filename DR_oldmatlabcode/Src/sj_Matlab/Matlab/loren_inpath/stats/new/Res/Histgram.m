% HISTGRAM: Modification of hist(), producing a histogram with no space
%           between bars and with bars filled with a color.  Optionally
%           plots multiple histograms, with identical axes, for two or more 
%           groups.
%             Sturges' value for number of bins = 1+ceil((10/3)*log10(N))
%
%     Syntax: [nplot,freqs] = histgram(X,{grps},{labels},{bins}, ...
%                               {[histrows,histcols]},{noarrow},{scale}, ...
%                               {interval},{C},{sqplot},{fontsize})
%
%           X =        vector of values for which histogram is to be plotted.
%           grps =     optional grouping variable, producing vertically stacked
%                        histograms, one per group (passed as null if no groups).
%           labels =   matrix of text labels (by row) for titles of multiple 
%                        plots.
%           bins =    bins [default = Sturges' value for a single plot]: 
%                        - number of intervals if bins is scalar and data are 
%                             real, or
%                        - max number of intervals if bins is scalar and data
%                             are integer, or
%                        - interval midpoints if bins is a vector.
%           [histrows,histcols] = vector of dimensions of matrix of histogram 
%                        plots, for multiple groups [default = single column];
%                        if the vector = [1 1], a separate figure window is 
%                        opened for each plot but all are to a common scale.
%           noarrow =  if true (=1), suppresses printing of arrow indicating
%                        position of mean of data [default = 0].
%           scale =    frequency scale: 'abs' (absolute [default])
%                        or 'rel' (relative).
%           interval = 2-element vector [minval maxval] specifying a restricted
%                        interval of input data to plot.
%           C -        color [default 'k'] : a single character string chosen
%                        from the list 'r','g','b','c','m','y','w','k',
%                        or an RGB row-vector triple, [r g b]
%                        ([.5 .5 .5] for medium gray).
%           sqrplot =  if true (=1), prints square plots or subplots; 
%                        if false (=0), prints rectangular plots or subplots 
%                       [default = 0 for single plots, =1 for subplots].
%           fontsize = font size for axis tick labels [default = 10].
%           ----------------------------------------------------------------------
%           nplot =    number of points plotted.
%           freqs =    for single histogram, 2-column matrix of bin midpoints and 
%                        heights of bars.  For multiple-group histograms, 
%                        3-column matrix of histogram number, bin midpoints, and 
%                        heights of bars.
%

% RE Strauss, 7/2/96
%   9/ 9/96 - changed calling sequence;
%             ignore missing values;
%             allow specification of restricted interval of data.
%  11/ 2/96 - improved multiple-group histogram displays.
%   3/27/97 - added arrow indicating position of mean;
%             added subplot title labels.
%   6/ 8/97 - added capability to specify midpoints of bins;
%             corrected bug in calling sequence.
%   4/ 2/98 - allow for integer data for single histogram.
%   5/ 9/98 - return positions and heights of bars for single histogram.
%   7/10/98 - allows subplots to be rectangular rather than square.
%             plots x-ticklabels only on bottom row of subplots.
%   8/19/98 - allow variable fontsize for axis labels.
%   9/18/98 - allow for integer data for multiple histograms.
%   9/19/99 - changes for Matlab v5.
%  11/25/99 - simplified by partitioning into functions.
%   2/26/00 - correct problems with arrows; changed use of arrow().
%   9/12/00 - allow for single bar for invariant data.
%  10/17/00 - correct axis-scaling problem for multiple groups.
%   2/ 8/01 - use putylab() rather than ylabel() for 'Frequency' title;
%             change 'sqplot' argument to 'sqrplot'.
%   5/22/01 - add option of new figure window for each group histogram.
%   1/ 2/02 - changed arrow() to putarrow().
%   6/25/02 - changed fontsize of 'Frequency' label for multiple plots.
%   3/19/03 - for multiple plots, put group labels inside plot.
%   6/11/03 -  change headsize in putarrow() to default value.

function [nplot,freqs] = histgram(X,grps,labels,bins,histdimens,noarrow,...
                                  scale,interval,C,sqrplot,fontsize)
                                
  if (~nargin) help histgram; return; end;                                
  
  if (nargin < 2)  grps = []; end;
  if (nargin < 3)  labels = []; end;
  if (nargin < 4)  bins = []; end;
  if (nargin < 5)  histdimens = []; end;
  if (nargin < 6)  noarrow = []; end;
  if (nargin < 7)  scale = []; end;
  if (nargin < 8)  interval = []; end;
  if (nargin < 9)  C = []; end;
  if (nargin < 10) sqrplot = []; end;
  if (nargin < 11) fontsize = []; end;

  freqs = [];
  get_freqs = 0;
  if (nargout>1)
    get_freqs = 1;
  end;

  if (isempty(grps))
    ngrps = 1;
  else
    if (size(grps) ~= size(X))
      error('  HISTGRAM: grouping vector must be same size as data vector');
    else
      grp_ids = uniquef(grps,1);
      ngrps = length(grp_ids);
    end;
  end;

  X = X(:);                               % Convert matrices to vectors
  grps = grps(:);

  default_C = 'k';
  default_scale = 'abs';
  default_fontsize = 10;

  if (~isempty(labels))
    if (isvector(labels))                 % Convert vector of labels to col vector
      labels = labels(:);
    end;
    if (~ischar(labels))
      labels = tostr(labels);
    end;
    if (size(labels,1) ~= ngrps)
      error('  HISTGRAM: number of group labels not equal to number of groups');
    end;
  end;

  histrows = ngrps;
  histcols = 1;

  if (~isempty(histdimens))
    if (length(histdimens) == 2)
      histrows = histdimens(1);
      histcols = histdimens(2);
    else
      disp('  HISTGRAM Warning: invalid histogram-display dimension');
    end;
  end;

  if (~isempty(interval))
    if (length(interval) ~= 2)
      disp('  HISTGRAM Warning: invalid restricted interval');
      interval = [];
    end;
  end;

  if (isempty(scale))
    scale = default_scale;
  end;
  if (isempty(noarrow))
    noarrow = 0;
  end;
  if (isempty(C))
    C = default_C;
  end;
  if (isempty(sqrplot))
    if (isempty(grps))
      sqrplot = 0;
    else
      sqrplot = 1;
    end;
  end;
  if (isempty(fontsize))
    fontsize = default_fontsize;
  end;

  indx = find(isfinite(X));               % Eliminate missing values
  X = X(indx);
  if (~isempty(grps))
    grps = grps(indx);
  end;

  if (~isempty(interval))                 % Restrict input to specified interval
    indx = find(X>=min(interval) & X<=max(interval));
    X = X(indx);
    if (~isempty(grps))
      grps = grps(indx);
    end;
  end;

  nplot = length(X);


  % <<< Single plot >>>

  if (ngrps == 1)
    [n,centers,intrvl,v] = histbins(X,[],bins); % Put obs into bins

    if (~noarrow)                           % Prepare for arrow
      shaftlen = 0.12;                        % Relative length of shaft
      [position,shaftlen,v] = histarw(X,n,centers,v,shaftlen);
    end;

    if (scale == 'rel')                     % Relative frequency scale
      sumn = sum(n);
      n = n/sumn;
      v(4) = v(4)/sumn;
      if (~noarrow)
        position([2 4]) = position([2 4])/sumn;
      end;
    end;

    hold on;                                % Begin plot
    axis(v);                                % Set axes
    set(gca,'box','on');
    if (sqrplot)
      axis('square');
    end;

    histpltb(n,centers,C);                  % Draw histogram bars
    if (~noarrow)
      putarrow(position(3:4),position(1:2),sqrplot);   % Draw arrow
    end;
    plot(v(1:2),[0 0],'k');                 % Re-draw plot base

    set(gca,'FontSize',fontsize);           % Set font size of tick labels
    if (scale == 'abs')
      ytick = get(gca,'Ytick');               % Omit non-integer tick values
      if (~isintegr(ytick))
        ntick = 8;
        pertick = ceil(v(4)./ntick);
        puttick([],[0:pertick:floor(max(ytick))]);
      end;
    end;

    if (scale == 'rel')                     % Ordinate label
      putylab('Relative Frequency');
    else
      putylab('Frequency');
    end;

    hold off;

    if (get_freqs)
      freqs = [centers' n'];
    end;
  end;

  % <<< Plots by group >>>

  if (ngrps > 1)
    if (histrows>1 | histcols>1)
      figure;
    end;

    ymax = 0;
    shaftlen = 0;
    s = 0;

    for g = 1:ngrps                       % Cycle thru grps to get y-axis range
      Xg = X(grps==grp_ids(g));             % Isolate data for current group
      [n,centers,intrvl,v] = histbins(Xg,X,bins); % Put obs into bins

      if (~noarrow)                         % Get arrow shaft length based on all data
        s = 0.12;                           % Relative length of shaft
        [position,s,v] = histarw(Xg,n,centers,v,s);
      end;

      if (scale == 'rel')                   % Relative frequency scale
        sc = sum(n);
        v(4) = v(4)/sc;
      end;

      if (v(4)>ymax)                        % Save max values of ymax
        ymax = v(4);
      end;
      if (s>shaftlen)                       %   and arrow shaft-length
        shaftlen = s;
      end;
    end;
    v(4) = 1.1*ymax;
    if (shaftlen < 1)
      shaftlen = 1;
    end;

    for g = 1:ngrps                       % Cycle thru grps to plot histograms
      if (histrows==1 & histcols==1)
        figure;
      else
        subplot(histrows,histcols,g);
      end;
      hold on;
      Xg = X(grps==grp_ids(g));             % Isolate data for current group
      [n,centers,intrvl] = histbins(Xg,X,bins); % Put current obs into bins

      if (~noarrow)                         % Prepare for arrow
        position = histarw(Xg,n,centers,v,shaftlen);
      end;

      if (scale == 'rel')                   % Relative frequency scale
        sc = sum(n);
        n = n/sc;
        if (~noarrow)
          position([2 4]) = position([2 4])/sc;
        end;
      end;

      hold on;                                % Begin plot
      set(gca,'box','on');
      axis(v);                                % Set axes
      if (sqrplot)
        axis('square');
      end;

      histpltb(n,centers,C);                  % Draw histogram bars

      if (~noarrow)
        putarrow(position(3:4),position(1:2),sqrplot);   % Draw arrow
      end;
      plot(v(1:2),[0 0],'k');                 % Re-draw plot base

      set(gca,'FontSize',fontsize);           % Set font size of tick labels
      if (scale == 'abs')
        ytick = get(gca,'Ytick');               % Omit non-integer tick values
        if (~isintegr(ytick))
          ntick = 8;
          pertick = ceil(v(4)./ntick);
          puttick([],[0:pertick:floor(max(ytick))]);
        end;
      end;
      putylab('Frequency',[],22/max([histrows,histcols]));

      if (histrows>1 | histcols>1)
        if (g <= (histrows-1)*histcols)         % Plot x ticklabels only on bottom row of plots
          set(gca,'XtickLabel',[]);
        end;
        set(gca,'FontSize',22/max([histrows,histcols]));
      end;

      if (~isempty(labels))
        if (histrows>1 | histcols>1)
%           text(v(1),v(4)+0.09*(v(4)-v(3)),labels(g,:));
          puttext(0.05,0.90,labels(g,:),[],13);
        else
          puttitle(labels(g,:));
        end;
      end;

      hold off;

      if (get_freqs)
        freqs = [freqs; g*ones(length(n),1) centers' n'];
      end;
    end;
  end;

  return;

