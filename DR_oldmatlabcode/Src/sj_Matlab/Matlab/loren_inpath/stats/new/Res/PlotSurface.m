% PlotSurface: Given a set of 3D point coordinates, plots a mesh or surface plot.
%              A colorized, smoothed surface plot is produced by default.
%
%     Usage: [X,Y,Z] = plotsurface(x,y,z,{plottype},{ngrid},{ptsymbol},{extrapts},...
%                        {contours},{color},{notsmoothed})
%
%         x,y,z =       corresponding column vectors of length n.
%         plottype =    optional flag indicating the type of plot to be produced:
%                         0: no plot;
%                         1: surface plot [default];
%                         2: mesh plot;
%                         3: waterfall plot.
%         ngrid =       optional 2-element vector of numbers of grid positions
%                         in the x- and y-directions [default = [40,40]]; if a single
%                         value is passed, it is used for both dimensions.
%         ptsymbol =    optional plot symbol for points to be superimposed onto 
%                         the surface plot (e.g., 'ko'), either the original points 
%                         (if 'extrapts' are not supplied) or extra points (if they
%                         are).
%         extrapts =    optional 3-column matrix of (x,y,z) coordinates to be superimposed
%                         onto the plot.
%         contours =    optional boolean flag indicating that a contour plot is to 
%                         be superimposed onto the floor of the surface plot
%                         [default = 0].
%         color =       optional flag indicating the color scheme for the plot:
%                         0: default colormap [default];
%                         1: gray-scale;
%                         2: black (for mesh or waterfall plots).
%         notsmoothed = optional boolean flag indicating that the surface is not to 
%                         be smoothed (via cubic interpolation) before being plotted
%                         [default = 0].
%         ---------------------------------------------------------------------------------
%         X,Y,Z =       corresponding square matrices in standard Matlab 3D format.
%

% RE Strauss, 11/1/02
%   11/12/02 - added ability to superimpose points onto surface plot;
%              rearranged input arguments; other miscellaneous changes.
%   2/25/03 -  allow for superimposition of additional points.
%   4/19/03 -  output coordinate matrices in standard Matlab 3D format.

function [X,Y,Z] = plotsurface(x,y,z,plottype,ngrid,ptsymbol,extrapts,...
                               contours,color,notsmoothed)
                             
  if (nargin <  4) plottype = []; end;
  if (nargin <  5) ngrid = []; end;
  if (nargin <  6) ptsymbol = []; end;
  if (nargin <  7) extrapts = []; end;
  if (nargin <  8) contours = []; end;
  if (nargin <  9) color = []; end;
  if (nargin < 10) notsmoothed = []; end;
  
  plotpoints = 0;
  if (~isempty(ptsymbol))
    plotpoints = 1;
  end;
  if (isempty(ngrid))
    ngrid = [40,40];
  end;
  if (length(ngrid)==1)
    ngrid = [ngrid,ngrid];
  end;
  
  if (isempty(plottype))
    plottype = 0;
  end;
  if (isempty(notsmoothed))
    smoothed = 0;
  end;
  if (isempty(contours))
    contours = 0;
  end;
  if (isempty(color))
    color = 0;
  end;
  
  switch (color)                          % Set colormap
    case 0,
    case 1,
      colormap('gray');
    case 2,
      colormap(zeros(size(colormap)));
    otherwise
      error('  PlotSurface: invalid value for color');  
  end;
  
  xi = linspace(min(x),max(x),ngrid(1));  % Create grid via Delaunay triangulation
  yi = linspace(min(y),max(y),ngrid(2));
  [X,Y] = meshgrid(xi,yi);   
  if (notsmoothed)
    Z = griddata(x,y,z,X,Y,'linear');
  else
    Z = griddata(x,y,z,X,Y,'cubic');
  end;
  
  switch (plottype)
    case 0,                               % No plot
    case 1,                               % Surface plot
      if (contours)                       
        surfc(X,Y,Z);
      else
        surf(X,Y,Z);
      end;
    case 1,                               % Mesh plot
      if (contours)
        meshc(X,Y,Z);
      else
        mesh(X,Y,Z);
      end;
    case 2,                               % Waterfall plot
      waterfall(X,Y,Z);
    otherwise
      error('  PlotSurface: invalid plot type');
  end;
  
  if (plotpoints)
    hold on;
    if (isempty(extrapts))
      plot3(x,y,z,ptsymbol);
    else
      plot3(extrapts(:,1),extrapts(:,2),extrapts(:,3),ptsymbol);
    end;
    hold off;
  end;
  
  return;
  