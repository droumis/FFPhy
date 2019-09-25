function hh = errorbar2(x, y, e, w, symbol, varargin)
%ERRORBAR2 Error bar plot.
%   ERRORBAR2(X,Y,E,W) plots the graph of vector X vs. vector Y with
%   error bars specified by the vector E, with a top line total 
%   length of W. 
%
%   If X and Y are vectors and E is a two row matrix, the error bar bottom and
%   top magnitudes are taken from E(:,1) and E(:,2) respectively
%   If X,Y, and E are matrices then each column produces a separate line.
%
%   ERRORBAR2(...,'LineSpec') uses the color and linestyle specified by
%   the string 'LineSpec'.  See PLOT for possibilities.
%
%   H = ERRORBAR2(...) returns a vector of line handles.
%
%   For example,
%      x = 1:10;
%      y = sin(x);
%      e = std(y)*ones(size(x));
%      errorbar2(x,y,e)
%   draws symmetric error bars of unit standard deviation.

if (nargin < 3)
    error('errorbar2 must be called with at least three input arguments');
end

if (nargin == 3)
    w = 0;
    symbol = '-';
elseif (nargin == 4)
    if (isstr(w))
	symbol = w;
	w = 0;
    else
	symbol = '-';
    end
end

plottype = [];
colorvector = [];
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'plottype'
                plottype = varargin{option+1};
                varargin{option+1} = [];
                varargin{option} = [];
            case 'Color'
                colorvector = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

if min(size(x))==1,
  npt = length(x);
  x = x(:);
  y = y(:);
  w = w(:);
  if (min(size(e)) == 2)
      l = e(1,:);
      u = e(2,:);
      l = l(:);
      u = u(:);
  else
      e = e(:);
      u = e;
      l = e;
   end
else
  [npt,n] = size(x);
  u = e;
  l = e;
end

[m,n] = size(y);

u = abs(u);
l = abs(l);
    
if isstr(x) | isstr(y) | isstr(u) | isstr(l)
    error('Arguments must be numeric.')
end

if ~isequal(size(x),size(y)) | ~isequal(size(x),size(l)) | ~isequal(size(x),size(u)),
  error('The sizes of X, Y, L and U must be the same.');
end

%tee = (max(x(:))-min(x(:)));  % make tee .02 x-distance for error bars
if (w > 0) 
    tee = w / 2;
elseif (length(x) > 1)
    tmp = sort(x);
    tee = mean(tmp(2:end) - tmp(1:end-1)) / 6;
else
    % assume that the points are separate by one unit and make the t half a
    % unit
    tee = .5;
end
xl = x - tee;
xr = x + tee;
ytop = y + u;
ybot = y - l;
n = size(y,2);

% Plot graph and bars
hold_state = ishold;
cax = newplot;
next = lower(get(cax,'NextPlot'));

% build up nan-separated vector for bars
xb = zeros(npt*9,n);
xb(1:9:end,:) = x;
xb(2:9:end,:) = x;
xb(3:9:end,:) = NaN;
xb(4:9:end,:) = xl;
xb(5:9:end,:) = xr;
xb(6:9:end,:) = NaN;
xb(7:9:end,:) = xl;
xb(8:9:end,:) = xr;
xb(9:9:end,:) = NaN;

yb = zeros(npt*9,n);
yb(1:9:end,:) = ytop;
yb(2:9:end,:) = ybot;
yb(3:9:end,:) = NaN;
yb(4:9:end,:) = ytop;
yb(5:9:end,:) = ytop;
yb(6:9:end,:) = NaN;
yb(7:9:end,:) = ybot;
yb(8:9:end,:) = ybot;
yb(9:9:end,:) = NaN;

[ls,col,mark,msg] = colstyle(symbol); if ~isempty(msg), error(msg); end
symbol = [ls mark col]; % Use marker only on data part
esymbol = ['-' col]; % Make sure bars are solid

if isempty(plottype) && isempty(colorvector)
    h = plot(xb,yb,esymbol, varargin{1:end}); hold on
elseif isempty(colorvector)
    h = eval([plottype,'(xb,yb,esymbol, varargin{1:end})']); hold on
elseif isempty(plottype)
    h = plot(xb,yb,esymbol,'Color',colorvector); hold on
else
    h = eval([plottype,'(xb,yb,esymbol, Color, [',num2str(colorvector),'])']); hold on
end
%h = [h;plot(x,y,symbol)]; 

if ~hold_state, hold off; end

if nargout>0, hh = h; end
