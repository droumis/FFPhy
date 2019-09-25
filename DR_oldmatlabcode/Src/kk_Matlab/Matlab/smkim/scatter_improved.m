function hh = scatter_improved(varargin)
%SCATTER_IMPROVED Scatter/bubble plot, using LINE instead of PATCH.
%   SCATTER_IMPROVED(X,Y,S,C) displays colored circles at the locations specified
%   by the vectors X and Y (which must be the same size).  
%
%   S determines the area of each marker (in points^2). S can be a
%   vector the same length a X and Y or a scalar. If S is a scalar, 
%   MATLAB draws all the markers the same size. If S is empty, the
%   default size is used.
%   
%   C determines the colors of the markers. When C is a vector the
%   same length as X and Y, the values in C are linearly mapped
%   to the colors in the current colormap. When C is a 
%   length(X)-by-3 matrix, it directly specifies the colors of the  
%   markers as RGB values. C can also be a color string. See ColorSpec.
%
%   SCATTER_IMPROVED(X,Y) draws the markers in the default size and color.
%   SCATTER_IMPROVED(X,Y,S) draws the markers at the specified sizes (S)
%   with a single color. This type of graph is also known as
%   a bubble plot.
%   SCATTER_IMPROVED(...,M) uses the marker M instead of 'o'.
%   SCATTER_IMPROVED(...,'filled') fills the markers.
%
%   SCATTER_IMPROVED(AX,...) plots into AX instead of GCA.
%
%   H = SCATTER_IMPROVED(...) returns handles to the scatter objects created.
%
%See also:
%   SCATTER (Matlab special graphing function)

[cax,args,nargs] = axescheck(varargin{:});
error(nargchk(1,inf,nargs,'struct'));
[pvpairs,args,nargs,msg] = parseargs(args);
error(msg); %#ok
error(nargchk(2,4,nargs,'struct'));

dataargs = datachk(args(1:nargs));
 
switch (nargs)
  case 2
    [x,y] = deal(dataargs{:});
    error(Lxychk(x,y)); %#ok
    [cax,parax] = localGetAxesInfo(cax);
    [ls,c,m] = nextstyle(cax); %#ok
    error(Lcchk(x,c));  %#ok
    s = get(cax,'defaultlinemarkersize')^2;
  case 3
    [x,y,s] = deal(dataargs{:});
    error(Lxychk(x,y));  %#ok
    error(Lschk(x,s));  %#ok
    [cax,parax] = localGetAxesInfo(cax);
    [ls,c,m] = nextstyle(cax); %#ok
    error(Lcchk(x,c));  %#ok
  case 4
    [x,y,s,c] = deal(dataargs{:});
    error(Lxychk(x,y));  %#ok
    error(Lschk(x,s));  %#ok
    if ischar(args{nargs}), c = args{nargs}; end
    error(Lcchk(x,c));  %#ok
    [cax,parax] = localGetAxesInfo(cax);
end

if length(x) ~= length(y) || ...
    length(x) ~= numel(x) || length(y) ~= numel(y)
  error(id('InvalidData'),'X and Y must be vectors of the same length.');
end

% Map colors into colormap colors if necessary.  Reshape so as to
% make the loop below easier to follow
[color, scaled] = MapColorsToColorMap(x, c);
colorSize = size(color);
if isequal(colorSize,[1 3]) || ischar(color),
    color = repmat(color,length(x),1);
elseif any(colorSize == 1),
    color = color(:);
end

% Scalar expand the marker size if necessary
if isempty(s)
  s = pi*5^2;
end
if length(s)==1
  s = repmat(s,length(x),1); 
elseif length(s)~=numel(s) || length(s)~=length(x)
  error(id('InvalidSData'),'S must be a scalar or a vector the same length as X.')
end

% create an invisible handle invisible axes for temporary parent
fig = ancestor(cax,'figure');
curax = get(fig,'currentaxes');
tax = axes('parent',fig,'visible','off','handlevisibility','off');
h = -1; h = h(ones(length(x),1));
for i = 1:length(x)
    specs = pvpairs;
    j = find(cellfun(@(c) isequal('MarkerFaceColor',c),pvpairs));
    if ~isempty(j) && (j+1 <= numel(pvpairs))
      specs{j+1} = color(i,:);
    end
    % Notice that marker 'o' gets overridden if specs{:} specifies a different
    % marker
    h(i) = line('parent',tax,'xdata',x(i),'ydata',y(i),...
        'Marker','o',specs{:},'linestyle','none', ...
        'markersize',sqrt(s(i)/pi),'markeredgecolor',color(i,:));
end
set(h,'parent',cax);
delete(tax);
set(fig,'currentaxes',curax);

if (nargout > 0)
  hh = h;
end

%--------------------------------------------------------------------------
function [cax,parax] = localGetAxesInfo(cax)

if isempty(cax) || isa(handle(cax),'hg.axes') 
    cax = newplot(cax);
    parax = cax;
else
    parax = cax;
    cax = ancestor(cax,'Axes');
end


%--------------------------------------------------------------------------
function [pvpairs,args,nargs,msg] = parseargs(args)
msg = '';
% separate pv-pairs from opening arguments
[args,pvpairs] = parseparams(args);
n = 1;
extrapv = {};
% check for 'filled' or LINESPEC or ColorSpec
while length(pvpairs) >= 1 && n < 4 && ischar(pvpairs{1})
  arg = lower(pvpairs{1});
  if arg(1) == 'f'
    pvpairs(1) = [];
    extrapv = {'MarkerFaceColor','filled','MarkerEdgeColor','none', ...
               extrapv{:}};
  else
    [l,c,m,tmsg]=colstyle(pvpairs{1});
    if isempty(tmsg)
      pvpairs(1) = [];
      if ~isempty(l) 
        extrapv = {'LineStyle',l,extrapv{:}};
      end
      if ~isempty(c)
        extrapv = {'CData',ColorSpecToRGB(c),extrapv{:}};
      end
      if ~isempty(m)
        extrapv = {'Marker',m,extrapv{:}};
      end
    end
  end
  n = n+1;
end
pvpairs = [extrapv pvpairs];
if isempty(args)
  msg.message = 'Must supply X and Y data as first arguments.';
  msg.identifier = id('NoDataInputs');
else
  msg = checkpvpairs(pvpairs);
end
nargs = length(args);

%--------------------------------------------------------------------------

function color = ColorSpecToRGB(s)
color=[];
switch s
 case 'y'
  color = [1 1 0];
 case 'm'
  color = [1 0 1];
 case 'c'
  color = [0 1 1];
 case 'r'
  color = [1 0 0];
 case 'g'
  color = [0 1 0];
 case 'b'
  color = [0 0 1];
 case 'w'
  color = [1 1 1];
 case 'k'
  color = [0 0 0];
end

%--------------------------------------------------------------------------

function [colorMap, scaled] = MapColorsToColorMap(x,c)
% Map colors into colormap colors if necessary.

scaled = false;
if ischar(c) || isequal(size(c),[1 3]); % string color or scalar rgb
    colorMap = c;
elseif length(c)==numel(c) && length(c)==length(x), % is C a vector?
    scaled = true;
    colorMap = c;
elseif isequal(size(c),[length(x) 3]), % vector of rgb's
    colorMap = c;
else
    error(id('InvalidCData'),...
    'C must be a single color, a vector the same length as X, or an M-by-3 matrix.')
end

%--------------------------------------------------------------------------
function msg = Lxychk(x,y)
msg = [];
% Verify {X,Y) data is correct size
if any([length(x) length(y) ...
        numel(x) numel(y) ] ~= length(x))
    msg = struct('identifier',id('InvalidXYData'),...
                 'message','X and Y must be vectors of the same length.');
end

 %--------------------------------------------------------------------------
function msg = Lcchk(x,c)
msg = [];
% Verify CData is correct size
if ischar(c) || isequal(size(c),[1 3]); 
    % string color or scalar rgb 
elseif length(c)==numel(c) && length(c)==length(x)
    % C is a vector
elseif isequal(size(c),[length(x) 3]), 
    % vector of rgb's
else
    msg = struct('identifier',id('InvalidCData'),...
                 'message',['C must be a single color, a vector the same length as X, ',...
           'or an M-by-3 matrix.']);
end

%--------------------------------------------------------------------------
function msg = Lschk(x,s)
msg = [];
% Verify correct S vector
if length(s) > 1 && ...
              (length(s)~=numel(s) || length(s)~=length(x))
    msg = struct('identifier',id('InvalidSData'),...
                 'message','S must be a scalar or a vector the same length as X.');
end

%--------------------------------------------------------------------------
function str = id(str)
str = ['MATLAB:scatter:' str];

%--------------------------------------------------------------------------
function msg = checkpvpairs(pvpairs,linestyle)
%CHECKPVPAIRS Check length of property value pair inputs
%   MSG = CHECKPVPAIRS(PAIRS) returns an error message if the
%   length of the supplied cell array of property value pairs is
%   incorrect. The message is tailored to accepting LINESPEC as an
%   option convenience input argument.
msg = '';
npvpairs = length(pvpairs)/2;
if nargin == 1, linestyle = true; end
if (length(pvpairs) == 1) && linestyle
  msg = struct('message','Error in color/linetype argument.',...
               'identifier',id('UnrecognizedLineStyleParameter'));
elseif npvpairs ~= fix(npvpairs)
  msg = struct('message','Incorrect number of inputs for property-value pairs.',...
               'identifier',id('EvenPropValuePairs'));
end

%--------------------------------------------------------------------------
function y = datachk(x)
%DATACHK Convert input to full, double data for plotting
%  Y=DATACHK(X) creates a full, double array from X and returns it in Y.
%  If X is a cell array each element is converted elementwise.
if iscell(x)
    y = cellfun(@datachk,x,'UniformOutput',false);
else
    y = full(double(x));
end

