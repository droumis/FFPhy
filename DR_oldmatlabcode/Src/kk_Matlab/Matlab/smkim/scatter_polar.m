function hh = scatter_polar(varargin)
%SCATTER_POLAR Scatter/bubble plot of complex values on polar axes.
%
%   SCATTER_POLAR(RLIM,Z,S,C) displays colored circles at the locations
%   specified by the vector Z containing complex values on polar axes.
%
%   RLIM sets the radius of the polar axes. RLIM can be a positive real finite
%   scalar, or the empty matrix []. If RLIM is empty, then the data range is
%   used to determine the axes radius.
%
%   S determines the area of each marker (in points^2). S can be a vector of the
%   same length as Z or a scalar. If S is empty, the default size is used.
%
%   C determines the colors of the markers. When C is a vector of the same
%   length as Z, the values in C are linearly mapped to the colors in the
%   current colormap. When C is a length(Z)-by-3 matrix, it directly specifies
%   the colors of the markers at RGB. C can also be a character color code.
%
%   SCATTER_POLAR(RLIM,Z) draws the markers with default size and color.
%
%   SCATTER_POLAR(RLIM,Z,S) draws the markers with specified size and default
%   color.
%
%   SCATTER_POLAR(...,M) uses the marker M instead of 'o'.
%
%   SCATTER_POLAR(...,'filled') fills the markers.
%
%   SCATTER_POLAR(AX,...) plots into AX instead of GCA.
%
%   H = SCATTER_POLAR(...) returns handles to the scatter objects created.
%
%   This is a wrapper around the function SCATTER_IMPROVED (written by SMK).
%
%Depends on :
%   SCATTER_IMPROVED (written by SMK)
%
%Written by SMK, 2010 March 26.
%

if (exist('scatter_improved') ~= 2)
  error(['SCATTER_POLAR depends on the m-file SCATTER_IMPROVED ' ...
      '(written by smk)']);
end

% Check whether the first argument is an axes handle
[cax,args,nargs] = axescheck(varargin{:});
error(nargchk(1,inf,nargs,'struct'));

% Pop the first two elements of args
rlim = args{1};
z = args{2};
args = args(3:end);

if ~isvector(z) || isreal(z) || ~isfloat(z)
  error('Z must be a vector of complex floating-point values');
end
if any(abs(z(:)) > 1)
  warning('One or more elements of Z has magnitude greater than one');
end
if isempty(rlim) || (isscalar(rlim) && isreal(rlim) && (rlim == +inf))
  rlim = max(abs(z(isfinite(z(:)))));
elseif ~isscalar(rlim) || ~isreal(rlim) || ~isfinite(rlim) || ~(rlim > 0)
  error('RLIM is not valid');
end

% get hold state
cax = newplot(cax);

next = lower(get(cax,'NextPlot'));
hold_state = ishold(cax);

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');
ls = get(cax,'gridlinestyle');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
fUnits  = get(cax, 'DefaultTextUnits');
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
    'DefaultTextFontName',   get(cax, 'FontName'), ...
    'DefaultTextFontSize',   get(cax, 'FontSize'), ...
    'DefaultTextFontWeight', get(cax, 'FontWeight'), ...
    'DefaultTextUnits','data')

% only do grids if hold is off
if ~hold_state

% set axes to be square
    hold(cax,'on');
    hhh=line([-rlim -rlim rlim rlim],[-rlim rlim rlim -rlim],'parent',cax);
    set(cax,'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto')
    v = [get(cax,'xlim') get(cax,'ylim')];
    delete(hhh);
% two radial ticks
    rticks = linspace(0,rlim,3);
    rticks(1) = [];

% define a circle
    th = 0:pi/50:2*pi;
    xunit = cos(th);
    yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = 1:(length(th)-1)/4:length(th);
    xunit(inds(2:2:4)) = zeros(2,1);
    yunit(inds(1:2:5)) = zeros(3,1);
% plot background if necessary
    if ~ischar(get(cax,'color')),
       patch('xdata',xunit*rlim,'ydata',yunit*rlim, ...
             'edgecolor',tc,'facecolor',get(cax,'color'),...
             'handlevisibility','off','parent',cax);
    end

% draw radial circles and label them
    c82 = cos(-45*pi/180);
    s82 = sin(-45*pi/180);

    for i = 1:numel(rticks)
        hhh = line(xunit*rticks(i),yunit*rticks(i), ...
            'linestyle',ls,'color',tc,'linewidth',1,...
                   'handlevisibility','off','parent',cax);
        text(rticks(i)*c82,rticks(i)*s82, ...
            ['  ' num2str(rticks(i))],'verticalalignment','middle',...
            'horizontalalignment','center','handlevisibility','off', ...
            'parent',cax)
    end
    set(hhh,'linestyle','-') % Make outer circle solid

% plot spokes
    th = (1:6)*2*pi/12;
    cst = cos(th); snt = sin(th);
    cs = [-cst; cst];
    sn = [-snt; snt];
    line(rlim*cs,rlim*sn,'linestyle',ls,'color',tc,'linewidth',1,...
         'handlevisibility','off','parent',cax)

% annotate spokes in degrees
    rt = 1.1*rlim;
    for i = 3:3:length(th)
        text(rt*cst(i),rt*snt(i),int2str(i*30),...
             'horizontalalignment','center',...
             'handlevisibility','off','parent',cax);
        if i == length(th)
            loc = int2str(0);
        else
            loc = int2str(180+i*30);
        end
        text(-rt*cst(i),-rt*snt(i),loc,'horizontalalignment','center',...
             'handlevisibility','off','parent',cax)
    end

% set view to 2-D
    view(cax,2);
% set axis limits
    axis(cax,rlim*[-1 1 -1 1]);
end

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle , ...
    'DefaultTextFontName',   fName , ...
    'DefaultTextFontSize',   fSize, ...
    'DefaultTextFontWeight', fWeight, ...
    'DefaultTextUnits',fUnits );


% plot data on top of grid
h = scatter_improved(cax,real(z),imag(z),args{:});

if ~hold_state
    set(cax,'dataaspectratio',[1 1 1]), axis(cax,'off'); set(cax,'NextPlot',next);
end
set(get(cax,'xlabel'),'visible','on')
set(get(cax,'ylabel'),'visible','on')

if ~isempty(h) && ~isdeployed
    makemcode('RegisterHandle',cax,'IgnoreHandle',h,'FunctionName','polar');
end

if (nargout > 0)
  hh = h;
end


