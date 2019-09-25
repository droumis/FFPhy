function [transform, input_points] = align_Utrack(rawpos)
%ALIGN_UTRACK Register video frame image to U track
%
%   This is a hack. The U-track dimensions are hard-coded into the m-file!
%
%Depends on:
%   CP2TFORM (MATLAB Image Processing Toolbox)
%   MAKETFORM (MATLAB Image Processing Toolbox)
%   TFORMFWD (MATLAB Image Processing Toolbox)
%   TFORMINV (MATLAB Image Processing Toolbox)
%   IMTRANSFORM (MATLAB Image Processing Toolbox)
%   FMINBND (MATLAB Optimization Toolbox)
%   TRICUBE (written by smk)
%   BISQUARE (written by smk)
%   IS_RAWPOS(written by smk)
%
%Written by smk 2009 October 30.
%

if (exist('cp2tform') ~= 2)
  error(['This function depends on the m-file CP2TFORM ' ...
      '(MATLAB Image Processing Toolbox)']);
end
if (exist('imtransform') ~= 2)
  error(['This function depends on the m-file IMTRANSFORM ' ...
      '(MATLAB Image Processing Toolbox)']);
end
if (exist('tformfwd') ~= 2)
  error(['This function depends on the m-file TFORMFWD ' ...
      '(MATLAB Image Processing Toolbox)']);
end
if (exist('tforminv') ~= 2)
  error(['This function depends on the m-file TFORMINV ' ...
      '(MATLAB Image Processing Toolbox)']);
end
if (exist('maketform') ~= 2)
  error(['This function depends on the m-file MAKETFORM ' ...
      '(MATLAB Image Processing Toolbox)']);
end
if (exist('fminbnd') ~= 2)
  error(['This function depends on the m-file FMINBND ' ...
      '(MATLAB Optimization Toolbox)']);
end
if (exist('tricube') ~= 2)
  error('This function depends on the m-file TRICUBE (written by smk)');
end
if (exist('bisquare') ~= 2)
  error('This function depends on the m-file BISQUARE (written by smk)');
end
if (exist('is_rawpos') ~= 2)
  error('This function depends on the m-file IS_RAWPOS (written by smk)');
end

if ~is_rawpos(rawpos) || ~isscalar(rawpos)
  error('RAWPOS must be a valid raw position data scalar struct');
end

% SEGMENTS is a struct array in which each element corresponds to a linear
% segment the Utrack. The user selects points along each edge, which are then
% interpolated/resampled to produce control points for image registration. The
% start->end order of traversal along each edge is important. The fields xstart,
% ystart, xend, yend are expressed in centimeters relative to the location of
% one of the food wells.
SEGMENTS = cell2struct({ ...
    0   , 147 , 0   , 0 ; ...
    0   , 0   , 44  , 0 ; ...
    44  , 0   , 44  , 147   }, ...
    {'xstart','ystart','xend','yend'},2);

figure('Units','normalized');
% make the mouse pointer an open circle so that you can see the click location
set(gcf,'Pointer','circle');
% draw a schematic diagram of the rectangle (and keep handles for later
% manipulation
diagram_axes = axes('Units','normalized','Position',[0.76 0 0.20 1], ...
    'DataAspectRatio',[1 1 1],'Visible','off');
for i = 1:numel(SEGMENTS)
  SEGMENTS(i).line_handle = line( ...
      [SEGMENTS(i).xstart SEGMENTS(i).xend], ...
      [SEGMENTS(i).ystart SEGMENTS(i).yend],'Color','k','LineWidth',4);
end
% text labels (initially invisible)
start_label_handle = text(0,0,'start','Color','m','Visible','off');
end_label_handle = text(0,0,'end','Color','m','Visible','off');
% draw the input image
image_axes = axes('Units','normalized','Position',[0.01 0 0.74 1]);
h = image('Parent',image_axes,'CData',rawpos.video_frame);
set(image_axes,'DataAspectRatio',[1 1 1],'YDir','normal','Visible','off');

% now the user clicks points along the edges
input_points = zeros(0,2);
base_points = zeros(0,2);
SPLINE_INTERP_FACTOR = 20;
for i = 1:numel(SEGMENTS)
  % highlight the current edge and annotate
  set(start_label_handle,'Visible','on', ...
      'Position',[SEGMENTS(i).xstart SEGMENTS(i).ystart 0]);
  set(end_label_handle,'Visible','on', ...
      'Position',[SEGMENTS(i).xend SEGMENTS(i).yend 0]);
  set(SEGMENTS(i).line_handle,'Color','m');
  % user specifies control points along edge
  [u_cp, v_cp] = trace_edge(image_axes,SPLINE_INTERP_FACTOR);
  % revert to normal color and delete labels
  set(SEGMENTS(i).line_handle,'Color','k');
  set(start_label_handle,'Visible','off');
  set(end_label_handle,'Visible','off');
  % compute cumulative arc length along the interpolating spline that connects
  % the points
  numcp = numel(u_cp);
  cp_interp = interp1((1:numcp)',[u_cp, v_cp], ...
      linspace(1,numcp,1+SPLINE_INTERP_FACTOR*(numcp-1)),'spline');
  d_interp = [0; cumsum(hypot(diff(cp_interp(:,1)),diff(cp_interp(:,2))))];
  % the distances of the control points along the interpolating spline
  d = d_interp(1+SPLINE_INTERP_FACTOR*(0:(numcp-1)));
  % we *assume* that the transformation approximately preserves the relative
  % distances of the control points along each edge 
  x_cp = SEGMENTS(i).xstart + ...
      d(:) / d(end) * (SEGMENTS(i).xend - SEGMENTS(i).xstart);
  y_cp = SEGMENTS(i).ystart + ...
      d(:) / d(end) * (SEGMENTS(i).yend - SEGMENTS(i).ystart);
  input_points = [input_points; horzcat(u_cp,v_cp)];
  base_points = [base_points; horzcat(x_cp,y_cp)];
end
delete(gcf);

% infer the transformation
cp2tform_args = {'polynomial',2};
try
  transform = cp2tform(input_points,base_points,cp2tform_args{:});
catch
  error(['CP2TFORM failed with the given arguments. Perhaps you ' ...
      'specified transform_args incorrectly?']);
end

% view the transformed image
figure();
subplot(1,2,1);
[image_out, x_out, y_out] = imtransform(rawpos.video_frame,transform, ...
    'XYScale',0.5);
image(x_out,y_out,image_out);
for i = 1:numel(SEGMENTS)
  SEGMENTS(i).line_handle = line( ...
      [SEGMENTS(i).xstart SEGMENTS(i).xend], ...
      [SEGMENTS(i).ystart SEGMENTS(i).yend],'Color','m');
end
set(gca,'YDir','normal','DataAspectRatio',[1 1 1]);

% Check whether transform has a defined forward_fcn
if isa(transform.forward_fcn,'function_handle')
  direct_transform = true;
else
  direct_transform = false;
end

nanflag = isnan(rawpos.xfront) | isnan(rawpos.yfront) | ...
    isnan(rawpos.xback) | isnan(rawpos.yback);
u = double(rawpos.xfront(~nanflag));
v = double(rawpos.yfront(~nanflag));

if (direct_transform)
  [x, y] = tformfwd(transform,u,v);
else
  % construct lookup table in the inverse direction, then interpolate in the
  % forward direction
  [xgrid, ygrid] = meshgrid(-100:144,-100:244);
  [ugrid, vgrid] = tforminv(transform,xgrid,ygrid);
  % x and y are expressed in real-world coordinates with units of centimeters
  x = griddata(ugrid,vgrid,xgrid,u,v,'cubic');
  y = griddata(ugrid,vgrid,ygrid,u,v,'cubic');
end

% resample uniformly at 1 cm spacing along path length
path_length = [0; cumsum(hypot(diff(x),diff(y)))];
% remove points at which path_length has duplicate values, otherwise interp1
% will fail
zero_length = find(diff(path_length) < max(eps(path_length)));
path_length(zero_length) = [];
x(zero_length) = [];
y(zero_length) = [];
path_length_resampled = (0:ceil(path_length(end)))';
x_resampled = interp1(path_length,x,path_length_resampled,'linear','extrap');
y_resampled = interp1(path_length,y,path_length_resampled,'linear','extrap');
x_diff = [0; diff(x_resampled)];
y_diff = [0; diff(y_resampled)];

% estimate x_offset from the path segments when the y-component of motion is
% large (i.e. motion along the vertical arms of the U-track)
x_offset = find_best_x_offset(x_resampled,tricube(1-abs(y_diff)));
% estimate y_offset from the path segments when the x-component of motion is
% large (i.e. motion along the horizontal "crossbar" arm of the U-track)
y_offset = find_best_y_offset(y_resampled,tricube(1-abs(x_diff)));
% compose this offset with the camera distortion correction transform; the
% composition is not commutative so be careful to do it in the correct order
corners = [0, 0; 1, 1];
transform = maketform('composite', ...
    maketform('box',corners,corners+repmat([x_offset y_offset],[2 1])), ...
    transform);

% apply the transformation to the input and view the result
if (direct_transform)
  [x, y] = tformfwd(transform,double(rawpos.xfront),double(rawpos.yfront));
else
  % construct lookup table in the inverse direction, then interpolate in the
  % forward direction
  [xgrid, ygrid] = meshgrid(-100:144,-100:244);
  [ugrid, vgrid] = tforminv(transform,xgrid,ygrid);
  % x and y are expressed in real-world coordinates with units of centimeters
  x = griddata(ugrid,vgrid,xgrid, ...
      double(rawpos.xfront),double(rawpos.yfront),'cubic');
  y = griddata(ugrid,vgrid,ygrid, ...
      double(rawpos.xfront),double(rawpos.yfront),'cubic');
end
subplot(1,2,2);
line(x,y,'Color','k');
for i = 1:numel(SEGMENTS)
  SEGMENTS(i).line_handle = line( ...
      [SEGMENTS(i).xstart SEGMENTS(i).xend], ...
      [SEGMENTS(i).ystart SEGMENTS(i).yend],'Color','m');
end
set(gca,'YDir','normal','DataAspectRatio',[1 1 1]);

end % end main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u,v] = trace_edge(image_axes,SPLINE_INTERP_FACTOR)
%TRACE_EDGE User-interactive selection of points along edge of box
  POINT_DENSITY = 1/20;
  HIT_RADIUS = 5;
  finished = false;
  disp(['click on start and end vertices of the edge in the indicated ' ...
      'order (you can change your selection later']);
  while ~finished
    [endpoints_x, endpoints_y] = ginput(2);
    dist = hypot(diff(endpoints_x),diff(endpoints_y));
    numvertices = max(4,round(dist * POINT_DENSITY));
    % generate a control_points that goes between the endpoints
    control_points = line('Parent',image_axes,'Visible','off','Marker','o', ...
        'MarkerSize',8,'MarkerEdgeColor','c','MarkerFaceColor','none', ...
        'Color','c','LineStyle',':', ...
        'XData',linspace(endpoints_x(1),endpoints_x(end),numvertices)', ...
        'YData',linspace(endpoints_y(1),endpoints_y(end),numvertices)');
    % make the vertices of this control_points draggable
    set(control_points,'ButtonDownFcn', ...
        @(varargin) intercept_click(image_axes,control_points));
    % reveal the control_points
    set(control_points,'Visible','on','HitTest','on');
    disp(['drag vertices to adjust or press [y] to commit or [n] to ' ...
        'clear and redo']);
    while true
      pause(0.01);
      if strcmp('n',get(get(image_axes,'Parent'),'CurrentCharacter'))
        % clear and redo
        delete(control_points);
        % set 'CurrentCharacter' property to clear keypress history
        set(get(image_axes,'Parent'),'CurrentCharacter',char(0));
        break;
      elseif strcmp('y',get(get(image_axes,'Parent'),'CurrentCharacter'))
        % compute the spline interpolation that passes through the points
        u = get(control_points,'XData');
        u = u(:);
        v = get(control_points,'YData');
        v = v(:);
        assert(isequal(numel(u),numel(v),numvertices));
        tmp_interp = interp1((1:numvertices)',[u, v], ...
            linspace(1,numvertices,SPLINE_INTERP_FACTOR*(numvertices-1)+1), ...
            'spline');
        u_interp = tmp_interp(:,1);
        v_interp = tmp_interp(:,2);
        delete(control_points);
        line('XData',u,'YData',v,'Color','c','LineStyle','none', ...
            'Marker','o','MarkerSize',8,'MarkerFaceColor','none');
        line('XData',u_interp,'YData',v_interp,'Color','c', ...
            'LineStyle',':');
        finished = true;
        % set 'CurrentCharacter' property to clear keypress history
        set(get(image_axes,'Parent'),'CurrentCharacter',char(0));
        break;
      else
        % user is clicking/dragging the impoly object or pressing irrelevant
        % keys; do not interfere
      end
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Nested subfunctions for mouse callbacks
  function intercept_click(image_axes,control_points)
    if strcmp(get(gcf,'SelectionType'),'normal')
      cp = get(image_axes,'CurrentPoint');
      [mindist, nearest_idx] = min(hypot( ...
          cp(1,1) - get(control_points,'XData'), ...
          cp(1,2) - get(control_points,'YData')));
      if ~isempty(nearest_idx) && (mindist < HIT_RADIUS)
        set(control_points,'ButtonDownFcn','');
        set(gcf,'WindowButtonMotionFcn', ...
            @(varargin) drag_vertex(image_axes,control_points,nearest_idx));
        set(gcf,'WindowButtonUpFcn', ...
            @(varargin) stop_dragging(image_axes,control_points));
      end
    end
  end
  function drag_vertex(image_axes,control_points,idx)
    cp = get(image_axes,'CurrentPoint');
    xdata = get(control_points,'XData');
    ydata = get(control_points,'YData');
    xdata(idx) = cp(1,1);
    ydata(idx) = cp(1,2);
    set(control_points,'XData',xdata,'YData',ydata);
    drawnow;
  end
  function stop_dragging(image_axes,control_points)
    set(control_points,'ButtonDownFcn', ...
        @(varargin) intercept_click(image_axes,control_points));
    set(gcf,'WindowButtonMotionFcn','');
    set(gcf,'WindowButtonUpFcn','');
  end
end % end nested subfunction TRACE_EDGE

function x_offset = find_best_x_offset(x,weights)
  % this parameter corresponds to the width-scale of the trajectory along the
  % U-track
  EPSILON = 12;  
  % select the points with the smallest residuals from zero and use them to
  % compute an initial estimate
  residuals = min(abs(x - 0),abs(x - 44));
  idx = find(residuals <= quantile(residuals,0.1));
  x_offset = fminbnd(@(arg) sum(weights(idx) .* min( ...
      abs(x(idx) + arg - 0),abs(x(idx) + arg - 44))), -10,10);
  % now use this initial estimate as a seed for iterative fitting
  for n = 1:3
    % compute residuals and reweight accordingly
    residuals = min(abs(x + x_offset - 0),abs(x + x_offset - 44));
    weights = weights .* bisquare(residuals/EPSILON);
    x_offset = fminbnd(@(arg) sum(weights .* min( ...
        abs(x + arg - 0),abs(x + arg - 44))), -10,10);
  end
end % end subfunction FIND_BEST_X_OFFSET

function y_offset = find_best_y_offset(y,weights)
  % this parameter corresponds to the width-scale of the trajectory along the
  % U-track
  EPSILON = 12;  
  % select the points with the smallest residuals from zero and use them to
  % compute an initial estimate
  residuals = abs(y - 0);
  idx = find(residuals <= quantile(residuals,0.1));
  y_offset = fminbnd(@(arg) sum(weights(idx) .* abs(y(idx) + arg - 0)), -10,10);
  % now use this initial estimate as a seed for iterative fitting
  for n = 1:3
    % compute residuals and reweight accordingly
    residuals = abs(y + y_offset - 0);
    weights = weights .* bisquare(residuals/EPSILON);
    y_offset = fminbnd(@(arg) sum(weights .* abs(y + arg - 0)), -10,10);
  end
end % end subfunction FIND_BEST_Y_OFFSET


