function [transform, input_points] = align_rectangle(rawpos)
%ALIGN_CYLINDER Register video frame image to rectangular sleep chamber
%
%   This is a hack. The rectangle dimensions are hard-coded into the m-file!
%
%Depends on:
%   CP2TFORM (MATLAB Image Processing Toolbox)
%   IMTRANSFORM (MATLAB Image Processing Toolbox)
%   IS_RAWPOS (written by smk)
%
%Written by smk 2009 October 10.
%

if (exist('cp2tform') ~= 2)
  error(['This function depends on the m-file CP2TFORM ' ...
      '(MATLAB Image Processing Toolbox)']);
end
if (exist('imtransform') ~= 2)
  error(['This function depends on the m-file IMTRANSFORM ' ...
      '(MATLAB Image Processing Toolbox)']);
end
if (exist('is_rawpos') ~= 2)
  error('This function depends on the m-file IS_RAWPOS (written by smk)');
end

if ~is_rawpos(rawpos) || ~isscalar(rawpos)
  error('RAWPOS must be a valid raw position data scalar struct');
end

% EDGES is a struct array in which each element corresponds to a linear
% edge of the rectangle. The user selects points along each edge, which are then
% interpolated/resampled to produce control points for image registration. The
% start->end order of traversal along each edge is important. The fields xstart,
% ystart, xend, yend are expressed in centimeters relative to the center of the
% box.
% The floor is actually smaller than the dimesions given here, but the top edge
% of the rectangle is closer to the camera than the rat's head so these
% artificial dimensions give a better result.
EDGES = cell2struct({ ...
    -35 , 22  , 35  , 22  ; ...
    35  , 22  , 35  , -22 ; ...
    35  , -22 , -35 , -22 ; ...
    -35 , -22 , -35 , 22  }, ...
    {'xstart','ystart','xend','yend'},2);

figure('Units','normalized');
% make the mouse pointer an open circle so that you can see the click location
set(gcf,'Pointer','circle');
% draw a schematic diagram of the rectangle (and keep handles for later
% manipulation
diagram_axes = axes('Units','normalized','Position',[0.76 0 0.23 1], ...
    'DataAspectRatio',[1 1 1],'Visible','off');
for i = 1:numel(EDGES)
  EDGES(i).line_handle = line( ...
      [EDGES(i).xstart EDGES(i).xend], ...
      [EDGES(i).ystart EDGES(i).yend],'Color','k');
end
% text labels (initially invisible)
start_label_handle = text(0,0,'start','Color','m','Visible','off');
end_label_handle = text(0,0,'end','Color','m','Visible','off');

% draw the input image
image_axes = axes('Units','normalized','Position',[0.01 0 0.74 1]);
h = image('Parent',image_axes,'CData',rawpos.video_frame);
set(image_axes,'DataAspectRatio',[1 1 1],'YDir','normal','Visible','off');

% draw video tracking
line(rawpos.xfront,rawpos.yfront,'Color','c','Marker','.','LineStyle','none');

% now the user clicks points along the edges
input_points = zeros(0,2);
base_points = zeros(0,2);
SPLINE_INTERP_FACTOR = 20;
for i = 1:numel(EDGES)
  % highlight the current edge and annotate
  set(EDGES(i).line_handle,'Color','m');
  set(start_label_handle,'Visible','on', ...
      'Position',[EDGES(i).xstart EDGES(i).ystart 0]);
  set(end_label_handle,'Visible','on', ...
      'Position',[EDGES(i).xend EDGES(i).yend 0]);
  % user specifies control points along edge
  [u_cp, v_cp] = trace_edge(image_axes,SPLINE_INTERP_FACTOR);
  % revert to normal color and delete labels
  set(EDGES(i).line_handle,'Color','k');
  set(start_label_handle,'Visible','off');
  set(end_label_handle,'Visible','off');
  % compute cumulative arc length along the interpolating spline that connects
  % the points
  numcp = numel(u_cp);
  cp_interp = interp1((1:numcp)',[u_cp(:), v_cp(:)], ...
      linspace(1,numcp,1+SPLINE_INTERP_FACTOR*(numcp-1)),'spline');
  d_interp = [0; cumsum(hypot(diff(cp_interp(:,1)),diff(cp_interp(:,2))))];
  % the distances of the control points along the interpolating spline
  d = d_interp(1+SPLINE_INTERP_FACTOR*(0:(numcp-1)));
  % we *assume* that the transformation approximately preserves the relative
  % distances of the control points along each edge 
  x_cp = EDGES(i).xstart + ...
      d(:) / d(end) * (EDGES(i).xend - EDGES(i).xstart);
  y_cp = EDGES(i).ystart + ...
      d(:) / d(end) * (EDGES(i).yend - EDGES(i).ystart);
  input_points = [input_points; horzcat(u_cp,v_cp)];
  base_points = [base_points; horzcat(x_cp,y_cp)];
end
delete(gcf);

% infer the transformation
cp2tform_args = {'affine'};
try
  transform = cp2tform(input_points,base_points,cp2tform_args{:});
catch
  error(['CP2TFORM failed with the given arguments. Perhaps you ' ...
      'specified transform_args incorrectly?']);
end

% apply the transformation to the input image and view the result
figure();
[image_out, x_out, y_out] = imtransform(rawpos.video_frame,transform, ...
    'XYScale',0.25);
image(x_out,y_out,image_out);
set(gca,'DataAspectRatio',[1 1 1],'YDir','normal');

end %end main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u,v] = trace_edge(image_axes,SPLINE_INTERP_FACTOR)
%TRACE_EDGE User-interactive selection of points along edge of box
  POINT_DENSITY = 1/40;
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

