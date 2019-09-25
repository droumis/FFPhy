function [transform, input_points] = align_cylinder(rawpos)
%ALIGN_CYLINDER Register video frame image to cylindrical sleep chamber
%
%   This is a hack. The cylinder dimensions are hard-coded into the m-file!
%
%Depends on:
%   CircleFitByTaubin (Nikolai Chernov, MATLAB Central File Exchange #22678)
%   CP2TFORM (MATLAB Image Processing Toolbox)
%   IMTRANSFORM (MATLAB Image Processing Toolbox)
%   IS_RAWPOS (written by smk)
%
%Written by smk 2009 October 10.
%

if (exist('CircleFitByTaubin') ~= 2)
  error(['This function depends on the m-file CircleFitByTaubin ' ...
      '(MATLAB Central File Exchange #22678)']);
end
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

% Radius of the cylinder in centimeters. The transform converts pixel
% coordinates to coordinates in centimeters relative to the center of the
% cylinder, under the assumption of no camera distortion and isotropic symmetry.
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   Camera distortion is not corrected!
%   The absolute orientation of the cylinder is ignored!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% the floor radius is actually 14 centimeters, but the top edge of the cylinder
% is closer to the camera than the rat's head so 18 centimeters gives a better
% result
CYLINDER_RADIUS = 18; 

figure('Units','normalized');
% make the mouse pointer an open circle so that you can see the click location
set(gcf,'Pointer','circle');
% draw a schematic diagram of the cylinder
diagram_axes = axes('Units','normalized','Position',[0.76 0 0.23 1], ...
    'DataAspectRatio',[1 1 1],'Visible','off');
t = linspace(-pi,pi,500);
line(CYLINDER_RADIUS*cos(t),CYLINDER_RADIUS*sin(t),'Color','k');

% draw the input image
image_axes = axes('Units','normalized','Position',[0.01 0 0.74 1]);
h = image('Parent',image_axes,'CData',rawpos.video_frame);
set(image_axes,'DataAspectRatio',[1 1 1],'YDir','normal','Visible','off');

% draw video tracking
line(rawpos.xfront,rawpos.yfront,'Color','c','Marker','.','LineStyle','none');

% now the user clicks points along the edges
disp(['Click points along the edge of the cylinder in either a clockwise ' ...
    'or counter-clockwise sequence. Do not click points out of order. ' ...
    'Press return when finished (you can adjust the points then)']);
finished = false;
[u, v] = ginput();
control_points = line('XData',u,'YData',v,'Marker','o','LineStyle','none', ...
    'MarkerSize',8,'MarkerEdgeColor','c','MarkerFaceColor','none');
set(control_points,'ButtonDownFcn', ...
    @(varargin) intercept_click(image_axes,control_points));
disp('drag the points to adjust and press [y] to commit');
% set 'CurrentCharacter' property to clear keypress history
set(get(image_axes,'Parent'),'CurrentCharacter',char(0));
while ~finished
  pause(0.01);
  if strcmp('y',get(get(image_axes,'Parent'),'CurrentCharacter'))
    finished = true;
  end
end
delete(gcf);

% fit a circle
tmp = CircleFitByTaubin([u v]);
u_center = tmp(1);
v_center = tmp(2);
uv_radius = tmp(3);

% I don't know how to generate a valid MATLAB Image Processing Toolbox TFORM
% structure by hand, so I do the following circuitous procedure: transform each
% (u,v) point to corresponding (x,y) points, then use CP2TFORM to obtain the
% TFORM that achieves this transformation
input_points = [u, v];
base_points = [CYLINDER_RADIUS/uv_radius*(u - u_center), ...
    CYLINDER_RADIUS/uv_radius*(v - v_center)];

transform = cp2tform(input_points,base_points,'nonreflective similarity');

% apply the transformation to the input image and view the result
figure();
[image_out, x_out, y_out] = imtransform(rawpos.video_frame,transform, ...
    'XYScale',0.25);
image(x_out,y_out,image_out);
set(gca,'DataAspectRatio',[1 1 1],'YDir','normal');

end %end main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested subfunctions for mouse callbacks
function intercept_click(image_axes,control_points)
  HIT_RADIUS = 5;
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


