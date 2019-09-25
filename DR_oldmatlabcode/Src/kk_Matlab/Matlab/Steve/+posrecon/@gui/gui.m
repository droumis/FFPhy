classdef gui < hgsetget

  properties (SetAccess = public, GetAccess = public)
    % these public properties are linked to set methods that automatically
    % update the display
    picture
    mask
    frame_number
    timestamp
    front_epsilon
    back_epsilon
    % we specify default values for front_coordinates and back_coordinates so
    % that initial set calls for front_epsilon and back_epsilon go smoothly
    front_coordinates = [NaN NaN];
    back_coordinates = [NaN NaN];
    threshold;
    % automode is an enum (see controller for details of the modes)
    %   1: playback only
    %   2: make automated assignments when possible
    %   3: clear all assignments in frames encountered during playback
    automode = 1;
    enable_input = true;
    progress_fraction

    debug
  end

  properties (SetAccess = private, GetAccess = public)
    % dummy variables for parametric drawing of circles
    dummy_x = cos(linspace(-pi,pi,100));
    dummy_y = sin(linspace(-pi,pi,100));
    % cdata element for the pixel selection mask
    mask_color = uint8([intmax('uint8') 0 0]);
    % opacity of the pixel selection mask
    mask_opacity = uint8(0.75*intmax('uint8'));
    % colors for the markup indicating the current position estimate
    front_color = 'c';
    back_color = 'y';
    line_color = 'm';
    % parent figure
    figure
    % the axes that holds picture_container as well as any overlays
    picture_axes
    % the image object that displays the current picture
    picture_container
    % the image object that displsy the current mask
    mask_container
    % a struct containing handles to line objects for drawing markers to
    % indicate the estimated position of the rat
    markers_container
    % a vector of handles to uimenu objects (menu buttons)
    menu_handles
    % an 'edit' style uicontrol object for displaying a frame number and
    % allowing the user to input a new value
    frame_counter
    % a (static) 'text' style uicontrol object that displays a timestamp
    clock
    % a (static) 'text' style uicontrol object that reports the coordinates of
    % the front marker
    front_coordinates_display
    % a (static) 'text' style uicontrol object that reports the coordinates of
    % the back marker
    back_coordinates_display
    % an 'edit' style uicontrol object for displaying the scale of the front
    % marker and allowing the user to input a new value
    front_epsilon_display
    % an 'edit' style uicontrol object for displaying the scale of the back
    % marker and allowing the user to input a new value
    back_epsilon_display
    % an 'edit' style uicontrol object for displaying the current threshold and
    % allowing the user to input a new value
    threshold_display
    % a uibuttongroup object containing radio buttons for toggling between
    % playback/assignment modes
    automode_radiobox
    % a vector of handles to the uicontrol radio buttons that are children of
    % automode_radiobox
    radio_handles
    % a (static) 'text' style uicontrol that reports the proportion of frames
    % that have assigned position estimates
    progress_meter
    % vector of handles to overlay graphics objects, which are plotted on
    % picture_axes. the public methods add_overlay, update_overlay,
    % remove_overlay allow users to plot arbitrary overlays on the axes
    overlays = zeros([1000 1]);
    % upon creation (by calling the add_overlay method), overlays are assigned a
    % simple integer identifier which is used to index into the overlays vector.
    % next_overlay_identifier is the value of the unique identifier to be issued
    % on the next add_overlay request.
    next_overlay_identifier = 1;
  end

  events
    user_event
  end

  methods (Access = public)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constructor
    function self = gui(picture,mask,pixel_aspect_ratio,frame_number, ...
        timestamp,front_coordinates,back_coordinates, ...
        front_epsilon,back_epsilon,threshold, ...
        progress_fraction,automode,enable_input)
      % Check inputs
      if ~isa(picture,'uint8') || ~isreal(picture) || ...
          ~(ndims(picture) == 3) || ~(size(picture,3) == 3)
        error(['picture argument must be a 3-dimensional real uint8 array ' ...
            'whose third dimension is length 3']);
      end
      if ~islogical(mask) || ~(ndims(mask) == 2) || ...
          ~(size(mask,1) == size(picture,1)) || ...
          ~(size(mask,2) == size(picture,2))
        error(['mask argument must be a 2-dimensional logical array ' ...
            'whose size along the first and second dimensions equals ' ...
            'the size along the corresponding dimensions of picture']);
      end
      if ~isnumeric(pixel_aspect_ratio) || ~isscalar(pixel_aspect_ratio) || ...
          ~isreal(pixel_aspect_ratio) || ~isfinite(pixel_aspect_ratio) || ...
          ~(pixel_aspect_ratio > 0)
        error('pixel_aspect_ratio argument must be a positive real');
      end
      if ~isnumeric(frame_number) || ~isscalar(frame_number) || ...
          ~isreal(frame_number) || ~isfinite(frame_number) || ...
          ~(frame_number >= 1) || ~(round(frame_number) == frame_number)
        error('frame_number argument must be a positive real integer');
      end
      if (exist('ts2str') ~= 2)
        error('this gui requires the m-file TS2STR (written by smk)');
      end
      try
        ts2str(timestamp);
      catch
        error('timestamp argument must be a valid uint32 timestamp');
      end
      if ~isfloat(front_coordinates) || ~isvector(front_coordinates) || ...
          ~(numel(front_coordinates) == 2) || ~all(~isinf(front_coordinates))
        error(['front_coordinates must be a 2-element vector of non-Inf ' ...
            'floating-point values']);
      end
      if ~isfloat(back_coordinates) || ~isvector(back_coordinates) || ...
          ~(numel(back_coordinates) == 2) || ~all(~isinf(back_coordinates))
        error(['back_coordinates must be a 2-element vector of non-Inf ' ...
            'floating-point values']);
      end
      if ~isfloat(front_epsilon) || ~isscalar(front_epsilon) || ...
          ~isreal(front_epsilon) || ~isfinite(front_epsilon) || ...
          ~(front_epsilon > 0)
        error('front_epsilon must be a positive real scalar');
      end

      if isnan(back_epsilon)
        
      elseif ~isfloat(back_epsilon) || ~isscalar(back_epsilon) || ...
          ~isreal(back_epsilon) || ~isfinite(back_epsilon) || ...
          ~(back_epsilon > 0)
        error('back_epsilon must be a non-negative real scalar');
      end
      if ~isa(threshold,'uint8') || ~isscalar(threshold) || ...
          ~isreal(threshold)
        error('threshold must be a real uint8 scalar');
      end
      if ~isfloat(progress_fraction) || ~isscalar(progress_fraction) || ...
          ~isreal(progress_fraction) || ~(progress_fraction >= 0) || ...
          ~(progress_fraction <= 1)
        error(['progress_fraction must be a positive real scalar between ' ...
            '0 and 1']);
      end
      if ~isnumeric(automode) || ~isscalar(automode) || ...
          ~any(automode == [0 1 2])
        error('automode argument must assume a scalar value of {0,1,2}');
      end
      if ~islogical(enable_input) || ~isscalar(enable_input)
        error('enable_input argument must be a logical scalar');
      end
      % Populate the static gui elements. 'HandleVisibility' property of the
      % main figure must be 'on' for add_overlay method to work.
      self.figure = figure('Name','posrecon','Visible','off', ...
          'NextPlot','add','HandleVisibility','on', ...
          'DockControls','off','MenuBar','none', ...
          'Toolbar','none','Pointer','crosshair','Renderer','OpenGL', ...
          'Units','normalized','Color','k', ...
          'DefaultUipanelHitTest','off', ...
          'DefaultLineHitTest','off', ...
          'DefaultImageHitTest','off', ...
          'DefaultUipanelForegroundcolor','w', ...
          'DefaultUipanelBackgroundcolor','k', ...
          'DefaultUipanelHighlightcolor','k', ...
          'DefaultUipanelShadowcolor','k', ...
          'DefaultUicontrolForegroundcolor','w', ...
          'DefaultUicontrolBackgroundcolor','k', ...
          'KeyPressFcn',@(varargin) self.notify('user_event', ...
          posrecon.user_event(varargin{:})), ...
          'CloseRequestFcn',@(varargin) self.notify('user_event', ...
          posrecon.user_event(varargin{:})));
      button_labels = { ...
          'Quit', ...
          'Export video frame', ...
          'Export assignments', ...
          'Save assignments', ...
          'Load assignments', ...
          'Summary' };
      self.menu_handles = zeros(size(button_labels));
      for i = 1:numel(button_labels)
        self.menu_handles(i) = uimenu('Parent',self.figure, ...
            'Label',button_labels{i});
        set(self.menu_handles(i),'Callback',@(varargin) ...
            self.notify('user_event',posrecon.user_event(varargin{:})));
      end
      control_panel = uipanel('Parent',self.figure, ...
          'Units','normalized','Position',[0.0 0.90 1.0 0.1]);
      self.frame_counter = uicontrol('Parent',control_panel, ...
          'Tag','frame counter','Style','edit','Units','normalized', ...
          'Position',[0.01 0.5 0.18 0.4], ...
          'Callback',@(varargin) self.notify('user_event', ...
          posrecon.user_event(varargin{:})));
      self.clock = uicontrol('Parent',control_panel, ...
          'Style','text','Units','normalized', ...
          'Position',[0.01 0.0 0.18 0.4]);
      self.front_epsilon_display = uicontrol('Parent',control_panel, ...
          'Tag','front epsilon','Style','edit','Units','normalized', ...
          'Position',[0.21 0.5 0.08 0.4], ...
          'String',num2str(front_epsilon), ...
          'Callback',@(varargin) self.notify('user_event', ...
          posrecon.user_event(varargin{:})));
      self.back_epsilon_display = uicontrol('Parent',control_panel, ...
          'Tag','back epsilon','Style','edit','Units','normalized', ...
          'Position',[0.21 0.0 0.08 0.4], ...
          'String',num2str(back_epsilon), ...
          'Callback',@(varargin) self.notify('user_event', ...
          posrecon.user_event(varargin{:})));
      self.front_coordinates_display = uicontrol('Parent',control_panel, ...
          'Style','text','Units','normalized', ...
          'Position',[0.31 0.5 0.18 0.4],'HorizontalAlignment','left', ...
          'ForegroundColor',self.front_color, ...
          'String',sprintf('front (%.1f,%.1f)', ...
          front_coordinates(1),front_coordinates(2)));
      self.back_coordinates_display = uicontrol('Parent',control_panel, ...
          'Style','text','Units','normalized', ...
          'Position',[0.31 0.0 0.18 0.4],'HorizontalAlignment','left', ...
          'ForegroundColor',self.back_color, ...
          'String',sprintf('front (%.1f,%.1f)', ...
          back_coordinates(1),back_coordinates(2)));
      uicontrol('Parent',control_panel,'Style','text', ...
          'Units','normalized','Position',[0.51 0.5 0.13 0.4], ...
          'String','threshold')
      self.threshold_display = uicontrol('Parent',control_panel, ...
          'Tag','threshold','Style','edit','Units','normalized', ...
          'Position',[0.51 0.0 0.13 0.4], ...
          'String',num2str(threshold), ...
          'Callback',@(varargin) self.notify('user_event', ...
          posrecon.user_event(varargin{:})));
      self.automode_radiobox = uibuttongroup('Parent',control_panel, ...
          'Tag','automode radiobox','Units','normalized', ...
          'Position',[0.66 0.5 0.33 0.4]);
      mode_labels = { ...
          'none', ...
          'auto', ...
          'clear' };
      for i = 1:numel(mode_labels)
        self.radio_handles(i) = uicontrol('Parent',self.automode_radiobox, ...
            'Style','radiobutton','Units','normalized', ...
            'Position',[(i-1)/numel(mode_labels) 0 i/numel(mode_labels) 1], ...
            'Min',false,'Max',true,'String',mode_labels{i});
      end
      set(self.automode_radiobox,'SelectionChangeFcn',@(varargin) ...
          self.notify('user_event',posrecon.user_event(varargin{:})), ...
          'SelectedObject',self.radio_handles(1));
      self.progress_meter = uicontrol('Parent',control_panel, ...
          'Style','text','Units','normalized', ...
          'Position',[0.66 0.0 0.33 0.4],'HorizontalAlignment','right', ...
          'String',sprintf('%3.1f%% complete',progress_fraction));
      self.picture_axes = axes('Parent',self.figure,'NextPlot','add', ...
          'DrawMode','fast','Position',[0.01 0.01 0.98 0.88], ...
          'Layer','top','Color','k','Box','on','GridLineStyle','none', ...
          'XColor','w','YColor','w', ...
          'DefaultLineHitTest','off','DefaultImageHitTest','off');
      % All children drawn on picture_axes must have their 'HitTest' property
      % set to 'off' to avoid intercepting mouse clicks intended for
      % picture_axes. This is set through default properties of self.figure 
      self.picture_container = image('Parent',self.picture_axes, ...
          'XData',[1 size(picture,2)], ...
          'YData',[1 size(picture,1)]/pixel_aspect_ratio, ...
          'CData',picture,'AlphaData', ...
          intmax('uint8')*ones([size(picture,1) size(picture,2)],'uint8'));
      self.mask_container = image('Parent',self.picture_axes, ...
          'XData',[1 size(picture,2)], ...
          'YData',[1 size(picture,1)]/pixel_aspect_ratio, ...
          'CData',repmat(reshape(self.mask_color,[1 1 3]), ...
          [size(picture,1) size(picture,2) 1]), ...
          'AlphaData',self.mask_opacity*uint8(mask));
      % replace front_marker and back_marker with circles whose size is expresed
      % in the axes data units
      self.markers_container.front_marker = line('Parent',self.picture_axes, ...
          'XData',front_epsilon*self.dummy_x + front_coordinates(1), ...
          'YData',front_epsilon*self.dummy_y + front_coordinates(2), ...
          'LineStyle','-','Marker','none','Color',self.front_color, ...
          'LineWidth',1);
      self.markers_container.back_marker = line('Parent',self.picture_axes, ...
          'XData',back_epsilon*self.dummy_x + back_coordinates(1), ...
          'YData',back_epsilon*self.dummy_y + back_coordinates(2), ...
          'LineStyle','-','Marker','none','Color',self.back_color, ...
          'LineWidth',1);
      self.markers_container.line_segment = line('Parent',self.picture_axes, ...
          'XData',[front_coordinates(1) back_coordinates(1)], ...
          'YData',[front_coordinates(2) back_coordinates(2)], ...
          'LineStyle','-','Marker','none','Color',self.line_color, ...
          'LineWidth',1);
      % Scale the image axes to match the picture dimensions and configure to
      % intercept mouse clicks. Both 'Visible' and 'HitTest' properties must
      % equal 'on' for the axes to intercept mouse clicks. The 'ButtonDownFcn'
      % function is constructed in a peculiar way to pass an inline struct that
      % reports the click location and the type of click, because MATLAB does
      % not convey this information by default
      set(self.picture_axes,'Visible','on','HitTest','on', ...
          'PlotBoxAspectRatioMode','manual','TickDirMode','manual', ...
          'XTick',[],'YTick',[],'ZTick',[],'DataAspectRatio',[1 1 1], ...
          'XTickLabelMode','manual','YTickLabelMode','manual', ...
          'ZTickLabelMode','manual','ZLim',[0 1], ...
          'XLim',[1-0.5 size(picture,2)+0.5], ...
          'YLim',[1-0.5 size(picture,1)+0.5]/pixel_aspect_ratio, ...
          'ALim',uint8([0 intmax('uint8')]), ...
          'ButtonDownFcn',@(varargin) self.notify('user_event', ...
          posrecon.user_event(varargin{:})));
      % Now that we're finished populating the static gui elements, use set
      % methods (which have side effects) to achieve desired state.
      self.picture = picture;
      self.mask = mask;
      self.frame_number = frame_number;
      self.timestamp = timestamp;
      % order matters! front_epsilon and back_epsilon get assigned first so
      % that the gui knows how to draw the front and back markers
      self.front_epsilon = front_epsilon;
      self.back_epsilon = back_epsilon;
      self.front_coordinates = front_coordinates;
      self.back_coordinates = back_coordinates;
      self.threshold = threshold;
      self.automode = automode;
      self.enable_input = enable_input;
      % Reveal the gui
      set(self.figure,'Visible','on');
    end  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is just a wrapper around the line command to allow outsiders to
    % manipulate the line without exposing the handle or the parent axes
    function identifier = add_overlay(self,func,varargin)
      identifier = self.next_overlay_identifier;
      % We want the graphics to plot on self.picture_axes. Ideally we would
      % append arguments to varargin to specify self.picture_axes as the
      % Parent, but unfortunately line(...) does not take a 'Parent' argument.
      % This method will behave unpredictably if the caller includes a 'Parent'
      % argument that specifies some other axes.
      set(self.figure,'CurrentAxes',self.picture_axes);
      self.overlays(identifier) = func(varargin{:});
      self.next_overlay_identifier = identifier + 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is just a wrapper around the set method for the line overlay
    function update_overlay(self,identifier,varargin)
      set(self.overlays(identifier),varargin{:});
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is just a wrapper around the delete method for the line overlay
    function remove_overlay(self,identifier)
      delete(self.overlays(identifier));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Destructor
    function delete(self)
      delete(self.figure);
    end

  end % end public methods

  methods % set methods - be careful, no error checking!

    function set.picture(self,value)
      set(self.picture_container,'CData',value);
      self.picture = value;
    end

    function set.mask(self,value)
      set(self.mask_container,'AlphaData',self.mask_opacity*uint8(value));
      self.mask = value;
    end

    function set.frame_number(self,value)
      set(self.frame_counter,'String',num2str(value));
      self.frame_number = value;
    end

    function set.timestamp(self,value)
      set(self.clock,'String',ts2str(value));
      self.timestamp = value;
    end

    function set.front_coordinates(self,front_coords)
      set(self.markers_container.front_marker, ...
          'XData',self.front_epsilon*self.dummy_x + front_coords(1), ...
          'YData',self.front_epsilon*self.dummy_y + front_coords(2));
      set(self.markers_container.line_segment, ...
          'XData',[front_coords(1) self.back_coordinates(1)], ...
          'YData',[front_coords(2) self.back_coordinates(2)]);
      set(self.front_coordinates_display,'String', ...
          sprintf('front (%.1f,%.1f)',front_coords(1),front_coords(2)));
      self.front_coordinates = front_coords;
    end

    function set.back_coordinates(self,back_coords)
      set(self.markers_container.back_marker, ...
          'XData',self.back_epsilon*self.dummy_x + back_coords(1), ...
          'YData',self.back_epsilon*self.dummy_y + back_coords(2));
      set(self.markers_container.line_segment, ...
          'XData',[self.front_coordinates(1) back_coords(1)], ...
          'YData',[self.front_coordinates(2) back_coords(2)]);
      set(self.back_coordinates_display,'String', ...
          sprintf('back (%.1f,%.1f)',back_coords(1),back_coords(2)));
      self.back_coordinates = back_coords;
    end
    
    function set.front_epsilon(self,value)
      set(self.markers_container.front_marker, ...
          'XData',value*self.dummy_x + self.front_coordinates(1), ...
          'YData',value*self.dummy_y + self.front_coordinates(2));
      set(self.front_epsilon_display,'String',num2str(value));
      self.front_epsilon = value;
    end

    function set.back_epsilon(self,value)
      set(self.markers_container.back_marker, ...
          'XData',value*self.dummy_x + self.back_coordinates(1), ...
          'YData',value*self.dummy_y + self.back_coordinates(2));
      set(self.back_epsilon_display,'String',num2str(value));
      self.back_epsilon = value;
    end

    function set.threshold(self,value)
      set(self.threshold_display,'String',num2str(value));
      self.threshold = value;
    end

    function set.progress_fraction(self,value)
      set(self.progress_meter,'String',sprintf('%3.1f%% complete',100*value));
      self.progress_fraction = value;
    end

    function set.automode(self,value)
      h = self.radio_handles(value);
      set(h,'Value',get(h,'Max'));
      self.automode = value;
    end

    function set.enable_input(self,value)
      if value
        for i = 1:numel(self.menu_handles)
          set(self.menu_handles(i),'Enable','on');
        end
        for i = 1:numel(self.radio_handles)
          set(self.radio_handles(i),'Enable','on');
        end
        set(self.frame_counter,'Enable','on');
        set(self.front_epsilon_display,'Enable','on');
        set(self.back_epsilon_display,'Enable','on');
        set(self.threshold_display,'Enable','on');
        set(self.picture_axes,'HitTest','on');
      else
        for i = 1:numel(self.menu_handles)
          set(self.menu_handles(i),'Enable','off');
        end
        for i = 1:numel(self.radio_handles)
          set(self.radio_handles(i),'Enable','inactive');
        end
        set(self.frame_counter,'Enable','inactive');
        set(self.front_epsilon_display,'Enable','inactive');
        set(self.back_epsilon_display,'Enable','inactive');
        set(self.threshold_display,'Enable','inactive');
        set(self.picture_axes,'HitTest','off');
      end
      self.enable_input = value;
    end

  end % end set methods

end

