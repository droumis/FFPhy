function process_user_event(self,src,event)
% private method for mapping user input via gui to underlying events

  switch get(event.src,'Type')
  case 'axes'
    % If the source was an axes, then we are interested in mouse clicks. Some
    % trickery is required to extract the desired information about the
    % mouse click
    cp = get(event.src,'CurrentPoint');
    switch get(get(event.src,'Parent'),'SelectionType');
    case 'normal'
      self.model.assign_front(self.current_index,cp(1,1:2));
    case 'alt'
      self.model.assign_back(self.current_index,cp(1,1:2));
    case 'extend'
      % toggle playback
      if strcmp('on',self.playback_timer.Running)
        stop(self.playback_timer);
      else
        start(self.playback_timer);
      end
    end
  case 'uimenu'
    switch get(event.src,'Label')
    % these strings must exactly match the 'Label' property of the uimenu
    % objects!
    case 'Quit'
      self.quit();
    case 'Export video frame'
      self.export_picture();
    case 'Export assignments'
      self.export_assignments();
    case 'Save assignments'
      self.save_assignments();
    case 'Load assignments'
      self.load_assignments();
    case 'Summary'
      self.summarize();
    otherwise
      warning('Label property of uimenu object is not recognized');
    end
  case 'uipanel'
    if strcmp('SelectionChanged',event.data.EventName)
      switch get(event.data.NewValue,'String')
      case 'none'
        self.autoclear = false;
        self.autoassign = false;
      case 'auto'
        self.autoassign = true;
        % as a side effect, self.autoclear = false;
      case 'clear'
        self.autoclear = true;
        % as a side effect, self.autoassign = false;
      otherwise
        warning('String property of radio button is not recognized');
      end
    end
  case 'uicontrol'
    switch get(event.src,'Tag');
    case 'frame counter'
      value = str2num(get(event.src,'String'));
      if isnumeric(value) && isscalar(value) && isreal(value) && ...
          (value >= 1) && (value <= self.number_of_frames)
        % this set operation triggers side effects...
        self.current_index = self.start_index + value - 1;
      else
        set(event.src,'String', ...
            num2str(self.current_index - self.start_index + 1));
        disp('input is not a valid frame number');
      end
    case 'front epsilon'
      value = str2num(get(event.src,'String'));
      if isnumeric(value) && isscalar(value) && isreal(value) && ...
          (value > 0)
        % this set operation triggers side effects...
        self.model.epsilon = [value self.model.epsilon(2)];
      else
        set(event.src,'String',num2str(self.model.epsilon(1)));
        disp('input is not a valid proximity value');
      end
    case 'back epsilon'
      value = str2num(get(event.src,'String'));
      if isnumeric(value) && isscalar(value) && isreal(value) && ...
          (value > 0)
        % this set operation triggers side effects...
        self.model.epsilon = [self.model.epsilon(1) value];
      else
        set(event.src,'String',num2str(self.model.epsilon(2)));
        disp('input is not a valid proximity value');
      end
    case 'threshold'
      value = str2num(get(event.src,'String'));
      if isnumeric(value) && isscalar(value) && isreal(value) && ...
          (uint8(value) == value)
        % this set operation triggers side effects...
        self.luminance_threshold = uint8(value);
      else
        set(event.src,'String',num2str(self.luminance_threshold));
        disp('input is not a valid luminance threshold');
      end
    otherwise
      warning('Tag property of uicontrol object is not recognized');
    end
  case 'figure'
    if ~isempty(event.data)
      % if there is an event data struct the event was a Keystroke
      switch event.data.Key
      case 'q'
        if strcmp('control',event.data.Modifier)
          self.quit();
        end
      case 'space'
        % toggle playback
        if strcmp('on',self.playback_timer.Running)
          stop(self.playback_timer);
        else
          start(self.playback_timer);
        end
      case 'j'
        % jump to an unassigned frame
        if strcmp('off',self.playback_timer.Running)
          self.jump_to_unassigned();
        end
      case 'f'
        % play in forward direction
        self.playback_timer.TimerFcn = @(varargin) self.play_forward();
        if strcmp('off',self.playback_timer.Running)
          start(self.playback_timer);
        end
      case 'r'
        % play in backward direction
        self.playback_timer.TimerFcn = @(varargin) self.play_reverse();
        if strcmp('off',self.playback_timer.Running)
          start(self.playback_timer);
        end
      case 'n'
        self.current_index = self.current_index + 1;
      case 'b'
        self.current_index = self.current_index - 1;
      case 'a'
        % toggle autoassign mode
        self.autoassign = ~self.autoassign;
      case 'c'
        self.model.assign_front(self.current_index,[NaN NaN]);
        self.model.assign_back(self.current_index,[NaN NaN]);
      case 'add'
        % uint8 class automatically imposes range limits
        self.luminance_threshold = self.luminance_threshold + uint8(5);
      case 'equal'
        % uint8 class automatically imposes range limits
        self.luminance_threshold = self.luminance_threshold + uint8(5);
      case 'hyphen'
        % uint8 class automatically imposes range limits
        self.luminance_threshold = self.luminance_threshold - uint8(5);
      case 'subtract'
        % uint8 class automatically imposes range limits
        self.luminance_threshold = self.luminance_threshold - uint8(5);
      end
    else
      % If no event data is provided then we surmise that the figure window
      % received a Close request. This isn't an ideal approach, but I can't
      % figure out how to identify a Close request with certainty.
      self.quit();
    end
  otherwise
    warning('Type property of callback source is not recognized');
  end
end


