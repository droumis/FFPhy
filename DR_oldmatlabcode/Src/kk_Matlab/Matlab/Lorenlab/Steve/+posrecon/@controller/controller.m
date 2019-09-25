classdef controller < handle
% The controller is responsible for listening to all event notifications and
% invoking the appropriate methods of data, model, and gui.

  properties (GetAccess = public, SetAccess = private)
    % interface to the video data, of class video_reader
    data
    % maintains the estimate of the rat's position and has access to data
    model
    % thin stateless client that displays to the user and receives user input;
    % emits event notifications that are processed by controller, which in turn
    % invokes gui methods
    gui

    % timer object for scheduling playback of video in the gui
    playback_timer;
    playback_interval = 0.001; % some small positive number

    % these are self-explanatory
    mpeg_stream_filename;
    mpeg_index_filename;
    timestamps_filename;
    % range of timestamps of interest
    timerange
    % indices into the video frames
    start_index
    end_index
    current_index
    % number_of_frames == numel(start_index:end_index)
    number_of_frames
    % threshold for detecting bright pixels
    luminance_threshold = intmax('uint8');
    % baseline for background subtraction
    background_luminance
    % flag for whether the position model should attempt to make automatic
    % position assignments
    autoassign = false;
    % flag for whether the controller should automatically clear assignments for
    % all video frames that appear during playback
    autoclear = false;
    % autoassign and autoclear are mutually exclusive; the controller must never
    % be in a state such that (autoassign & autoclear) == true
  end

  properties (GetAccess = public, SetAccess = private, SetObservable = true)
    dirty = true
  end

  methods (Access = public)

    % Constructor
    function self = controller(mpeg_stream_filename,mpeg_index_filename, ...
        timestamps_filename,timerange)
      if (nargin ~= 4)
        error(['constructor must be called with 4 arguments: ' ...
            'mpeg_stream_filename, mpeg_index_filename, ' ...
            'timestamps_filename, timerange']);
      end
      if (exist('video_reader') ~= 2)
        error('Can not find the VIDEO_READER class (written by SMK)');
      end
      if (exist('ts2str') ~= 2)
        error('Can not find the m-file TS2STR (written by SMK)');
      end
      if isa(timerange,'uint32') && isvector(timerange) && ...
          (size(timerange,1) == 1) && (size(timerange,2) == 2)
        self.timerange = timerange;
      else
        error(['timerange argument must be a 2-element row vector of ' ...
            'uint32 timestamps in increasing order']);
      end
      self.mpeg_stream_filename = mpeg_stream_filename;
      self.mpeg_index_filename = mpeg_index_filename;
      self.timestamps_filename = timestamps_filename;
      % attach interface to video data
      self.data = video_reader(mpeg_stream_filename,mpeg_index_filename, ...
          timestamps_filename);
      self.start_index = find(self.data.timestamps < self.timerange(1), ...
          1,'last');
      if isempty(self.start_index)
        error('video timestamps do not cover the start time %s', ...
            ts2str(self.timerange(1)));
      end
      self.end_index = find(self.data.timestamps > self.timerange(2), ...
          1,'first');
      if isempty(self.end_index)
        error('video timestamps do not cover the end time %s', ...
            ts2str(self.timerange(2)));
      end
      self.number_of_frames = self.end_index - self.start_index + 1;
      % compute background image for noise suppression
      disp('computing background luminance; this may take a long while');
      self.compute_background_luminance(1000);
      % attach gui
      picture = self.data.get_picture(self.start_index);
      mask = picture(:,:,1) - self.background_luminance > ...
          self.luminance_threshold;
      self.gui = posrecon.gui(picture,mask, ...
          self.data.aspect_ratio_numerator/ ...
          self.data.aspect_ratio_denominator, ...
          1,self.data.timestamps(self.start_index), ...
          [NaN NaN],[NaN NaN],1,1, ...
          self.luminance_threshold,0.0,1,true);
      % attach position model, which reads data and draws its own stuff on gui
      self.model = posrecon.model(self.data,self.gui);
      % timer for interruptible playback. 'ExecutionMode' property must be set to
      % 'fixedSpacing' to avoid lockups. The timer is linked to the gui to
      % prevent the user for manually seeking to a frame during playback
      self.playback_timer = timer('TimerFcn', ...
          @(varargin) self.play_forward(), ...
          'Period',self.playback_interval,'ExecutionMode','fixedSpacing', ...
          'StartFcn',@(varargin) set(self.gui,'enable_input',false), ...
          'StopFcn',@(varargin) set(self.gui,'enable_input',true));
      set(self.playback_timer,'BusyMode','error', ...
          'ErrorFcn',@(varargin) stop(self.playback_timer));
      % set current gui state
      self.current_index = self.start_index;
      % listen to the gui and model
      addlistener(self.model,'assignment_event', ...
          @(src,event) self.process_assignment_event(src,event));
      addlistener(self.gui,'user_event', ...
          @(src,event) self.process_user_event(src,event));
    end  

    % this computes a background image for noise suppression
    compute_background_luminance(self,n)

    % this handles input from the gui (external m-file)
    process_user_event(self,src,event)
  
    % this responds whenever the model makes an assignment
    function process_assignment_event(self,src,event)
      self.dirty = true;
      % we only need to update gui.front_coordinates and gui.back_coordinates
      % if it is currently displaying the video frame for which the assignment
      % was made
      match = find(event.index == self.current_index);
      if ~isempty(match)
        self.gui.front_coordinates = event.front_coordinates(match,:);
        self.gui.back_coordinates = event.back_coordinates(match,:);
      end
      % update gui.progress_fraction
      self.gui.progress_fraction = ...
          nnz(self.model.assigned(self.start_index : self.end_index)) / ...
          self.number_of_frames;
    end

    function delete(self)
      if isobject(self.gui)
        delete(self.gui);
      end
      if isobject(self.model)
        delete(self.model);
      end
      if isobject(self.data)
        delete(self.data);
      end
      delete(self.playback_timer);
      % After this point, MATLAB automatically calls the delete method of the
      % video_reader superclass
    end

    % export currently displayed picture to workspace
    function h = export_picture(self)
      cdata = self.data.get_picture(self.current_index);
      h = export2wsdlg({'Export current image to workspace as:'},{'cdata'}, ...
          {cdata},'Export current image to workspace',[true]);
      ch = get(h,'Children');
      % export2wsdlg appends a trailing number to the specified variable name if
      % a variable with that name already exists. Replace with the original
      % variable name to allow overwrite (user will be prompted).
      set(ch(3),'String','cdata');
      waitfor(h);
    end

    % export current position assignments as a struct
    function h = export_assignments(self)
      % compress from double- to single-precision
      index_range = self.start_index : self.end_index;
      rawpos = struct( ...
          'timestamp',{self.data.timestamps(index_range)}, ...
          'xfront',{single(self.model.xfront(index_range))}, ...
          'yfront',{single(self.model.yfront(index_range))}, ...
          'xback',{single(self.model.xback(index_range))}, ...
          'yback',{single(self.model.yback(index_range))}, ...
          'sources',{{self.mpeg_stream_filename, ...
          self.mpeg_index_filename, ...
          self.timestamps_filename}}, ...
          'units',{'cm'});
      h = export2wsdlg( ...
          {'Export current assignments to workspace as:'},{'rawpos'}, ...
          {rawpos},'Export position assignments to workspace',[true]);
      ch = get(h,'Children');
      % export2wsdlg appends a trailing number to the specified variable name if
      % a variable with that name already exists. Replace with the original
      % variable name to allow overwrite (user will be prompted).
      set(ch(3),'String','rawpos');
      waitfor(h);
    end

    function save_assignments(self)
      [filename, pathname] = uiputfile('*.mat', ...
          'Save position assignments to file');
      index_range = self.start_index : self.end_index;
      rawpos = struct( ...
          'timestamp',{self.data.timestamps(index_range)}, ...
          'xfront',{single(self.model.xfront(index_range))}, ...
          'yfront',{single(self.model.yfront(index_range))}, ...
          'xback',{single(self.model.xback(index_range))}, ...
          'yback',{single(self.model.yback(index_range))});
      if isequal(filename,0)
        return;
      end
      try
        save([pathname '/' filename],'rawpos');
        self.dirty = false;
      catch
        error('could not save assignments to file %s/%s', ...
            pathname,filename);
      end
    end

    function load_assignments(self)
      [filename, pathname] = uigetfile('*.mat');
      flag = false;
      if isequal(filename,0)
        return;
      end
      try
        load([pathname '/' filename]);
      catch
        error('could not load assignments from file %s/%s', ...
            pathname,filename);
        return;
      end
      index_range = self.start_index : self.end_index;
      if ~isstruct(rawpos) || ~isscalar(rawpos) || ~all(isfield(rawpos, ...
          {'timestamp','xfront','yfront','xback','yback'})) || ...
          ~isequal(rawpos.timestamp, ...
          self.data.timestamps(index_range)) || ...
          ~all(cellfun(@(c) isfloat(c) && isvector(c) && ...
          isreal(c) && ~any(isinf(c)) && numel(c) == self.number_of_frames, ...
          {rawpos.xfront, rawpos.yfront, rawpos.xback, rawpos.yback}))
        error('data loaded from file is not valid');
      else
        % expand vectors from single to double to take advantage of MATLAB's
        % acceleration with double values.
        self.model.assign_front(index_range, ...
            [double(rawpos.xfront) double(rawpos.yfront)]);
        self.model.assign_back(index_range, ...
            [double(rawpos.xback) double(rawpos.yback)]);
      end
    end

    function summarize(self)
      index_range = self.start_index : self.end_index;
      disp(sprintf('%f seconds (%s to %s)',double(diff(self.timerange))/1e4, ...
          ts2str(self.timerange(1)),ts2str(self.timerange(2))));
      disp(sprintf('%d frames with mean inter-frame interval %f ms', ...
          self.number_of_frames, ...
          mean(diff(double(self.data.timestamps(index_range))))/10));
      n = nnz(self.model.assigned(index_range));
      if (n > 0)
        disp(sprintf('%d frames assigned (%f%% complete)',n, ...
            100*n/self.number_of_frames));
        gaps = diff(self.data.timestamps(self.model.assigned));
        first_assigned = find(self.model.assigned,1,'first');
        last_assigned = find(self.model.assigned,1,'last');
        if (first_assigned > self.start_index)
          disp('video frame preceding the start time has no assignment');
          gaps = [ ...
              self.data.timestamps(first_assigned) - self.timerange(1); ...
              gaps];
        end
        if (last_assigned < self.end_index)
          disp('video frame following the end time has no assignment');
          gaps = [gaps; ...
              self.timerange(2) - self.data.timestamps(last_assigned)];
        end
        disp(sprintf('largest gap in coverage is %f ms',max(double(gaps)/10)));
      else
        disp('no position assignments have been made');
      end

    end

    function quit(self)
      % dialog box to ask user to confirm
      response = questdlg('Do you really want to quit?','Quit', ...
          'Yes','No','No');
      if strcmp('Yes',response)
        if self.dirty
          %self.save_assignments();
          waitfor(self.export_assignments());
        end
        % call destructor
        delete(self);
      end
    end

    function play_forward(self)
      target = self.current_index + 1;
      % if we are using the autoassign feature, then skip ahead to the first
      % frame that does not have an assignment already
      if self.autoassign && (target <= self.end_index)
        target = self.current_index + find( ...
            ~self.model.assigned(target:self.end_index),1,'first');
        if isempty(target)
          stop(self.playback_timer);
          target = self.end_index;
        end
      end
      self.current_index = target;
    end

    function play_reverse(self)
      target = self.current_index - 1;
      % if we are using the autoassign feature, then skip back to the last
      % frame that does not have an assignment already
      if self.autoassign && (target >= self.start_index)
        target = self.start_index - 1 + find( ...
            ~self.model.assigned(self.start_index:target),1,'last');
        if isempty(target)
          stop(self.playback_timer);
          target = self.start_index;
        end
      end
      self.current_index = target;
    end

  end % end public methods


  methods % reserved block for get/set methods

    function self = set.current_index(self,value)
      if (value < self.start_index) || (value > self.end_index)
        disp('Can not go further');
        % stop the playback_timer if it is running; this also has the side
        % effect of changing the state of the gui through callbacks of the
        % playback_timer
        stop(self.playback_timer);
        return;
      end
      % self.current_index should be updated first, because other components may
      % wish to reference its current value
      self.current_index = value;
      picture = self.data.get_picture(value);
      mask = picture(:,:,1) - self.background_luminance > self.luminance_threshold;
      self.gui.mask = mask;
      self.gui.frame_number = value - self.start_index + 1;
      self.gui.picture = picture;
      self.gui.timestamp = self.data.timestamps(value);
      self.gui.front_coordinates = ...
          [self.model.xfront(value) self.model.yfront(value)];
      self.gui.back_coordinates = ...
          [self.model.xback(value) self.model.yback(value)];
      if self.autoassign && ~self.model.auto_assign(value,mask)
        % note that the model will emit an assignment_event notification if the
        % auto-assignment is successful, which triggers update of the gui
        stop(self.playback_timer);
      elseif self.autoclear
        % erase the assignments if we are in autoclear mode
        self.model.assign_front(value,[NaN NaN]);
        self.model.assign_back(value,[NaN NaN]);
      end
    end

    function self = set.luminance_threshold(self,value)
      % force update of the gui mask
      self.gui.mask = ...
          self.gui.picture(:,:,1) - self.background_luminance > value;
      self.gui.threshold = value;
      self.luminance_threshold = value;
    end

    function self = set.autoassign(self,value)
      if value
        self.autoclear = false;
        self.gui.automode = 2;
      end
      self.autoassign = value;
      self.gui.automode = 1 + 1*self.autoassign + 2*self.autoclear;
    end

    function self = set.autoclear(self,value)
      if value
        self.autoassign = false;
      end
      self.autoclear = value;
      self.gui.automode = 1 + 1*self.autoassign + 2*self.autoclear;
    end

  end

end

