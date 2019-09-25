classdef video_reader < handle
% Opaque wrapper class around video_reader_mex

  properties (GetAccess = public, SetAccess = private)
    mpeg_stream_filename
    mpeg_index_filename
    timestamps_filename
    width
    height
    aspect_ratio_numerator
    aspect_ratio_denominator
    number_of_frames
    timestamps
  end

  properties (GetAccess = private, SetAccess = private)
    current_index
  end


  methods (Access = public)

    function self = video_reader(mpeg_stream_filename, ...
        mpeg_index_filename,timestamps_filename);
      if (nargin ~= 3)
        error(['constructor must be called with 3 arguments: ' ...
            'mpeg_stream_filename, mpeg_index_filename, ' ...
            'timestamps_filename']);
      end
      if (exist(mpeg_stream_filename) ~= 2)
        error('mpeg stream filename %s does not exist',mpeg_stream_filename);
      end
      if (exist(mpeg_index_filename) ~= 2)
        error('mpeg index filename %s does not exist',mpeg_index_filename);
      end
      if (exist(timestamps_filename) ~= 2)
        error('timestamps filename %s does not exist',timestamps_filename);
      end
      if mislocked('video_reader_mex')
        error('video_reader_mex MEX-file must be cleared from memory first');
      end
      try
        video_reader_mex('open',mpeg_stream_filename, ...
            mpeg_index_filename,timestamps_filename);
      catch
        error(['Could not open handle to video reader object. This class ' ...
            'depends on the private function video_reader_mex. Check ' ...
            'that this mex function is compiled and located in the ' ...
            'private folder for this class.']);
      end
      % set immutable properties
      self.mpeg_stream_filename = mpeg_stream_filename;
      self.mpeg_index_filename = mpeg_index_filename;
      self.timestamps_filename = timestamps_filename;
      self.timestamps = video_reader_mex('all_timestamps');
      self.number_of_frames = video_reader_mex('number_of_frames');
      self.width = video_reader_mex('width');
      self.height = video_reader_mex('height');
      self.aspect_ratio_numerator = ...
          video_reader_mex('aspect_ratio_numerator');
      self.aspect_ratio_denominator = ...
          video_reader_mex('aspect_ratio_denominator');
    end

    function delete(self)
      video_reader_mex('close');
      clear('video_reader_mex');
    end

    % given a scalar index in range 1:number_of_frames, return the corresponding
    % video frame as a 3-dimensional uint8 array
    function cdata = get_picture(self,index)
      if isnumeric(index) && isscalar(index) && isreal(index) && ...
          (round(index) == index) && (index >= 1) && ...
          (index <= self.number_of_frames)
        if (index ~= self.current_index)
            self.seek(index,'beg');
        end
        cdata = video_reader_mex('current_picture');
      else
        error('index argument is not valid');
      end
    end

  end

  methods (Access = private)

    function seek(self,offset,origin)
      % offset must be an integer. note that offset uses the MATLAB 1:end
      % indexing convention
      %
      %origin can be specified as string whose legal values are
      %  'beg' -1: first frame of stream
      %  'cur'  0: current frame
      %  'end'  1: end frame of stream
      %
      if ischar(origin)
        switch origin
        case 'beg'
          new_index = offset;
        case 'cur'
          new_index = self.current_index + offset;
        case 'end'
          new_index = self.number_of_frames + offset;
        otherwise
          error('%s is not a valid origin string',origin);
        end
      elseif isnumeric(origin) && isscalar(origin) && isreal(origin) && ...
          (round(origin) == origin)
        switch origin
        case -1
          new_index = offset;
        case 0
          new_index = self.current_index + offset;
        case 1
          new_index = self.number_of_frames + offset;
        otherwise
          error('%s is not a valid origin integer',num2str(origin));
        end
      else
        error(['origin must be either an integer {-1, 0, 1} or a string ' ...
            '{"beg", "cur", "end"}']);
      end
      if (new_index > self.number_of_frames) || (new_index < 1)
        error('can not seek to frame %d; out of range',new_index);
      end
      if (new_index ~= self.current_index)
        video_reader_mex('seek',new_index-1,-1);
      end
      % IMPORTANT: update current_index property
    end

  end

  % methods block reserved for dynamic property get
  methods

    function value = get.current_index(self)
      value = video_reader_mex('current_index');
    end

  end

end

