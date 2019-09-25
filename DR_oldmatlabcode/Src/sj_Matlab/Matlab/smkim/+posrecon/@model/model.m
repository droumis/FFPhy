classdef model < handle

  properties (GetAccess = public, SetAccess = public, SetObservable = true)
    % a 2-element vector of the scale (in pixels) of the front and back markers
    epsilon 
    % sum of the elements of epsilon
    sum_epsilon
    % number of samples to look forward and back for estimation
    halfwidth
  end

  properties (GetAccess = public, SetAccess = private, SetObservable = true)
    % current assignments
    xfront
    yfront
    xback
    yback
    % logical flag for whether a frame has complete (non-NaN) assignments
    assigned
  end
 
  properties (GetAccess = private, SetAccess = private)
    data
    gui
    % coordinates of the pixel centers, adjusted for pixel aspect ratio
    xgrid
    ygrid
    % complex-valued representation of the line segment between the front and
    % back markers. abs(z) is the distance between the front and back, and
    % angle(z) is the direction
    z
  end 

  events
    assignment_event
  end

  methods (Access = public)

    function self = model(data,gui)
      if (exist('kmeans') ~= 2)
        error('posrecon.model depends on the m-file KMEANS (MATLAB Statistics Toolbox');
      end
      if isa(data,'video_reader')
        self.data = data;
      else
        error('DATA argument must be a video_reader object');
      end
      if isa(gui,'posrecon.gui')
        self.gui = gui;
      else
        error('GUI argument must be a posrecon.gui object');
      end
      n = data.number_of_frames;
      self.xfront = nan([n 1]);
      self.yfront = nan([n 1]);
      self.xback = nan([n 1]);
      self.yback = nan([n 1]);
      self.z = complex(self.xfront-self.xback,self.yfront-self.yback);
      self.assigned = false([n 1]);
      [self.xgrid, self.ygrid] = meshgrid((1:data.width), ...
          data.aspect_ratio_denominator/ ...
          data.aspect_ratio_numerator*(1:data.height));
      % silence warnings related to kmeans
      warning('off','stats:kmeans:FailedToConverge');
      warning('off','stats:kmeans:EmptyCluster'); 
      self.epsilon = [5 4];
      self.halfwidth = 3;
    end
     
    function assign_front(self,index,front_coords)
      % index must a vector and front_coords must be a numel(index)x2 array
      self.xfront(index) = front_coords(:,1);
      self.yfront(index) = front_coords(:,2);
      % preserve back assignments
      back_coords = [self.xback(index) self.yback(index)];
      self.z(index) = complex(front_coords(:,1) - back_coords(:,1), ...
        front_coords(:,2) - back_coords(:,2));
      self.assigned(index) = ~isnan( ...
          prod(front_coords,2) .* prod(back_coords,2));
      notify(self,'assignment_event', ...
          posrecon.assignment_event(index,front_coords,back_coords));
    end

    function assign_back(self,index,back_coords)
      % index must a vector and back_coords must be a numel(index)x2 array
      self.xback(index) = back_coords(:,1);
      self.yback(index) = back_coords(:,2);
      % preserve front assignments
      front_coords = [self.xfront(index) self.yfront(index)];
      self.z(index) = complex(front_coords(:,1) - back_coords(:,1), ...
          front_coords(:,2) - back_coords(:,2));
      self.assigned(index) = ~isnan( ...
          prod(back_coords,2) .* prod(front_coords,2));
      notify(self,'assignment_event', ...
          posrecon.assignment_event(index,front_coords,back_coords));
    end

    % Method defined in a separate m-file
    value = auto_assign(self,index,mask)
    
  end

  methods % reserved block for get/set methods

    function self = set.epsilon(self,value)
      self.gui.front_epsilon = value(1);
      self.gui.back_epsilon = value(2);
      self.epsilon = value;
      self.sum_epsilon = sum(value);
    end

  end

end

