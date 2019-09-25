classdef assignment_event < event.EventData
% event subclass to notify of updates to the position model

  properties
    index % which video frame?
    front_coordinates = [NaN NaN];
    back_coordinates = [NaN NaN];
  end

  methods
    % event constructor
    function self = assignment_event(index,front_coords,back_coords)
      self.index = index;
      self.front_coordinates = front_coords;
      self.back_coordinates = back_coords;
    end
  end

end

