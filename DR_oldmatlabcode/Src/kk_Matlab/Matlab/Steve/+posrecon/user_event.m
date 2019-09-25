classdef user_event < event.EventData
% event subclass for user interaction with gui

% properties inherited from event.EventData superclass:
%   EventName
%   Source

  properties
    % user_event.src is the local source of the keypress/mouseclick/callback
    % that triggered the event notification (for example, a gui widget).
    % This is not always the same thing as user_event.Source, which is the
    % object whose notify method was invoked.
    src
    % data contains any event-specific data that was passed in
    data
    % extra_args contain any additional data after the data argument; this is
    % usually an empty cell array
    extra_args
  end

  methods
    function self = user_event(src,data,varargin);
      self.src = src;
      self.data = data;
      self.extra_args = varargin;
    end
  end

end

