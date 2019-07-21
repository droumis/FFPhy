

function out = trim2win(data, srate, pwin, varargin)
% data dim 2 must be time
dsamp = 1; % 1 = not downsampled
if ~isempty(varargin)
   assign(varargin{:}); 
end

midpoint = floor(size(data,2)/2); %get middle index of window
sampsrate = (srate/floor(dsamp));
out = data(:,midpoint-pwin(1)*sampsrate:midpoint+pwin(2)*sampsrate,:);
end