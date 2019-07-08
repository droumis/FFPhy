

function out = trim2win(inData, srate, pwin, varargin)

dsamp = 1; % 1 = not downsampled
if ~isempty(varargin)
   assign(varargin{:}); 
end

midpoint = floor(size(inData,2)/2); %get middle index of window
sampsrate = (srate/floor(dsamp));
out = inData(:,midpoint-pwin(1)*sampsrate:midpoint+pwin(2)*sampsrate);
end