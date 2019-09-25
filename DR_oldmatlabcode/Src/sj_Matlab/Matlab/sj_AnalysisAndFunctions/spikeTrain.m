function hh = spikeTrain(tSpikes, baseline, height, varargin)
%function hh = spikeTrain(tSpikes, baseline, height, varargin)
hh = [];
for i = 1:length(tSpikes)
   hh(i) = line([tSpikes(i) tSpikes(i)],[baseline (baseline+height)], ...
      varargin{:});
end
