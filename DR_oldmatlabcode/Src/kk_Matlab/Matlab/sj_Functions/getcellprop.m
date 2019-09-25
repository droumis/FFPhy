function out = cellprop(ind, excludeperiods, spikes, varargin)
% function out = cellprop(ind, excludeperiods, spikes, cellinfo, varargin)
% Returns various spikes waveform and spike rate characterstics
% 
% the output has the following fields:
% mean_waveform_width  mean_rate  csi  prop_bursts theta_mod_depth  
%   peak_theta_phase gamma_mod_depth  peak_gamma_phase
% 
%
% Options:
%   'burstlen', burstlength    specify the burst length to use in ms 
%                              (default 10 ms)
%   'appendind', 1 or 0 -- set to 1 to append the cell ind to the
%   output [tetrode cell value].  Default 0.
%

appendind = 0;
burstlength = 10;
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendind'
                appendind = varargin{option+1};
            case 'burstlength'
                burstlength = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

% create a filtered set of spikes and get their amplitudes
s = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)};
if (~isempty(excludeperiods))
    goodspikes = ~isExcluded(s.data(:,1), excludeperiods);
else
    goodspikes = 1:length(s.data(:,1));
end

width = s.spikewidth;

if ~isempty(excludeperiods)
    totalexclude = sum(excludeperiods(:,2) - excludeperiods(:,1));
else
    totalexclude = 0;
end
ontime = diff(spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.timerange)/10000;
totalontime = ontime-totalexclude;

meanrate = length(goodspikes) / ontime;

[csi propb] = computecsi(s.data(goodspikes,1), s.data(goodspikes,6), ...
			burstlength);

[thetamod peakthetaphase] = modulation(s.data(:,5));
[gammamod peakgammaphase] = modulation(s.gammaphase);


% mean_waveform_width  mean_rate  csi  prop_bursts theta_mod_depth  
%   peak_theta_phase gamma_mod_depth  peak_gamma_phase

if (appendind)
    out = [ind width meanrate csi propb thetamod peakthetaphase gammamod peakgammaphase];
else
    out = [width meanrate csi propb thetamod peakthetaphase gammamod peakgammaphase];
end
