function [out] = DFAsj_getripparams(index, excludeperiods, ripples, varargin)
% out = getriptimes(index, excludeperiods, ripples, options)
%
% Shantanu: Sep2011. 
% To return ripple detection and size parameters
%
%   index [day epoch]
%   tetlist is a list of tetrodes to include in the analysis
%
%   options are
%	'minenergy', E
%		     specifies the minimum energy of a valid ripple event
%   'minthresh', minthresh
%		     specifies a minimum threshold in stdev units for a valid
%			ripple event  (default 0)
%   'numtetrodes'
%           specifies number of tetrodes a ripple must be recorded on to be
%           included in analysis, default 1
%   'proptetrodes', examples: 1, 0.5, 0.25,
%           proportion of tetrodes a ripple must be recorded on to be
%           included in analysis
%   'appendindex' , 1 or 0, default 0
%           set to 1 to append the cell index to the output [day epoch
%           value]
%
%   out is [rate proportiontime]
%       proprotiontime is proportion of included time during which ripples
%   were recorded
%       rate is number ripples/sec during included time

% assign the options
appendindex = 0;
for option = 1:2:length(varargin)-1
    switch varargin{option}  
        case 'appendindex'
            appendindex = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

day = unique(index(:,1));
epoch = unique(index(:,2));

tetlist = index(:,3);
r = ripples{index(1,1)}{index(1,2)}{tetlist(1)}; % For time range for current epoch

for t = 1:length(tetlist)
    tmprip = ripples{index(1,1)}{index(1,2)}{tetlist(t)}; 
    % mean, std and threshold of ripple envelope used for each tet in current epoch
    % Return a mean acroos all tets for comparing acroos tets later
    % Also, sj_calcriprate can figure out which tet is assigning ripsize,
    % and return its baseline, std and threshold. % Not implented now
    rip.baseline(t) = tmprip.baseline;
    rip.std(t) = tmprip.std;
    rip.thresh(t) = tmprip.threshold;
end

% Return mean, std and threshold of ripple envelope used for
% current epoch
out.rip_baseline = mean(rip.baseline);
out.rip_std = mean(rip.std);
out.rip_thresh = mean(rip.thresh);
if appendindex == 1
    out.index = index;
end


