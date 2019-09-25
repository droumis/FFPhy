function [out] = getripactivprob(index, excludeperiods, spikes, ripples, task, cellinfo, cellfilter, varargin)
% out = getripactivprob(index, excludeperiods, spikes, ripples, task, cellinfo, cellfilter, options)
%
%   index [day epoch tetrode cell] for each cell
%   tetlist save in animdirectory: ie bartetlist02, tetlist{2}{4} =
%   tetrodes to use for analysis
%
%   options are
%     'minnumspikes', default 1
%           min number of spikes in the entire epoch (within and outside ripples)
%           to include in analysis
%	  'appendindex' , 1 or 0, default 0
%           set to 1 to append the cell index to the output [day epoch
%           value]
%
%   out is [activationprob]
%       #ripples containing at leat one spike/total # ripples

% assign the options
minnumspikes = 1;
numtetrodes = 1;
minenergy = 0;
proptetrodes = [];
appendindex = 0;
appendnumrip = 0;
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'appendindex'
            appendindex = varargin{option+1};
        case 'minnumspikes'
            minnumspikes = varargin{option+1};
        case 'appendnumrip'
            appendnumrip = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

%get tetrode with most cells
tet = gettetmaxcell(cellinfo, task, index(1,1), cellfilter, 2);

r = ripples{index(1,1)}{index(1,2)}{tet(1)};
s = spikes{index(1,1)}{index(1,2)}{index(1,3)}{index(1,4)};

if size(s.data, 1) > minnumspikes
    activcount = 0;
    numrip = 0;
    for i = 1:length(r.starttime) %for each ripple
        if ~isExcluded(r.starttime(i), excludeperiods) %if ripple is not excluded by exclude periods, proceed with activ coundt
            numrip = numrip + 1;
            if  sum(isExcluded(s.data(:,1),[r.starttime(i) r.endtime(i)]))>0
                activcount = activcount+1;
            end
        end
    end
    if numrip>0
        activprob = activcount/numrip;
    else
        activprob = NaN;
    end
else
    activprob = NaN; %if not enough spike
end
%calculate ripple rate
if appendindex == 0
    out = activprob;
elseif appendindex ==1
    out = [index activprob];
end

if appendnumrip == 1
    out = [out numrip];
end

end

