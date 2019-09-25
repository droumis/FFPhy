function [out] = DFA_getripcoactivprob_linfields(index, excludeperiods, spikes, ripples, task,linfields, cellinfo, cellfilter, varargin)
% Shantanu - Renamed from getripcoactiveprob.m
% out = getripcoactiveprob(index, excludeperiods, spikes, ripples, task, linfields,.. options)
% Gets co-activ prob during ripples and overlap for pairs of cells
% - Ripples from tetrode with max cells - Can use tetfilter. cellinfo and cellfilter not used
% - overlap from calcoverlap (prob expects calcoverlap_ver0): expects linfields data structure which has saved trajdata 
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
%   out is [coactivationprob overlap]
%       #ripples containing at leat one spike/total # ripples

% assign the options
minnumspikes = 1;
numtetrodes = 1;
minenergy = 0;
proptetrodes = [];
appendindex = 0;
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'appendindex'
            appendindex = varargin{option+1};
        case 'minnumspikes'
            nimnumspikes = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

%get tetrode with most cells
tet = gettetmaxcell(cellinfo, task, index(1,1), cellfilter, 2);

r = ripples{index(1,1)}{index(1,2)}{tet(1)};
s1 = spikes{index(1,1)}{index(1,2)}{index(1,3)}{index(1,4)};
s2 = spikes{index(1,1)}{index(1,2)}{index(1,5)}{index(1,6)};

if size(s1.data, 1) > minnumspikes  & size(s2.data, 1) > minnumspikes
    activcount = 0;
    numrip = 0;
    for i = 1:length(r.starttime) %for each ripple
        if ~isExcluded(r.midtime(i), excludeperiods) %if ripple is not excluded by exclude periods, proceed with activ count
            numrip = numrip+1;
            if sum(isExcluded(s1.data(:,1),[r.starttime(i) r.endtime(i)])) > 0 & sum( isExcluded(s2.data(:,1),[r.starttime(i) r.endtime(i)]))>0;
                activcount = activcount+1;
            end
        end
    end
    activprob = activcount/numrip;
else
    activprob = NaN; %if not enough spike
end

lf1 = linfields{index(1)}{index(2)}{index(3)}{index(4)}; %lf1 = all traj
lf2 = linfields{index(1)}{index(2)}{index(5)}{index(6)};
overlap = calcoverlap(lf1,lf2);

%calculate ripple rate
if appendindex == 0
    out = [activprob overlap];
elseif appendindex ==1
    out = [index activprob overlap];
end

end

