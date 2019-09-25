
function out = DFAsj_calcoverlap(index, excludetimes, linfields, mapfields, varargin)

% Takes in index for pairs of cells and calculates overlap
% No excludetimes in here - linfields has no time information

%function overlap = calcoverlap(trajdata1,trajdata2, varargin)
% overlap = calcoverlap(trajdata1,trajdata2, varargin)
% Trajdata had all trajectories unlike calcoverlap2 which expects only one trajectory as input
%
% trajdata1{traj}
% Compute overlap for all traj in trajdata
% If more than half bins are exclude, overlap will not be .computed
%
% options
%   Normalize, 0 or 1, default 0
%       if 0 calculates overlap, if 1 calculates normalized overlap
%   thresh, minimim peak to compute overlap
%   MinBins, value between 0 and 1, default 0.5
%       proportion of bins that must be defined to calculate overlap,
%       otherwise overlap = NaN

normalize = 0;
thresh = 3;
minbins = 0.5;

if ~isempty(excludetimes)
    excludetimes = [];
end

for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'normalize'
                normalize = varargin{option+1};
            case 'thresh'
                thresh = varargin{option+1};
            case 'minbins'
                minbins = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    end
end


lf1 = linfields{index(1)}{index(2)}{index(3)}{index(4)}; %lf1 = all traj
lf2 = linfields{index(1)}{index(2)}{index(5)}{index(6)};
[overlap, peakcomb, trajpeakcomb, peak1, peak2, trajpeak1, trajpeak2] = sj_calcoverlap(lf1,lf2,...
    'normalize',normalize,'thresh',thresh,'minbins',minbins);


% if (appendindex) 
%     out.index = index;
% end

out.index = index;

out.overlap = overlap; 
out.peakcomb = peakcomb;
out.trajpeakcomb = trajpeakcomb;

out.peak1 = peak1;
out.peak2 = peak2;
out.trajpeak1 = trajpeak1;
out.trajpeak2 = trajpeak2;

out.trajdata1 = lf1;
out.trajdata2 = lf2;
out.mapdata1 = mapfields{index(1)}{index(2)}{index(3)}{index(4)};
out.mapdata2 = mapfields{index(1)}{index(2)}{index(5)}{index(6)};