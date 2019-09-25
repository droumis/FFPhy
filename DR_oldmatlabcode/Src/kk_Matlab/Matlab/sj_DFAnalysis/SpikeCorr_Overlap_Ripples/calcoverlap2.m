function overlap = calcoverlap2(traj1,traj2, varargin)
% overlap = calcoverlap2(traj1,traj2, varargin)
% Expects single trajectory as input unlike calcoverlap2
% TRAJ1 and TRAJ2 are vectors of occupancy normalized firing rates, column
% 5 of trajdata{d}{e}{t}{c}{traj}, NaN are okay and will be excluded from
% analysis.  If more than half bins are exclude, overlap will not be
% computed
%
% options
%   Normalize, 0 or 1, default 0
%       if 0 calculates overlap, if 1 calculates normalized overlap
%   MinBins, value between 0 and 1, default 0.5
%       proportion of bins that must be defined to calculate overlap,
%       otherwise overlap = NaN
%
normalize = 0;
thresh = 0;
minbins = 0.5;
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'normalize'
                normalize = varargin{option+1};
            case 'minbins'
                minbins = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

overlap = [];
if ~isempty(traj1) & ~isempty(traj2)
    trajlengths = [length(traj1) length(traj2)];
    if diff(trajlengths)>2
        warning('Trajectory lengths differ by more than 2 bins')
        overlap = -5;
    else
        trajlength = min(trajlengths);
        if normalize==1
            area1 = nansum(traj1(1:trajlength));
            area2 = nansum(traj2(1:trajlength));
            if area1 == 0
                area1 = 1;
            end
            if area2 == 0
                area2 = 1;
            end
            ratediff = abs((traj1(1:trajlength)/area1) - (traj2(1:trajlength)/area2));
            totalrates = ((traj1(1:trajlength)/area1) + (traj2(1:trajlength)/area2));
        elseif normalize ==0
            ratediff = abs((traj1(1:trajlength)) - (traj2(1:trajlength)));
            totalrates = ((traj1(1:trajlength)) + (traj2(1:trajlength)));
        end
        if sum(isnan(totalrates)) > minbins * length(totalrates)
            overlap = -2; %low occ
            %    warning('half of bins excluded because of low occupancy, overlap is NaN')
        else
            ratediff = ratediff(~isnan(ratediff)); %exclude NaN bins
            totalrates = totalrates(~isnan(totalrates));
            if (sum(totalrates) > 0)
                overlap = (sum(totalrates)-sum(ratediff))/(sum(totalrates)); %overlap;
            else
                overlap = -4; %no firing in either traj
            end
        end
    end
else
    overlap = -3;
    % warning('Empty trajectory, overlap is NaN')
end

if isnan(overlap)
    keyboard
end

if isempty(overlap) %half bins empty or no firing in either traj
    overlap = NaN;
end