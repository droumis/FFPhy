function overlap = calc2doverlap(rate1,rate2, varargin)
%overlap = calc2doverlap(rate1,rate2, varargin)
%
% rate1 and rate2 are matrices of occupancy normalized firing rates. NaN
% are okay and will be excluded from analysis.  If more than half bins are
% excluded, overlap will not be computed
%
% options
%   Normalize, 0 or 1, default 0
%       if 0 calculates overlap, if 1 calculates normalized overlap
%   MinBins, value between 0 and 1, default 0.5
%       proportion of bins that must be defined to calculate overlap,
%       otherwise overlap = NaN
%
normalize = 0;
minbins = 0.25;
binsize = 1;
for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'normalize'
                normalize = varargin{option+1};
            case 'minbins'
                minbins = varargin{option+1};
            case 'binsize'
                binsize = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

overlap = [];
if ~isempty(rate1) && ~isempty(rate2)
    %First figure out area of intersection
    area12 = min(rate1,rate2);
    area12(isnan(rate1)|isnan(rate2)) = NaN;
    area12 = area12.*binsize.*binsize;
    area12 = nansum(nansum(area12));
    
    %Now figure out area of each place field
    area1 = rate1.*binsize.*binsize;
    area1(isnan(rate1)|isnan(rate2)) = NaN;
    area1 = nansum(nansum(area1));
    area2 = rate2.*binsize.*binsize;
    area2(isnan(rate1)|isnan(rate2)) = NaN;
    area2 = nansum(nansum(area2));
    if sum(sum(isnan(rate1+rate2))) < minbins * size(rate1,1) * size(rate2,2)
        overlap = -1; %low occupancy
    else
        %overlap is 2 times the intersection / by the sum of the areas        
        overlap = 2*area12/(area1+area2);
    end

end