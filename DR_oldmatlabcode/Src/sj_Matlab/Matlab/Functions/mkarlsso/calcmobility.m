function out = calcmobility(index, excludetimes, pos, varargin)
% out = calcmobility(index, excludetimes, options)
% Calculates the proportion of time spent in the designated mobility state
%
% Options:
%   'mintime' - minimum immobility time
%   'maxtime' - maximum immobility time
%   'appendindex' - whether or not to append the epoch index to the output
%



mintime = -1;
maxtime = inf;
appendindex = 0;

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'mintime'
                mintime = varargin{option+1};
            case 'maxtime'
                maxtime = varargin{option+1};  
            case 'appendindex'
                appendindex = varargin{option+1};  
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

if ~isempty(excludetimes)
    goodtimeind = find(~isExcluded(pos{index(1)}{index(2)}.data(:,1), excludetimes));
else
    goodtimeind = [1:length(pos{index(1)}{index(2)}.data(:,1))]';
end
tmpvelocity = pos{index(1)}{index(2)}.data(:,5);
timestep = (pos{index(1)}{index(2)}.data(2,1) - pos{index(1)}{index(2)}.data(1,1));

%immobilitytime = cumsum_reset(tmpvelocity<(mean(tmpvelocity)*.05))*timestep;
immobilitytime = cumsum_reset(tmpvelocity<(2))*timestep;

for i = 1:length(maxtime)
    %out(1,i) = sum((immobilitytime(goodtimeind) > mintime(i)) & (immobilitytime(goodtimeind) < maxtime(i)))/length(goodtimeind);
    out(1,i) = sum((immobilitytime(goodtimeind) > mintime(i)) & (immobilitytime(goodtimeind) < maxtime(i))) * timestep;
end
    
if (appendindex)
    out = [index out]; %append the cell index to the value    
end
%--------------------------------------------------------------------------

function out = cumsum_reset(inputvect)

resetpoints = find(diff(inputvect) < 0)+1;
out = cumsum(inputvect);
for i = 1:length(resetpoints)
    out(resetpoints(i):end) = out(resetpoints(i):end)-out(resetpoints(i)-1);
end




