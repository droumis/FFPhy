function out = getspikewidth(index, excludeperiods, spikes, varargin)


appendindex = 0;
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

if isfield(spikes{index(1)}{index(2)}{index(3)}{index(4)},'spikewidth')
    spikewidth = spikes{index(1)}{index(2)}{index(3)}{index(4)}.spikewidth;
else 
    spikewidth = nan;
end

if (appendindex)
    out = [index spikewidth];
else
    out = spikewidth;
end