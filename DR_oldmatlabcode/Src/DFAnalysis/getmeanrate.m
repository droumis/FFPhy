function out = getmeanrate(index, excludeperiods, spikes, varargin)


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


try
    epochtime = diff(spikes{index(1)}{index(2)}{index(3)}{index(4)}.timerange)/10000;
    numspikes = size(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data,1);
    meanrate = numspikes/epochtime;
catch
    meanrate = nan;
end

if (appendindex)
    out = [index meanrate];
else
    out = meanrate;
end