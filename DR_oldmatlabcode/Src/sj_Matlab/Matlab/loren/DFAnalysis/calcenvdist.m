function out = calcenvdist(index, excludeperiods, theta, varargin)


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

% assign a temporary variable for the theta
t = theta{index(1)}{index(2)}{index(3)};
clear theta;

% the envelope is the third field of the data element
env = t.data(:,3);

% create a list of times for the theta data
tme = t.starttime + [0:(length(env)-1)] / t.samprate;

% apply the exclude filter
env = env(find(~isExcluded(tme, excludeperiods)));

out.bins = 0:round(max(env)/100):(max(env)+1);
out.envhist = histc(env, out.bins);

if (appendindex)
    out.index = index;
end
