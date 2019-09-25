function out = calcspeedautocorrelation(index, excludetimes, pos,varargin)

for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'bin'
                bin = varargin{option+1};
            end
    else
        error('Options must be strings, followed by the variable');
    end
end

speed = pos{index(1)}{index(2)}.data(:,8);
timestep = median(diff(pos{index(1)}{index(2)}.data(:,1)));

max_lag = ceil(bin/timestep);

[c lags] = xcorr(speed,speed,max_lag,'coeff');

out.xcorr = c;
out.lags = lags;

end
