function out = calcgammatail(index, excludetimes, lowgamma, ripple, ripples,cellinfo,varargin)

%Define sampling rate in Hz
minthresh = 3;
out = [];
%Set options
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'minthresh'
                mintrhesh = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

if length(index)==3
    index = [index index(3)];
end

%Find valid ripples
riptimes = getripples([index(1) index(2)], ripples, cellinfo, 'cellfilter', ...
    '(isequal($area, ''CA1''))','excludeperiods', excludetimes,'minstd',minthresh);

%Exclude ripples occuring close together for preceding modulationg
valid_ripples = [1000; riptimes(2:end,2)-riptimes(1:end-1,1)];
valid_ripples = valid_ripples > 1;
riptimes = riptimes(valid_ripples,:);

if ~isempty(riptimes)
    %Initialize lowgamma and ripples
    l = double(lowgamma{index(1)}{index(2)}{index(3)}.data(:,3));
    times = geteegtimes(lowgamma{index(1)}{index(2)}{index(3)});
    r = double(ripple{index(1)}{index(2)}{index(4)}.data(:,3));
    samprate = ripple{index(1)}{index(2)}{index(4)}.samprate;
    clear lowgamma valid_ripples ripples cellinfo mintrhesh excludetimes ripple
    
    %Define offset for each trace
    triggers = riptimes(:,2)-times(1);
    time = -0.4+1/1500:1/1500:0.4;
    while triggers(1) < time(end)
        triggers(1) = [];
    end
    while triggers(end) > (times(end)-times(1)-time(end))
        triggers(end) = [];
    end
    ldata = createdatamatc(l,triggers,1500,[0.4 0.4]);
    rdata = createdatamatc(r,triggers,1500,[0.4 0.4]);

    
    % Determine when each gamma trace reaches a baseline value
    time_zero = lookup(0,time);
    baseline = (ldata-mean(l))./std(l);
    
    for i = 1:length(triggers)
        [maximum minimum] = peakdet(baseline(:,i),0.5,time);
        out = [out minimum(lookup(0,minimum(:,1)),1)];
    end
end

end
