function out = calccrossfrequencycoupling(index,excludetimes,ripples,cellinfo,ripple,lowgamma,varargin)
%function out = calccrossfrequencycoupling(index,excludetimes,ripple,lowgamma,varargin)

%Set options
cellfilter = [];
minthresh = 3;

%Process options
for option = 1:2:length(varargin)-1   
    if ischar(varargin{option})
        switch(varargin{option})
            case 'cellfilter'
                cellfilter = varargin{option+1};         
            case 'minthresh'
                minthresh = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

if length(index) == 3
    index = [index index(3)];
end
%Find valid ripples
if isempty(cellfilter)
    riptimes = getripples([index(1,1) index(1,2)], ripples, cellinfo, 'cellfilter', ...
        '(isequal($area, ''CA1''))','excludeperiods', excludetimes,'minstd',minthresh);
else
    riptimes = getripples([index(1,1) index(1,2)], ripples, cellinfo, 'cellfilter', ...
        cellfilter,'excludeperiods', excludetimes,'minstd',minthresh);
end

if ~isempty(riptimes)
    %Exclude ripples occuring close together for preceding modulationg
    valid_ripples = [1000; riptimes(2:end,2)-riptimes(1:end-1,1)];
    valid_ripples = valid_ripples > 1;
    riptimes = riptimes(valid_ripples,:);
    
    %Define ripple amplitude
    r = double(ripple{index(1)}{index(2)}{index(3)}.data(:,3));
    r = (r-mean(r))./std(r);    %Compute z-score
    
    %Define gamma signal and phase
    g = double(lowgamma{index(1)}{index(2)}{index(4)}.data(:,[1 2]));
    
    %Define valid times
    gtime = geteegtimes(lowgamma{index(1)}{index(2)}{index(4)});
    rtime = geteegtimes(ripple{index(1)}{index(2)}{index(3)});
    if length(gtime)>length(rtime)
        temp = lookup(rtime,gtime);
        g = g(temp,:);
        temp = logical(isExcluded(rtime,riptimes(:,[1 2])));
        g = g(temp,:);
        r = r(temp);
    elseif length(rtime)>length(gtime)
        temp = lookup(gtime,rtime);
        r = r(temp);
        temp = logical(isExcluded(gtime,riptimes(:,[1 2])));
        g = g(temp,:);
        r = r(temp);
    elseif length(gtime)==length(rtime)
        temp = logical(isExcluded(gtime,riptimes(:,[1 2])));
        g = g(temp,:);
        r = r(temp);
    end
    clear temp
    
    out.gamma = g(:,1);
    out.phase = g(:,2);
    out.amp = r;
else
    out.gamma = [];
    out.phase = [];
    out.amp = [];
end