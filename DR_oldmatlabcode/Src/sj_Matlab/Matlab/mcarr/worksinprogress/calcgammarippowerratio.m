function out  =  calcgammarippowerratio(index,excludetimes,ripples,cellinfo,ripple,lowgamma,varargin)

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
    
    %Define gamma amplitude
    g = double(lowgamma{index(1)}{index(2)}{index(4)}.data(:,3));

    %Define valid times
    gtime = geteegtimes(lowgamma{index(1)}{index(2)}{index(4)});
    rtime = geteegtimes(ripple{index(1)}{index(2)}{index(3)});

	%Triggers
	r_triggers = [lookup(riptimes(:,1),rtime) lookup(riptimes(:,2),rtime)];
    g_triggers = [lookup(riptimes(:,1),gtime) lookup(riptimes(:,2),gtime)];

	ratio = zeros(size(riptimes(:,[1 2])));
    for i = 1:size(riptimes,1)
		ratio(i,1) = max(g(g_triggers(i,1):1:g_triggers(i,2)));
		ratio(i,2) = max(r(r_triggers(i,1):1:r_triggers(i,2)));
    end
	ratio = ratio(:,1)./ratio(:,2);
    out = ratio;

else
	out = [];
end

end
