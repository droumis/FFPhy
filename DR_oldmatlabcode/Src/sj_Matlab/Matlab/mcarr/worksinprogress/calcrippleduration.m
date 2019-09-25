function out = calcrippleduration(index, excludeperiods, ripples,cellinfo, varargin)
% Computes the length of ripples detected using the exculdedperiods

%Set options
cellfilter = [];
minthresh = 3;
out = [];

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

%Find valid ripples
if isempty(cellfilter)
    riptimes = getripples([index(1,1) index(1,2)], ripples, cellinfo, 'cellfilter', ...
        '(isequal($area, ''CA1''))',...
        'excludeperiods', excludeperiods,'minstd',minthresh);
else
    riptimes = getripples([index(1,1) index(1,2)], ripples, cellinfo, 'cellfilter', ...
        cellfilter,'excludeperiods', excludeperiods,'minstd',minthresh);
end


if ~isempty(riptimes)
    valid_ripples = [1000; riptimes(2:end,2)-riptimes(1:end-1,1)];
    valid_ripples = valid_ripples > 1;
    riptimes = riptimes(valid_ripples,:);

    if ~isempty(riptimes)
        out = (riptimes(:,2)-riptimes(:,1))';
    end
end

end
