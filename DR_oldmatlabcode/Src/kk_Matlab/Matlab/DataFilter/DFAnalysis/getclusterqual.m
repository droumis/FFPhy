function out = getclusterqual(index, excludeperiods, clustqual, varargin)


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

% set out to include the lratio and isolation distance
out = [clustqual{index(1)}{index(2)}{index(3)}{index(4)}.lratio ...
       clustqual{index(1)}{index(2)}{index(3)}{index(4)}.isoldist];

if (appendindex)
    out = [index out];
end
