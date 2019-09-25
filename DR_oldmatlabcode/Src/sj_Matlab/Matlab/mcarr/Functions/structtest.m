function out = structtest(index, excludeperiods, spikes, linpos, varargin)



model = [];
filter = [];
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'model'
                model = varargin{option+1};
            case 'filter'
                filter = varargin{option+1};    
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end


out.phrase = 'hello';