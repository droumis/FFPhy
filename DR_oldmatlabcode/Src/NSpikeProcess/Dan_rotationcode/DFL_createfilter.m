function f = JY_createfilter(days,varargin)

f.animal = [];
f.epochs = [];
f.data = [];
f.excludetime = [];
f.iterator = [];
f.function = [];
f.output = [];

for option = 1:2:length(varargin)-1
    
    switch varargin{option}
        case 'animal'
            f = DFL_setfilteranimal(f,varargin{option+1});
        case 'days'
            days = str2num(varargin{option+1});
        case 'epochs'
            f = JY_setfilterepochs(f,varargin{option+1},days);          
        case 'excludetime'
            f = JY_setfiltertime(f,varargin{option+1});
        case 'excludetimefilter'
            f = JY_setfiltertime(f,varargin{option+1});
        case 'excludetimelist'
            f = setexcludetime(f,varargin{option+1});
        case 'cells'
            f = setfiltercells(f,varargin{option+1});
        case 'cellpairs'
            f = setfiltercellpairs(f,varargin{option+1});
        case 'eegtetrodes'
            f = setfiltereegtetrodes(f, varargin{option+1});
        case 'eegtetrodepairs'
            f = setfiltereegtetrodepairs(f, varargin{option+1});
        case 'iterator'
            f = setfilteriterator(f, varargin{option+1});
        case 'filterfunction'
            f = setfilterfunction(f, varargin{option+1}{:});    
        otherwise
            error(['Input ''', varargin{option}, ''' not defined']);
    end
end