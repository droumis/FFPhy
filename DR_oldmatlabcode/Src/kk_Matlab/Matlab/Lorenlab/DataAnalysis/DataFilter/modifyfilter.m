function f = modifyfilter(f,varargin)

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'animal'
            [f.animal] = deal([]);
            [f.epochs] = deal([]);
            [f.excludetime] = deal([]);
            [f.data] = deal([]);
            [f.output] = deal([]);
            f = setfilteranimal(f,varargin{option+1});
        case 'epochs'
            [f.epochs] = deal([]);
            [f.excludetime] = deal([]);
            [f.data] = deal([]);
            [f.output] = deal([]);
            f = setfilterepochs(f,varargin{option+1});          
        case 'excludetime'
            [f.excludetime] = deal([]);
            [f.output] = deal([]);
            f = setfiltertime(f,varargin{option+1});
        case 'excludetimefilter'
            [f.excludetime] = deal([]);
            [f.output] = deal([]);
            f = setfiltertime(f,varargin{option+1});
        case 'dio'
            [f.output] = deal([]);
            f = setdiofilter(f,varargin{option+1});
        case 'stim'
            [f.output] = deal([]);
            f = setstimfilter(f,varargin{option+1});
        case 'cells'
            [f.data] = deal([]);
            [f.output] = deal([]);
            f = setfiltercells(f,varargin{option+1});
        case 'cellpairs'
            [f.data] = deal([]);
            [f.output] = deal([]);
            f = setfiltercellpairs(f,varargin{option+1});
        case 'eegtetrodes'
            [f.data] = deal([]);
            [f.output] = deal([]);
            f = setfiltereegtetrodes(f, varargin{option+1});
        case 'tetrodes'
            [f.data] = deal([]);
            [f.output] = deal([]);
            f = setfiltertetrodes(f, varargin{option+1});
        case 'eegtetrodepairs'
            [f.data] = deal([]);
            [f.output] = deal([]);
            f = setfiltereegtetrodepairs(f, varargin{option+1});
        case 'iterator'
            [f.iterator] = deal([]);
            [f.output] = deal([]);
            f = setfilteriterator(f, varargin{option+1});
        case 'filterfunction'
            [f.filterfunction] = deal([]);
            [f.output] = deal([]);
            f = setfilterfunction(f, varargin{option+1}{:});    
        otherwise
            error(['Input ''', varargin{option}, ''' not defined']);
    end
end

