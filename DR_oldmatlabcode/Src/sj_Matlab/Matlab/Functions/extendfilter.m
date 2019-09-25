function f = extendfilter(f, varargin)
% function f = EXTENDFILTER(f, varargin)
% Take an existing filter and change only specified fields.  
% If you want to move the output to a new fieldname and clear output, use the 
% 'keepoutput' option as below.  Why do this?  You can run multiple analyses, 
% for example several sorts of shuffling, from one datafilter script and have the 
% outputs around for plotting, etc.  
%
% 'keepoutput'    'output_name'  This moves .output to .output_name
%                                AND clears the .output field.  
%
% Otherwise use options exactly as in createfilter.
%
% see also RUN_AND_EXTEND_FILTER, SETFILTERKEEPOUTPUT.  
% anathe may2009

for option = 1:2:length(varargin)-1
    
    switch lower (varargin{option})
    case 'animal'
        f = setfilteranimal(f,varargin{option+1});
    case 'epochs'
        f = setfilterepochs(f,varargin{option+1});          
    case 'excludetime'
        f = setfiltertime(f,varargin{option+1});
    case 'excludetimefilter'
        f = setfiltertime(f,varargin{option+1});
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
    case 'keepoutput'
	f = setfilterkeepoutput(f, varargin{option+1});
    otherwise
        error(['Input ''', varargin{option}, ''' not defined']);
    end

end







