function out = loadriptriggeredspectrum(index, excludetimes, spectrum, varargin)
% function out = loadriptriggeredspectrum(index, excludetimes, spectrum, varargin)
%
%  Loads the rip triggered spectrum for tetrodes defined by tetfilter and
%  organizes them nicely.
%
%   out is a structure with the following fields:
%       spectrum-- normalized ripple triggered spectrum.
%       time-- Time vector
%       frequency-- Frequency vector
%       index-- Only if appendindex is set to 1

% set options
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

%Apply excludetimes
included = ~isExcluded(spectrum{index(1)}{index(2)}{index(3)}.ripples, excludetimes);


out.spectrum = spectrum{index(1)}{index(2)}{index(3)}.spectrum(:,:,included);
out.time = spectrum{index(1)}{index(2)}{index(3)}.time;
out.frequency = spectrum{index(1)}{index(2)}{index(3)}.frequency;
out.ripples = spectrum{index(1)}{index(2)}{index(3)}.ripples(included);

if (appendindex)
    out.index = index;
end
