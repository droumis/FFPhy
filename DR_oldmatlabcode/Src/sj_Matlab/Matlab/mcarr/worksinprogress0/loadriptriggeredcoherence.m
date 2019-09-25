function out = loadriptriggeredcoherence(index, excludetimes, coherence, tetinfo, varargin)
% function out = loadriptriggeredcoherence(index, excludetimes, coherence, varargin)
%
%  Loads the rip triggered coherence for tetrodes defined by tetfilter and
%  organizes them nicely.
%
%   out is a structure with the following fields:
%       coherence--  ripple triggered coherence
%       time-- Time vector
%       frequency-- Frequency vector
%       index-- Only if appendindex is set to 1

% set options
appendindex = 0;
pair = 0;
bin = -pi+pi/40:pi/20:pi-pi/40;

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'pair'
                pair = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

%Determine which pairs of tetrodes are of interest
if isstr(pair)
    valid = evaluatefilter(tetinfo{index(1)}{index(2)},pair);
    if any(valid == index(3))
        valid(valid==index(3)) = [];
    end
else
    error(['Coherence requires definition of tetrode pair']);
end

%Apply excludetimes
if ~isempty(valid)
    included = ~isExcluded(coherence{index(1)}{index(2)}{index(3)}{valid(1)}.ripples, excludetimes);
    
    out.time = coherence{index(1)}{index(2)}{index(3)}{valid(1)}.time;
    out.frequency = coherence{index(1)}{index(2)}{index(3)}{valid(1)}.frequency;
    
    p = zeros(length(out.time),length(bin));  c = [];
    for v = 1:length(valid)
        try size(coherence{index(1)}{index(2)}{index(3)}{valid(v)}.coherence);
            if isempty(c)
                c = im2double(coherence{index(1)}{index(2)}{index(3)}{valid(v)}.coherence(:,:,included));
            else
                c = c + im2double(coherence{index(1)}{index(2)}{index(3)}{valid(v)}.coherence(:,:,included));
            end
        catch
            valid(v) = [];
        end
        for i = 1:length(out.time)
            p(i,:) = p(i,:) + hist(squeeze(mean(coherence{index(1)}{index(2)}{index(3)}{valid(v)}.phase(i,2:8,included),2)),bin);
        end
    end

    out.coherence = c./length(valid);
    out.gamma_phase = p;
    out.ripples = coherence{index(1)}{index(2)}{index(3)}{valid(1)}.ripples(included);
else
    out.time = [];
    out.frequency = [];
    out.coherence = [];
    out.gamma_phase = [];
    out.ripples = [];
end

if (appendindex)
    out.index = index;
end
