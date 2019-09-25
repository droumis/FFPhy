function out = loadriptriggeredcoherence(index, excludetimes, coherence, tetinfo, varargin)
% function out = loadriptriggeredcoherence(index, excludetimes, coherence, varargin)
%
%  Loads the rip triggered coherence for tetrodes defined by tetfilter,
%  computes the delta gamma coherence and delta gamma phase locking
%
%   out is a structure with the following fields:
%       coherence-- delta gamma coherence
%       phase-- delta gamma phase locking
%       index-- Only if appendindex is set to 1

% set options
appendindex = 1;
pair = 0;
hemisphere = 0;

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'pair'
                pair = varargin{option+1};
            case 'hemisphere'
                hemisphere = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

out.time = [];
out.frequency = [];
out.coherence = [];
out.coherence_baseline = [];
out.gamma_phase = [];
out.gamma_phase_baseline = [];
out.ripples = [];
    
%Determine which pairs of tetrodes are of interest
if isstr(pair)
    valid = evaluatefilter(tetinfo{index(1)}{index(2)},pair);
    if any(valid == index(3))
        valid(valid==index(3)) = [];
    end
else
    error(['Coherence requires definition of tetrode pair']);
end

%Apply hemisphere exclusion
if hemisphere
    tmp = evaluatefilter(tetinfo{index(1)}{index(2)},sprintf('~isequal($hemisphere,''%s'')',tetinfo{index(1)}{index(2)}{index(3)}.hemisphere));
    valid = intersect(valid,tmp);
end
%Apply excludetimes
if ~isempty(valid) && ~isempty(coherence)
    included = ~isExcluded(coherence{index(1)}{index(2)}{index(3)}{valid(1)}.ripples, excludetimes);
    valid_rip = [1000; diff(coherence{index(1)}{index(2)}{index(3)}{valid(1)}.ripples)];
    valid_rip = valid_rip > 1;
    included = included & valid_rip;    
    
    out.time = coherence{index(1)}{index(2)}{index(3)}{valid(1)}.time;
    out.frequency = coherence{index(1)}{index(2)}{index(3)}{valid(1)}.frequency;
    out.ripples = coherence{index(1)}{index(2)}{index(3)}{valid(1)}.ripples(included);

    p = [];  c = [];
    for v = 1:length(valid)
        try size(coherence{index(1)}{index(2)}{index(3)}{valid(v)}.coherence);
            if isempty(c)
                c = squeeze(mean(im2double(coherence{index(1)}{index(2)}{index(3)}{valid(v)}.coherence(:,3:9,included)),2));
                p = squeeze(mean(coherence{index(1)}{index(2)}{index(3)}{valid(v)}.phase(:,3:9,included),2));
            else
                c = c + squeeze(mean(im2double(coherence{index(1)}{index(2)}{index(3)}{valid(v)}.coherence(:,3:9,included)),2));
                p = cat(3,p,squeeze(mean(coherence{index(1)}{index(2)}{index(3)}{valid(v)}.phase(:,3:9,included),2)));
            end
        catch
            valid(v) = [];
        end
    end
    c = c./length(valid); p = mean(p,3);
    
    if ~isempty(c)
        %Determine delta coherence
        baseline = squeeze(mean(c(1:6,:),1));
        c = c(6:10:end,:) - repmat(baseline,size(c(6:10:end,:),1),1);
        coherence_baseline = mean(baseline);
        out.coherence = c;
        out.coherence_baseline = coherence_baseline;
    end
    
    if ~isempty(p)
        %Determine delta phase locking
        tmp_baseline = p(1:6,:); baseline = zeros(6,1);
        for j = 1:size(tmp_baseline,1)
            [m r] = anglemean(tmp_baseline(j,:));
            baseline(j) = r;
        end
        tmp_phase = p(6:10:end,:); phase = zeros(size(p(6:10:end,1)));
        for j = 1:size(tmp_phase,1)
           [m r] = anglemean(tmp_phase(j,:));
           phase(j) = r;
        end
        out.gamma_phase = phase- mean(baseline);
        out.gamma_phase_baseline = mean(baseline);
    end
end

if (appendindex)
    out.index = index;
end
