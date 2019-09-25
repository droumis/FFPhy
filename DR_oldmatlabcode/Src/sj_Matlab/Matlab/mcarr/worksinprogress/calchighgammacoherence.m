function out = calchighgammacoherence(index, excludetimes, rip, varargin)
% function out = calcripsstats(index, excludetimes, rip, varargin)
%
%  Loads the ripc structure and computes the high gamma coherence between CA1 and CA3 

% set options
appendindex = 1;
min_separation = 1;
frequency = [60 90];
lowrip = [100 130];

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'min_separation'
                min_separation = varargin{option+1};
            case 'frequency'
                frequency = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

rip = rip{index(1)}{index(2)};
%Apply excludetimes
if ~isempty(rip.starttime)
    included = ~isExcluded(rip.starttime, excludetimes);
    valid = [1000; diff(rip.starttime)];
    valid = valid > min_separation;
    included = included & valid;
    out.starttime = rip.starttime;
    out.time = rip.time(6:10:end);

    freq = lookup(frequency,rip.frequency);
    out.ca1_ca3 = squeeze(mean(rip.ca1_ca3_coherence(6:10:end,freq(1):freq(2),included),2));
    out.ca1_ca3_baseline = mean(squeeze(mean(rip.ca1_ca3_coherence(1:6,freq(1):freq(2),included),2)));
    out.ca1_ca3_ripple = squeeze(mean(rip.ca1_ca3_coherence(6:10:end,26:43,included),2));
    out.ca1_ca3_ripple_baseline = mean(squeeze(mean(rip.ca1_ca3_coherence(1:6,26:43,included),2)));
    lowrip = lookup(lowrip,rip.frequency);
    out.ca1_ca3_lowripple = squeeze(mean(rip.ca1_ca3_coherence(6:10:end,lowrip(1):lowrip(2),included),2));
    out.ca1_ca3_lowripple_baseline = mean(squeeze(mean(rip.ca1_ca3_coherence(1:6,lowrip(1):lowrip(2),included),2)));
    
else
    out.starttime = [];
    out.time = [];
    out.ca1_ca3 = [];
    out.ca1_ca3_baseline = [];
    out.ca1_ca3_ripple = [];
    out.ca1_ca3_ripple_baseline = [];
    out.ca1_ca3_lowripple = [];
    out.ca1_ca3_lowripple_baseline = [];
end

if (appendindex)
    out.index = index;
end

end