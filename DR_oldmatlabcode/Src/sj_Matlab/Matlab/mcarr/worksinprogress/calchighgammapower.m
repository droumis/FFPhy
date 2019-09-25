function out = calchighgammapower(index, excludetimes, rip, varargin)
% function out = calcripsstats(index, excludetimes, rip, varargin)
%
%  Loads the rips structure and computes the high gamma power in CA1 and CA3 

% set options
appendindex = 1;
min_separation = 1;
frequency = [60 90];
lowrip = [100 130];

out.starttime = [];
out.time = [];
out.ca1_lowgamma = [];
out.ca3_lowgamma = [];
out.ca1_highgamma = [];
out.ca3_highgamma = [];
out.ca1_lowrip = [];
out.ca3_lowrip = [];

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
    lowrip = lookup(lowrip,rip.frequency);
    out.ca1_lowgamma = squeeze(mean(rip.ca1_spectrum(6:10:end,3:9,included),2));
    out.ca3_lowgamma = squeeze(mean(rip.ca3_spectrum(6:10:end,3:9,included),2));
    
    out.ca1_highgamma = squeeze(mean(rip.ca1_spectrum(6:10:end,freq(1):freq(2),included),2));
    out.ca3_highgamma = squeeze(mean(rip.ca3_spectrum(6:10:end,freq(1):freq(2),included),2));
    
    out.ca1_lowrip = squeeze(mean(rip.ca1_spectrum(6:10:end,lowrip(1):lowrip(2),included),2));
    out.ca3_lowrip = squeeze(mean(rip.ca3_spectrum(6:10:end,lowrip(1):lowrip(2),included),2));
    
else

end

if (appendindex)
    out.index = index;
end

end