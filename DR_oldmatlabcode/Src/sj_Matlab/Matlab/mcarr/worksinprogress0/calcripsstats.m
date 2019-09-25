function out = calcripsstats(index, excludetimes, rip, varargin)
% function out = calcripsstats(index, excludetimes, rip, varargin)
%
%  Loads the rip structure and computes the average CA1 spectrum and the
%  average CA3 spectrum

% set options
appendindex = 1;
min_separation = 1;

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'min_separation'
                min_separation = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end
rip = rip{index(1)}{index(2)};
%Apply excludetimes
if ~isempty(rip)
    included = ~isExcluded(rip.starttime, excludetimes);
    valid = [1000; diff(rip.starttime)];
    valid = valid > min_separation;
    included = included & valid;
    out.starttime = rip.starttime;
    out.frequency = rip.frequency;
    out.time = rip.time;

    rip_ca1 = rip.ca1_spectrum(:,:,included);
    rip_ca3 = rip.ca3_spectrum(:,:,included);

    %Sum across valid ripples for average spectrum
    out.ca1_spectrum = sum(rip_ca1,3);
    out.ca3_spectrum = sum(rip_ca3,3);
    out.nrips = sum(included);

else
    out.starttime = [];
    out.frequency = [];
    out.time = [];
    out.ca1_spectrum = [];
    out.ca3_spectrum = [];
    out.nrips = [];
end

if (appendindex)
    out.index = index;
end

end