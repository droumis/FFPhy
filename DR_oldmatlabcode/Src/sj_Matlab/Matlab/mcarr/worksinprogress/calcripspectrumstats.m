function out = calcripspectrumstats(index, excludetimes, spectrum, varargin)
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
appendindex = 1;
gam = [20 50];
rip = [150 250];

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'gamma_frequency'
                gam = varargin{option+1};
            case 'ripple_frequency'
                rip = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

%Apply excludetimes
if ~isempty(spectrum)
    if length(index) == 4
        if length(spectrum{index(1)}{index(2)}{index(3)}.ripples)==...
                length(spectrum{index(1)}{index(2)}{index(4)}.ripples)
            included = ~isExcluded(spectrum{index(1)}{index(2)}{index(3)}.ripples, excludetimes);
            s = spectrum{index(1)}{index(2)}{index(3)}.spectrum(:,:,included);
            srip = spectrum{index(1)}{index(2)}{index(4)}.spectrum(:,:,included);
        elseif length(spectrum{index(1)}{index(2)}{index(3)}.ripples)<...
                length(spectrum{index(1)}{index(2)}{index(4)}.ripples)
            valid = lookup(spectrum{index(1)}{index(2)}{index(3)}.ripples,spectrum{index(1)}{index(2)}{index(4)}.ripples);
            included = ~isExcluded(spectrum{index(1)}{index(2)}{index(3)}.ripples, excludetimes);
            s = spectrum{index(1)}{index(2)}{index(3)}.spectrum(:,:,included);
            srip = spectrum{index(1)}{index(2)}{index(4)}.spectrum(:,:,valid);
            srip = srip(:,:,included);
        else
            valid = lookup(spectrum{index(1)}{index(2)}{index(4)}.ripples,spectrum{index(1)}{index(2)}{index(3)}.ripples);
            included = ~isExcluded(spectrum{index(1)}{index(2)}{index(3)}.ripples(valid), excludetimes);
            s = spectrum{index(1)}{index(2)}{index(3)}.spectrum(:,:,valid);
            s = s(:,:,included);
            srip = spectrum{index(1)}{index(2)}{index(4)}.spectrum(:,:,included);
            
        end
    else
        included = ~isExcluded(spectrum{index(1)}{index(2)}{index(3)}.ripples, excludetimes);
        s = spectrum{index(1)}{index(2)}{index(3)}.spectrum(:,:,included);
        srip = s;
    end
    out.time = spectrum{index(1)}{index(2)}{index(3)}.time;
    out.freq = spectrum{index(1)}{index(2)}{index(3)}.frequency;

    %Sum across valid ripples for average spectrum
    out.total = sum(s,3);
    out.nrips = sum(included);
    
    %Compute average gamma power before, during, and after peak of ripple
    %power and compute cross correlation between ripple and gamma power
    out.peak = zeros(size(s,3),1);
    out.power = zeros(size(s,3),3);
    rip_ind = lookup(rip,out.freq);
    gam_ind = lookup(gam,out.freq);

    for i = 1:size(s,3);
        out.peak(i) = max(max(srip(:,rip_ind(1):rip_ind(end),i),[],2));
        peak_ind = find(out.peak(i) == max(srip(:,rip_ind(1):rip_ind(end),i),[],2));
        out.power(i,:) = mean(s([1 peak_ind end],gam_ind(1):gam_ind(end),i),2);
    end
else
    out.time = [];
    out.freq = [];
    out.total = [];
    out.nrips = [];
    out.peak = [];
    out.power = [];
end

if (appendindex)
    out.index = index;
end

end