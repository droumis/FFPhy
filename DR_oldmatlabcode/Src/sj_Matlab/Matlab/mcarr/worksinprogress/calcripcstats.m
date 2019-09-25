function out = calcripcstats(index, excludetimes, rip, varargin)
% function out = calcripcstats(index, excludetimes, rip, varargin)
%
%  Loads the ripc structure and computes the average CA1-CA1, CA1-CA3, and
%  CA3-CA3 coherence.
%

% set options
appendindex = 1;
min_separation = 1;
bin = -pi+pi/40:pi/20:pi-pi/40;

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
if ~isempty(rip) & isfield(rip,'frequency')
    included = ~isExcluded(rip.starttime, excludetimes);
    valid = [1000; diff(rip.starttime)];
    valid = valid > min_separation;
    included = included & valid;
    out.starttime = rip.starttime;
    out.frequency = rip.frequency;
    out.time = rip.time;

  	rip_ca11 = rip.ca1_ca1_coherence(:,:,included);
    rip_ca13 = rip.ca1_ca3_coherence(:,:,included);
    rip_ca33 = rip.ca3_ca3_coherence(:,:,included);
        
    %Sum across valid ripples for average coherence
    out.ca1_ca1_coherence = sum(rip_ca11,3);
    out.ca1_ca3_coherence = sum(rip_ca13,3);
    out.ca3_ca3_coherence = sum(rip_ca33,3);
    out.nrips = sum(included);
    
    
    %Determine peak coherence frequency
    starttime = lookup(0,out.time);
    out.peak11 = zeros(size(rip_ca11,3),1);
    out.peak13 = zeros(size(rip_ca11,3),1);
    out.peak33 = zeros(size(rip_ca11,3),1);
    for i = 1:size(rip_ca11,3)
        val = max(rip_ca11(starttime:end,:,i));
        [val tmp] = max(val);
        out.peak11(i) = out.frequency(tmp);
        
        val = max(rip_ca13(starttime:end,:,i));
        [val tmp] = max(val);
        out.peak13(i) = out.frequency(tmp);
        
        val = max(rip_ca33(starttime:end,:,i));
        [val tmp] = max(val);
        out.peak33(i) = out.frequency(tmp);
    end
    %Determine change in phase locking
    if ~isempty(rip.ca1_ca1_phase)
        p = rip.ca1_ca1_phase(:,included);
        tmp = zeros(size(p,1),2);
        for i = 1:size(p,1)
             [m r] = anglemean(p(i,:));
             tmp(i,:) = [m r];
        end
        out.ca1_ca1_phase_std = tmp;
    else
        out.ca1_ca1_phase_std = [];
    end
    
    if ~isempty(rip.ca1_ca3_phase)
        p = rip.ca1_ca3_phase(:,included);
        tmp = zeros(size(p,1),2);
        for i = 1:size(p,1)
            [m r] = anglemean(p(i,:));
            tmp(i,:) = [m r];
        end
        out.ca1_ca3_phase_std = tmp;
    else
        out.ca1_ca3_phase_std = [];
    end
    
     if ~isempty(rip.ca1_ca3_phase)
        p = rip.ca3_ca3_phase(:,included);
        tmp = zeros(size(p,1),2);
        for i = 1:size(p,1)
            [m r] = anglemean(p(i,:));
            tmp(i,:) = [m r];
        end
        out.ca3_ca3_phase_std = tmp;
    else
        out.ca3_ca3_phase_std = [];
    end
    
else
    out.starttime = [];
    out.frequency = [];
    out.time = [];
    out.ca1_ca1_coherence = [];
    out.ca1_ca3_coherence = [];
    out.ca3_ca3_coherence = [];
	out.nrips = [];
    out.peak11 = [];
    out.peak13 = [];
    out.peak33 = [];
    out.ca1_ca1_phase_std = [];
    out.ca1_ca3_phase_std = [];
    out.ca3_ca3_phase_std = [];
end

if (appendindex)
    out.index = index;
end

end