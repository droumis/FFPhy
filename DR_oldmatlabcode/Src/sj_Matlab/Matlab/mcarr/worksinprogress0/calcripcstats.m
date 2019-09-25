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
if ~isempty(rip)
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
    
    %Determine change in variance for phase
    p11 = rip.ca1_ca1_phase(:,included);
    p13 = rip.ca1_ca3_phase(:,included);
    p33 = rip.ca3_ca3_phase(:,included);
    
    tmp11 = zeros(size(p11,1),1);
    tmp13 = zeros(size(tmp11));
    tmp33 = zeros(size(tmp11));
    for i = 1:size(p11,1)
        tmp11(i) = anglestd(p11(i,:));
        tmp13(i) = anglestd(p13(i,:));
        tmp33(i) = anglestd(p33(i,:));
    end
    
    out.ca1_ca1_phase_std = tmp11;
    out.ca1_ca3_phase_std = tmp13;
    out.ca3_ca3_phase_std = tmp33;

   % bin = -pi+pi/40:pi/20:pi-pi/40;
    %Could histogram the phase over time to make plots.
else
    out.starttime = [];
    out.frequency = [];
    out.time = [];
    out.ca1_ca1_coherence = [];
    out.ca1_ca3_coherence = [];
    out.ca3_ca3_coherence = [];
    out.nrips = [];
    out.ca1_ca1_phase_std = [];
    out.ca1_ca3_phase_std = [];
    out.ca3_ca3_phase_std = [];
end

if (appendindex)
    out.index = index;
end

end