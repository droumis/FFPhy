function out = calcripplecoherence(index,excludetimes,ripples,cellinfo,ripple,lowgamma,varargin)
%function out = calcripplecoherence(index,excludetimes,ripple,lowgamma,varargin)

%   This function goes through for each tetrode pair and computes the
%   correlation between ripple oscillations and gamma oscillation during
%   ripples.

%   out is structure with the following fields:
%       - ripple_corr
%       - gamma_corr


out.ripple_corr = []; out.gamma_corr = [];
%Define ripple times
riptimes = getripples(index([1 2]), ripples, cellinfo, 'cellfilter','(isequal($area, ''CA1''))');

e1 = ripple{index(1)}{index(2)}{index(3)}.data(:,1);
e1time = geteegtimes(ripple{index(1)}{index(2)}{index(3)});
ind1 = [lookup(riptimes(:,1),e1time) lookup(riptimes(:,2),e1time)];

e2 = ripple{index(1)}{index(2)}{index(4)}.data(:,1);
e2time = geteegtimes(ripple{index(1)}{index(2)}{index(4)}); clear ripple
ind2 = [lookup(riptimes(:,1),e2time) lookup(riptimes(:,2),e2time)];

for rip = 1:size(ind1,1)
	if length(e1(ind1(rip,1):ind1(rip,2)))== length(e2(ind2(rip,1):ind2(rip,2)))
    	tmpcorr = corr(double(e1(ind1(rip,1):ind1(rip,2))),double(e2(ind2(rip,1):ind2(rip,2))));
        out.ripple_corr = [out.ripple_corr tmpcorr];
    end
end


e1 = lowgamma{index(1)}{index(2)}{index(3)}.data(:,1);
e1time = geteegtimes(lowgamma{index(1)}{index(2)}{index(3)});
ind1 = [lookup(riptimes(:,1),e1time) lookup(riptimes(:,2),e1time)];

e2 = lowgamma{index(1)}{index(2)}{index(4)}.data(:,1);
e2time = geteegtimes(lowgamma{index(1)}{index(2)}{index(4)}); clear lowgamma
ind2 = [lookup(riptimes(:,1),e2time) lookup(riptimes(:,2),e2time)];

for rip = 1:size(ind1,1)
	if length(e1(ind1(rip,1):ind1(rip,2)))== length(e2(ind2(rip,1):ind2(rip,2)))
    	tmpcorr = corr(double(e1(ind1(rip,1):ind1(rip,2))),double(e2(ind2(rip,1):ind2(rip,2))));
        out.gamma_corr = [out.gamma_corr tmpcorr];
    end
end

end