%{
meant to be an updated version on calcxcorrmeasures
that is compatible with my changes to singlecellanal (which is actually
feeding in multiple su).. 
.. i.e. the way data is fed in as a varargin

also see work i had done for mu xcorr in dfa_perripspikingcorr

also looks like i did something similar with: dfa_lickBoutSpikeCorr..
except that i took out the coactive z stuff. also it looks like i was using
that to run MU 

%}

function out = dfa_suCoactiveXcorr(idx, excludetimes, varargin)
%function out = calcxcorrmeasures(index, excludetimes, spikes, varargin)
% Calculates the excess correlation and RMS time lag for the specified cell
% pairs using only spikes not excluded by the excludetimes
% 
%
% Options:
%   'bin', n 	binsize in sec. Default 0.002 (2 ms)
%   'tmax', n	maximum time for cross correlation. Default 1 sec.
%   'edgespikes', 0 or 1  
%		  0 will enforce the exclude periods and will thus
%   	    	    exclude all spikes outside of the exclude windows. (default)
%		  1 will include all spikes where at least one of the spikes is
%		    not excluded 
%   'sw1', n	The smaller smoothing width for the excess correlation measure
%		Default 0.005 (5 ms)
%   'sw2', n	The larger smoothing width for the excess correlation measure
%		Default 0.250 (250 ms)
%   'rmstmax', n The maximum absolute time to be used for the rms calculation
%		Default 0.1 (100 ms)
%   'rmsmincounts', n The minimum number of counts in the histogram for the rms
%  		measurement to be carried out. 
%		Default 50. 
%   'calclinfields', 1 or 0 -- set to 1 to calculate linearized place fields
%   'calctrajxcorr', 1 or 0 -- set to 1 to calculate one cross correlation per
%   				trajectory (default 0);
%

calclinfields = 0;
calctrajxcorr = 0;
appendindex = 0;
bin = 0.002;
tmax = 1;
edgespikes = 0;
sw1 = 0.005;
sw2 = 0.250;
rmstmax = 0.1;
rmsmincounts = 50;

if ~isempty(varargin)
    assign(varargin{:})
end

% set the defaults
out.time = [];
out.c1vsc2 = [];
out.ac1 = [];
out.ac2 = [];
out.index = idx;
out.ec = NaN;
out.rms = NaN;
out.probcoa = NaN;
out.coactivez_numsp = NaN;
out.coactivez = NaN;
out.lf1 = [];
out.lf2 = [];

% for each cell we calculate the cross correlation 
try
    t1 = spikes{idx(1)}{idx(2)}{idx(3)}{idx(4)}.data;
    t2 = spikes{idx(1)}{idx(2)}{idx(5)}{idx(6)}.data;
catch
    % if either of those produced an error, we return NaNs 
    return;
end

%apply the exclude rule
t1inc = [];
t2inc = [];
if (length(t1))
    t1inc = t1(find(~isExcluded(t1(:,1), excludetimes)),1);
else
    return
end
if (length(t2))
    t2inc = t2(find(~isExcluded(t2(:,1), excludetimes)),1);
else
    return
end

ac1 = spikexcorr(t1inc, t1inc, bin, tmax);
ac2 = spikexcorr(t2inc, t2inc, bin, tmax);
out.ac1 = ac1.c1vsc2;
out.ac2 = ac2.c1vsc2;
out.time = ac1.time;

out.index = idx;
if (~calctrajxcorr)
    xc = spikexcorr(t1inc, t2inc, bin, tmax);
    % if we want to include edge spikes, we need to add in the correlation of 
    % the excluded t1 spikes with the included t2spikes
    if ((edgespikes) & (~isempty(xc.time)))
	t1ex = t1(find(isExcluded(t1(:,1), excludetimes)));
	if (~isempty(t1ex))
	    tmpxc = spikexcorr(t1ex, t2inc, bin, tmax);
	    % add these values to the original histogram
	    xc.c1vsc2 = xc.c1vsc2 + tmpxc.c1vsc2;
	end
    end
    out.c1vsc2 = xc.c1vsc2;
    % compute the excess correlation at 0 lag
    out.ec = excesscorr(xc.time, xc.c1vsc2, xc.nspikes1, xc.nspikes2, sw1, sw2);
    out.rms = xcorrrms(xc.time, xc.c1vsc2, rmstmax, rmsmincounts);

    % get the probability that both were coactive in the included intervals
    n1 = nInInterval(t1inc, excludetimes, 0);
    n2 = nInInterval(t2inc, excludetimes, 0);
    out.probcoa = mean(logical(n1) & logical(n2));

    % get the zscore for this probabliity of coactivity
    out.coactivez = coactivezscore(n1, n2);
    out.coactivez_numsp = coactivezscore(n1, n2, 'anyspikes', 0);
else
    % calculated individual trajectory cross correlations 
    statematrix = linpos{idx(1)}{idx(2)}.statematrix;
    statevector = statematrix.traj;
    statevector(find(isExcluded(statematrix.time, excludetimes))) = -1;
    % look up the state for each spike
    posindexfield = 7;
    t1 = t1(:,[1 posindexfield]);
    t2 = t2(:,[1 posindexfield]);
    t1(:,3) = statevector(t1(:,2));
    t2(:,3) = statevector(t2(:,2));
    % get rid of -1 trajectories and excluded times 
    t1 = t1(find(t1(:,3) ~= -1),:);
    t2 = t2(find(t2(:,3) ~= -1),:);
    t2 = t2(find(~isExcluded(t2(:,1), excludetimes)), :);
    for i = 1:max(statevector)
	t1tmp = t1(find(t1(:,3) == i), 1);
	t1inc = t1tmp(find(~isExcluded(t1tmp, excludetimes)));
	t2tmp = t2(find(t2(:,3) == i), 1);
	t2inc = t2tmp(find(~isExcluded(t2tmp, excludetimes)));
	xc = spikexcorr(t1inc, t2inc, bin, tmax);
	if ((edgespikes) & (~isempty(xc.time)))
	    t1ex = t1tmp(find(isExcluded(t1tmp, excludetimes)));
	    if (~isempty(t1ex))
		tmpxc = spikexcorr(t1ex, t2inc, bin, tmax);
		% add these values to the original histogram
		xc.c1vsc2 = xc.c1vsc2 + tmpxc.c1vsc2;
	    end
	end
	out.c1vsc2{i} = xc.c1vsc2;
	% compute the excess correlation at 0 lag
	out.ec(i) = excesscorr(xc.time, xc.c1vsc2, xc.nspikes1, xc.nspikes2, sw1, sw2);
	out.rms(i) = xcorrrms(xc.time, xc.c1vsc2, rmstmax, rmsmincounts);

	% get the probability that both were coactive in the included intervals
	n1 = nInInterval(t1inc, excludetimes, 0);
	n2 = nInInterval(t2inc, excludetimes, 0);
	out.probcoa(i) = mean(logical(n1) & logical(n2));

	% get the zscore for this probabliity of coactivity
    out.coactivez(i) = coactivezscore(n1, n2);
    out.coactivez_numsp(i) = coactivezscore(n1, n2, 'anyspikes', 0);
    end
end

if (calclinfields)
    % calculate the linear fields 
    out.lf1 = filtercalclinfields(idx(1:4), excludetimes, spikes, linpos, 'phasedist', 1);
    out.lf2 = filtercalclinfields([idx(1:2) idx(5:6)], excludetimes, spikes, linpos, 'phasedist', 1);
end


