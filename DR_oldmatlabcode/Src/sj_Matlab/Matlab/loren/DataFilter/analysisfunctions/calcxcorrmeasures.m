function out = calcxcorrmeasures(ind, excludetimes, spikes, linfields, varargin)
%function out = calcxcorrmeasures(index, excludetimes, spikes, varargin)
% Calculates the excess correlation and RMS time lag for the specified cell
% pairs using only spikes not excluded by the excludetimes
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
%		Default 10.
%   'plotxcorr', 1 or 0, default 0
%       plot xcorrelation if plotxcorr = 1
%   'appendindex', 1 or 0 -- set to 1 to append the cell index to the
%   output [tetrode cell value].  Default 0.
%

%load /data/loren/mkarlsso/Ten/tenlinfields
%linfields = trajdata;

appendindex = 0;
bin = 0.002;
tmax = 1;
edgespikes = 0;
sw1 = 0.005;
sw2 = 0.250;
rmstmax = 0.1;
rmsmincounts = 10;
plotxcorr = 0;
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'bin'
                bin = varargin{option+1};
            case 'tmax'
                tmax = varargin{option+1};
            case 'edgespikes'
                edgespikes = varargin{option+1};
            case 'sw1'
                sw1 = varargin{option+1};
            case 'sw2'
                sw2 = varargin{option+1};
            case 'rmstmax'
                rmstmax = varargin{option+1};
            case 'plotxcorr'
                plotxcorr = varargin{option+1};
            case 'rmsmincounts'
                rmsmincounts = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

% for each cell we calculate the cross correlation
try
    t1 = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data(:,1);
    t2 = spikes{ind(1)}{ind(2)}{ind(5)}{ind(6)}.data(:,1);
catch
    % if either of those produced an error, we return NaN NaN
    out = [NaN NaN];
    return;
end

%apply the exclude rule
t1inc = t1(find(~isExcluded(t1, excludetimes)));
t2inc = t2(find(~isExcluded(t2, excludetimes)));

xc = spikexcorr(t1inc, t2inc, bin, tmax);
% if we want to include edge spikes, we need to add in the correlation of the
% excluded t1 spikes with the included t2spikes
if ((edgespikes) & (~isempty(xc.time)))
    t1ex = t1(find(isExcluded(t1, excludetimes)));
    if (~isempty(t1ex))
        tmpxc = spikexcorr(t1ex, t2inc, bin, tmax);
        % add these values to the original histogram
        xc.c1vsc2 = xc.c1vsc2 + tmpxc.c1vsc2;
    end
end

%Uncomment to plot each cross-correlation histogram

% compute the excess correlation at 0 lag
ec = excesscorr(xc.time, xc.c1vsc2, xc.nspikes1, xc.nspikes2, sw1, sw2);
rms = xcorrrms(xc.time, xc.c1vsc2, rmstmax, rmsmincounts);

% get the probability that both were coactive in the included intervals
n1 = nInInterval(t1inc, excludetimes, 0);
n2 = nInInterval(t2inc, excludetimes, 0);
probcoa = mean(logical(n1) & logical(n2));

% get the zscore for this probability of coactivity
coactivez = coactivezscore(n1, n2);

%compute overlap
lf1 = linfields{ind(1)}{ind(2)}{ind(3)}{ind(4)}; %lf1 = all traj
lf2 = linfields{ind(1)}{ind(2)}{ind(5)}{ind(6)};
overlap = calcoverlap(lf1,lf2);

%plot
if plotxcorr == 1;
    if (sum(xc.c1vsc2) > 10) %min number coincident event
        figure(1);
        g = gaussian(3, 18);
        plot(xc.time, smoothvect(xc.c1vsc2, g));
        title(sprintf('%d %d - %d %d vs. %d %d', ind))
        hold off
        figure(2)
%        lf1 = linfields{ind(1)}{ind(2)}{ind(3)}{ind(4)}; %lf1 = all traj
%        lf2 = linfields{ind(1)}{ind(2)}{ind(5)}{ind(6)};
        %len1 = min(length(lf1{1}(:,5)), length(lf2{1}(:,5)));
        %len2 = min(length(lf1{3}(:,5)), length(lf2{3}(:,5)));
        %tmp1 = [];
        %tmp2 = [];
        for i = 1:length(lf1) %for each trajectory
            %  if (i <= 2)
            %     l = len1;
            %  else
            %      l = len2;
            %  end
            % get rid of NaNs
            if ~isempty(lf1{i})
                lf1{i}(find(~isfinite(lf1{i}(:,5))), 5) = 0;
                lf2{i}(find(~isfinite(lf2{i}(:,5))), 5) = 0;
            end
            %tmp1 = [tmp1 ; lf1{i}(1:l,5)];
            %tmp2 = [tmp2 ; lf2{i}(1:l,5)];
        end
        for i = 1:length(lf1)
            if ~isempty(lf1{i})
                subplot(length(lf1),1,i)
                hold on
                plot(lf1{i}(:,1), lf1{i}(:,5));
                plot(lf2{i}(:,1), lf2{i}(:,5), 'r');
                subtitle(['overlap ', num2str(overlap), ' excess correlation ', num2str(ec)])
                hold off
            end
        end
        'pausing, strike key to continue'
        pause
        clf(1); %clear figures
        clf(2);
        %close all
    end
end


if (appendindex) 
    out = [ind ec rms probcoa coactivez overlap];
else
    out = [ec rms probcoa coactivez overlap];
end