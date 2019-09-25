function out = DFAsj_getrip_period_xcorr4(ind, excludetimes, spikes, ripples, tetinfo, pos, varargin)
% Version 4 - Sync with all
% Get crosscorrln for cell pairs during ripple periods.
% 1 sec around ripples is default, can specify time period around ripples
% Can also b e used to just return spiking data in ripple periods, for xcorr to be computed by script
% after combinging data across epochs
% Ripples defined as in DFAsj_getripalignspiking4
% Corrln as in DFAsj_getthetacrosscov, or previous xcorr code

% out = DFAsj_getrip_period_xcorr4(spike_index, excludeperiods, spikes, ripples, tetinfo, options)
% Called from DFSsj_HPexpt_Lags_ThetaAndRip_ver4

% Use tetinfo and tetfilter passed in, or redefine here to get riptets
% Then use ripples to getriptimes. Use inter-ripple-interval of 1 sec, and use a low-speed criterion.
% Then align spikes to ripples

appendindex = 1; % always appendindex


% Ripple parameters
% -----------------
tetfilter = '';
excludetimes = [];
maxcell = 0;
minstd = 3; 
lowsp_thrs = 5; %cm/sec
highsp_thrs = lowsp_thrs;
dospeed = 0;
minrip=1;
% For ripple period
% ------------------
pret=500/1000; postt=500/1000; % in secs
%cellcountthresh = 3;  % Can be used to parse ripples - Not used

% Corrln parameters
%------------------
bin = 0.01; % 10 ms bins            %bin = 0.002; % 2 ms bins            %
sw1 = bin*3;                      %sw1 = bin*3; % 6 ms smoothing      % for smoothing corrln. Necessary?
tmax = 0.5;
edgespikes = 0;


for option = 1:2:length(varargin)-1
    switch varargin{option}
        % Ripple parameters
        case 'tetfilter'
            tetfilter = varargin{option+1};
        case 'excludetimes'
            excludetimes = varargin{option+1};
        case 'minstd'
            minstd = varargin{option+1};
        case 'minrip'
            minrip = varargin{option+1};
        case 'maxcell'
            maxcell = varargin{option+1};
        case 'dospeed'
            dospeed = varargin{option+1};
        case 'lowsp_thrs'
            lowsp_thrs = varargin{option+1};
        % Corrln parameters
        case 'pret'
            pret = varargin{option+1};
        case 'postt'
            postt = varargin{option+1};
        case 'bin'
            bin = varargin{option+1};
        case 'tmax'
            tmax = varargin{option+1};
        case 'edgespikes'
            edgespikes = varargin{option+1};
        case 'sw1'
            sw1 = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

day = ind(1)
epoch = ind(2)

% Get riptimes
% -------------
if isempty(tetfilter)
    riptimes = sj_getripples_tetinfo(ind, ripples, tetinfo, 'tetfilter', '(isequal($descrip, ''riptet''))','minstd',minstd,'minrip',minrip);
else
    riptimes = sj_getripples_tetinfo(ind, ripples, tetinfo, 'tetfilter', tetfilter, 'minstd', minstd);
end
% Can opt to have a cellcountthresh for each event as in getpopulationevents2 or  sj_HPexpt_ripalign_singlecell_getrip4
% Not using as of now

% Get triggers as rip starttimes separated by at least 1 sec
% ----------------------------------------------------------
rip_starttime = 1000*riptimes(:,1);  % in ms

% Find ripples separated by atleast a second
% --------------------------------------------
iri = diff(rip_starttime);
keepidx = [1;find(iri>=1000)+1];
rip_starttime = rip_starttime(keepidx);

% Implement speed criterion - Keep in ver4. Try both
% ----------------------------------------
if dospeed
    absvel = abs(pos{day}{epoch}.data(:,5)); % Can also use field 9
    postime = pos{day}{epoch}.data(:,1); % in secs
    pidx = lookup(rip_starttime,postime*1000);
    speed_atrip = absvel(pidx);
    lowsp_idx = find(speed_atrip <= lowsp_thrs);
    highsp_idx = find(speed_atrip > highsp_thrs);
    
    rip_starttime = rip_starttime(lowsp_idx);
end


% Get total time for data
totaleptime = diff(spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.timerange)./10000; % in secs
excltime = sum(diff(excludetimes'));
T = totaleptime - excltime;


% Get the spike times for the cross correlation
% --------------------

try
    t1 = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data(:,1);
    t2 = spikes{ind(1)}{ind(2)}{ind(5)}{ind(6)}.data(:,1);
catch
    % if either of those produced an error (meaning no spikes in epoch), we return NaN
    out.index = ind;
    out.T = T;
    out.p1p2 = NaN;
    out.corrstruct = NaN;
    out.rawcorr = NaN;
    out.rawcorr_sm = NaN;
    out.normcorr = NaN;
    out.normcorr_sm = NaN;
    out.Neventscorr = 0;
    out.rip_starttime = rip_starttime;  
    out.pret = pret;
    out.postt = postt;
end


% Convert rip_starttime to secs
rip_starttime = rip_starttime./1000;


% Get the spike times in ripple period for both cells
% ----------------------------------------------------------
t1inc=[]; t2inc=[]; cntrip=0;

for i=2:length(rip_starttime)-1  % Skip first and last
    i;
    % Align to ripples
    % ------------------
    cntrip=cntrip+1;
    currrip = rip_starttime(i);
    %ripsize_cell(cntrip) = rip_sizes(i); ripsize_multi(cntrip) = rip_sizes(i);
    
    currspks =  t1(find( (t1>=(currrip-pret)) & (t1<=(currrip+postt)) ));
    t1inc = [t1inc; currspks];
    
    currspks =  t2(find( (t2>=(currrip-pret)) & (t2<=(currrip+postt)) ));
    t2inc = [t2inc; currspks];

end



% Get raw-cross-correlation
% --------------------------

if ~isempty(t1inc) && ~isempty(t2inc)
    
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
    
   
    % Expected probability
    p1 = xc.nspikes1/T; p2 = xc.nspikes2/T; % Fir rate in Hz
    exp_p = p1*p2; % per sec
    
%     % Crosscov
%     crosscov = (xc.c1vsc2 ./ (bin*T))-exp_p;
%     % Convert to Z-score
%     factor = sqrt((bin*T) / exp_p);
%     Zcrosscov = crosscov .* (factor);
    
    %Raw Corr  Normalize by geometric mean: Units will be coincidences/ spike
    rawcorr = xc.c1vsc2; 
    normcorr = xc.c1vsc2 ./ sqrt(xc.nspikes1 * xc.nspikes2);
    
    % Smooth values as well
    try
        nstd=round(sw1/(xc.time(2) - xc.time(1))); % will be 3 std
        %g1 = gaussian(nstd, 2*nstd+1);    %sigma=3, npoints=7 (+/- 3 bins     
        % Can also do a 3 point gaussian with sigma=3. Almost a boxcar.
        g1 = gaussian(nstd, 2*nstd+1); 
        
    catch
        keyboard;
    end
    rawcorr_sm = smoothvect(rawcorr, g1);
    normcorr_sm = smoothvect(normcorr, g1);
    
    % Plot to check
    %figure; hold on; plot(xc.time,normcorr,'r'); plot(xc.time, normcorr_sm, 'k','LineWidth',2);
    
    % Get Nevents in raw correlation from -200ms to 200ms
    bins = find(abs(xc.time)<=0.2);
    Neventscorr = sum(xc.c1vsc2(bins));
    
else
    
    out.index = ind;
    out.T = T;
    out.p1p2 = NaN;
    out.corrstruct = NaN;
    out.rawcorr = NaN;
    out.rawcorr_sm = NaN;
    out.normcorr = NaN;
    out.normcorr_sm = NaN;
    out.Neventscorr = 0;
    out.rip_starttime = rip_starttime;
    out.pret = pret;
    out.postt = postt;
    return;
end


% Output
% -------
out.index = ind;
% Corrln
% ------
out.T = T;
out.p1p2 = exp_p;
out.corrstruct = xc;
out.rawcorr = rawcorr;
out.rawcorr_sm = rawcorr_sm;
out.normcorr = normcorr;
out.normcorr_sm = normcorr_sm;
out.Neventscorr = Neventscorr;
% Rip
% ---
out.rip_starttime = rip_starttime;
out.pret = pret;
out.postt = postt;












