


function out = dfa_perripspikingcorr(index, excludeperiods, spikes, eventscons, varargin)
% Calculate per rip corr between the spikes from ntA, ntB
% exclude periods and other filters apply to the rips


%%% Inputs
%   index           [day epoch tetrode cell tetrode cell]
%   excludeperiods    [start end] time ranges to exclude
%   spikes          .data spiketimes
%   pos             .data [time x y ...]

%%% Options
%	'minthresh'-- specifies the minimum threshold of a valid event
%   'window'   -- 1x2 vector specifies the window before and after each included event.
%                   Default is 500 mseconds before and after
%                   event start time.
%   'TF'       -- tetfilter, use to filter for tetrodes for which to find
%                    consensus events

%%% Outputs
%     out.index     [D E T C], gives the identity of the cells for each
%                   column in out.out
%     out.time      vector of binedges used for histc
%     out.psth      Array of 1s and 0s, histc count of spikes around
%                   ripples for each ripple in 1 ms bins

%% init
bin = 0.002; % 2 ms for excess corr, corr of spikes, insta FR
frbin = 0.01; % 10 ms for FR
fprintf(sprintf('%d %d %d %d \n',index))
TF = 1;
an = ''; %animal as varargin from iterator
tmax = 1;
win = [0.5 0.5]; % [duration_before duration_after] in sec
sw1 = 0.005;
sw2 = 0.250;
rmstmax = 0.1;
rmsmincounts = 50;


if (~isempty(varargin))
    assign(varargin{:});
end
win = abs(win); % in case any durations are neg

day = index(1);
ep = index(2);
ntA = index(3);
ntB = index(4);
ec = eventscons{day}{ep}{TF};

wininds = win*ec.samprate;
%% process times
% if all i need the eventkons for is for times, i think i can just get
% times from the spikes data struct right? just stay with it for now

try
    times = ec.eegtimesvec_ref; %use lfp times in ec struct if it exists
catch
    times = ec.timerange(1):1/ec.samprate:ec.timerange(end);
end

% ones = included, zeros = excluded
includetimes = ~isExcluded(times, excludeperiods);
includetimes = includetimes(:);

if sum(includetimes)==0
    error('no includetimes %s d%de%d \n', an, day, ep);
end
% if isempty(ec)
%     error('empty eventcons %s d%de%d \n', an, day, ep);
% end

%% process events
% reconstitute the cons events
% use this instead of the eventcons directly because it's already been
% filtered to certain ripples according to the time filter options
eventstart = times([0 diff(includetimes')]==1);
eventend = times(diff(includetimes')==-1);

% throw out last event if interrupted by end of epoch
if (eventend(end)-eventstart(end) < 0)
    eventstart = eventstart(1:(end-1));
end
% throw out first event if occurred before beginning of epoch
if (eventend(1)-eventstart(1) < 0)
    eventend = eventend(2:end);
end
% throw out any event that begins less than window(2) (i.e. 0.5 seconds) from end of epoch
while (eventend(end) > times(end)-win(2))
    eventend(end) = [];
    eventstart(end) = [];
end

% choose event time to be the START of the event period
eventtimes = eventstart';

% % get indices of each event into the times
% eventStartIndices = lookup(eventstart', times);
% eventEndIndices = lookup(eventend', times);

%% Collect rasters and calculate per rip corr between ntA ntB

% histc bins, where center bin is centered time 0
time = (-win(1)-0.5*bin):bin:(win(2)+0.5*bin);
frtime = (-win(1)-0.5*frbin):frbin:(win(2)+0.5*frbin);
% retrieve spikes
try
    % combine spikes from all clusters on each tetrode
    % since these are all spike times, just collect them all and sort
    ntA_spiketimes = sort(cell2mat(cellfun(@(x) x.data, ...
        spikes{day}{ep}{ntA}, 'un',0)'));
    ntB_spiketimes = sort(cell2mat(cellfun(@(x) x.data, ...
        spikes{day}{ep}{ntB}, 'un',0)'));
catch
    error('couldnt load spikes %s d%de%d ntA:%d ntB:%d \n', an, day, ep);
    out = make_blank();
    return
end
% iterate through each event and calculate corrcoef of spike trains
% and 10 ms binned firing rate and instantaneous firing rate (1ms resolution)

% get indices of each event into the spike times 
ntA_eventstartinds = lookup(eventtimes-win(1), ntA_spiketimes);
ntA_eventendinds =  lookup(eventtimes+win(2), ntA_spiketimes);
ntB_eventstartinds = lookup(eventtimes-win(1), ntB_spiketimes);
ntB_eventendinds = lookup(eventtimes+win(2), ntB_spiketimes);

time = (-win(1)-0.5*bin):bin:(win(2)+0.5*bin);
for r=1:length(eventtimes)
    % there's an error if non unique spiketimes, so only use 1 spike per
    % time. mountainsort allows spikes with identical times -_-
    ntA_spikes_inripwin = unique(ntA_spiketimes(ntA_eventstartinds(r):ntA_eventendinds(r)));
    ntB_spikes_inripwin = unique(ntB_spiketimes(ntB_eventstartinds(r):ntB_eventendinds(r)));
    % rasters: spikes, firing rate, instantaneous firing rate
    ntA_spike_raster(r,:) = histc(ntA_spikes_inripwin , eventtimes(r) + time);
    ntB_spike_raster(r,:) = histc(ntB_spikes_inripwin , eventtimes(r) + time);
    ntA_instaFR_raster(r,:) = instantfr(ntA_spikes_inripwin, [eventtimes(r) + time]);
    ntB_instaFR_raster(r,:) = instantfr(ntB_spikes_inripwin, [eventtimes(r) + time]);
    ntA_10msFR_raster(r,:) = histc(ntA_spikes_inripwin , eventtimes(r) + frtime);
    ntB_10msFR_raster(r,:) = histc(ntB_spikes_inripwin , eventtimes(r) + frtime);
    
    %% xcorr and excess corr per rip
    xc = spikexcorr(ntA_spikes_inripwin, ntB_spikes_inripwin, bin, tmax);
    ntAB_xcorr(r,:) = xc.c1vsc2;
    % compute the excess correlation at 0 lag
	ntAB_excesscorr(r) = excesscorr(xc.time, xc.c1vsc2, xc.nspikes1, ...
        xc.nspikes2, sw1, sw2);
    % compute RMS
	ntAB_rms(r) = xcorrrms(xc.time, xc.c1vsc2, rmstmax, rmsmincounts);
    
    %% mathworks corr per rip
%     [ntAB_spike_xcorr, xcorr_lags] = xcorr(ntA_10msFR_raster(r,:) - ...
%         mean(ntA_10msFR_raster(r,:)), ntB_10msFR_raster(r,:) - ...
%         mean(ntB_10msFR_raster(r,:)), 50, 'coeff');
%     line(xcorr_lags,ntAB_spike_xcorr)

    ntAB_spike_corr(r) = corr(ntA_spike_raster(r,:)', ntB_spike_raster(r,:)');
%     ntAB_spike_corr(r) = diag(tmpcorr); % just keep the time-paired ABrasters' corr 
    ntAB_instaFR_corr(r) = corr(ntA_instaFR_raster(r,:)', ntB_instaFR_raster(r,:)');
%     ntAB_instaFR_corr(r) = diag(tmpcorr);
    ntAB_10msFR_corr(r) = corr(ntA_10msFR_raster(r,:)', ntB_10msFR_raster(r,:)');
%     ntAB_10msFR_corr(r) = diag(tmpcorr);
end

assert(length(ntAB_spike_corr) == length(eventtimes));

%% outputs %%%%%%%%%
out.ntA_spike_raster = ntA_spike_raster;
out.ntB_spike_raster = ntB_spike_raster;
out.ntA_10msFR_raster = ntA_10msFR_raster;
out.ntB_10msFR_raster = ntB_10msFR_raster;
out.ntA_instaFR_raster = ntA_instaFR_raster;
out.ntB_instaFR_raster = ntB_instaFR_raster;

out.ntAB_spike_corr = ntAB_spike_corr;
out.ntAB_instaFR_corr = ntAB_instaFR_corr;
out.ntAB_10msFR_corr = ntAB_10msFR_corr;

out.ntAB_xcorr = ntAB_xcorr;
out.ntAB_excesscorr = ntAB_excesscorr;
out.ntAB_rms = ntAB_rms;

out.time = time;
out.index = index;
out.win = win;
% out.eventStartIndices = eventStartIndices;
% out.eventEndIndices = eventEndIndices;
out.excludeperiods = excludeperiods;

out.xc_time = xc.time;
end

function out = make_blank(varargin)
out.ntA_spike_raster = [];
out.ntB_spike_raster = [];
out.ntA_10msFR_raster = [];
out.ntB_10msFR_raster = [];
out.ntA_instaFR_raster = [];
out.ntB_instaFR_raster = [];

out.ntAB_spike_corr = [];
out.ntAB_instaFR_corr = [];
out.ntAB_10msFR_corr = [];

out.ntAB_xcorr = [];
out.ntAB_excesscorr = [];
out.ntAB_rms = [];

out.time = [];
out.index = index;
out.win = [];

out.excludeperiods = [];

out.xc_time = [];
end