

function out = dfa_lickphaseSUclustering(index, excludeperiods, varargin)

warning('unfinished')

fprintf('%d %d %d %d\n',index)
reqData = {'lick', 'spikes', 'cellinfo'};
for s = 1:length(reqData)
    if ~any(cellfun(@(x) strcmp(x,reqData{s}), varargin(1:2:end), 'un', 1))
        error('missing required data')
    end
end

pconf = paramconfig;
eventName = 'lick';
TF = 1; % legacy tetfilter
% plotfigs = 1;
% displayplots = 0;
% saveplots = 1;
% savefigas = {'png', 'mfig'};
if ~isempty(varargin)
    assign(varargin{:});
end

out.index = index;
out.dims = {};
out.data = [];
out.mmVmags = [];
out.mmVangs = [];

day = index(1);
epoch = index(2);
ntrode = index(3);
cluster = index(4);
try
    spikeTimes = spikes{day}{epoch}{ntrode}{cluster}.data(:,1);
catch
    return
end

% get events
evvar = varargin{find(cellfun(@(x) ~isempty(x), strfind(varargin(1:2:end), ...
    eventName), 'un', 1))*2-1};
events = eval(evvar);
try
    lickTimes = events{day}{epoch}.starttime;
    id = events{day}{epoch}.id;
catch
    lickTimes = events{day}{epoch}{TF}.starttime;
    id = events{day}{epoch}{TF}.id;
end
% timefilter events
includetimes = ~isExcluded(lickTimes, excludeperiods); % result: 1 include, 0 exclude
lickTimes = lickTimes(includetimes);

% get the lick bout intervals
% criteria: gaps > 500 ms, bouts with at least 8 licks? how many are in error?
% get distribution of licks in unrewarded and rewarded trials

% for each lick bout, histc the spikes, for each spike in inter lick interval bin 
% get the proportion of the spike time from last lick of that spike's ILI. 
% spikeangles = (Lick1 - spiketime) / (Lick2 - Lick1) …
% Then do abs(mean(exp(1i*spikeradians))) to get the mean mrv magnitude
%%
lickgap = .5;
boutIntvStart = lickTimes(find(diff([diff(lickTimes) < lickgap]) == 1)+1);
boutIntvEnd = lickTimes(find(diff([diff(lickTimes) > lickgap]) == 1)+1);

while boutIntvStart(1) > boutIntvEnd(1)
    boutIntvEnd(1) = [];
    boutIntvStart(end) = [];
end
boutNum = 10;
% filter bouts with less than boutNum
licksInbout = logical(isExcluded(lickTimes, [boutIntvStart boutIntvEnd]));
N = histcounts(lickTimes(licksInbout), sort([boutIntvStart; boutIntvEnd]));
incintv = N(1:2:end) > boutNum;
boutIntvStart = boutIntvStart(incintv);
boutIntvEnd = boutIntvEnd(incintv);

%% get DIO reward output
isinput = cellfun(@(x) isequal(x.input,1), DIO{day}{epoch}, 'un', 1);
dioID = cellfun(@(x) str2double(regexp(x.original_id,'\d*','Match')), ...
    DIO{day}{epoch}, 'un', 1);
outputdios = task{day}{epoch}.outputdio;
outdioIdx = find(all([ismember(dioID, outputdios)' ~isinput'],2));
dDF = [];
for c = 1:length(outdioIdx)
    ch = outdioIdx(c);
    dioTimeCh = double(DIO{day}{epoch}{ch}.times);
    dioValsCh = double(DIO{day}{epoch}{ch}.values);
    % dios are Up (1) to Down (0)
    while ~isempty(dioValsCh) && dioValsCh(1) == 0
        dioValsCh(1) = []; dioTimeCh(1) = []; 
    end
    while ~isempty(dioValsCh) && dioValsCh(end) == 1
        dioValsCh(end) = []; dioTimeCh(end) = [];
    end
    if ~isempty(dioValsCh)
        dioVd = [1; find(abs(diff(dioValsCh)))+1];
        diotimes = dioTimeCh(dioVd);
        dioStEnd = [diotimes(1:2:end) diotimes(2:2:end)];
        isout = any(ismember(outputdios, ch));
        dDF = [dDF; ch*ones(length(dioStEnd(:,1)),1) dioStEnd ...
            isout*ones(length(dioStEnd(:,1)),1)];
    else
        continue
    end
end

figure(2)
%%
cla
scatsize = 2000;
subplot(2,1,1)
% plot spikes
scatter(spikeTimes, .9*ones(size(spikeTimes)), scatsize, 'r', 'Marker', '+');
hold on
% plot lick times
scatter(lickTimes, ones(size(lickTimes)), scatsize, 'b', 'Marker', '+')
% plot valid lick bout start and ends
scatter(boutIntvStart, 1.1*ones(size(boutIntvStart)), scatsize, 'g', 'Marker', '+')
scatter(boutIntvEnd, 1.1*ones(size(boutIntvEnd)), scatsize, 'm', 'Marker', '+')
% plot reward outputs
if ~isempty(dDF)
    l = scatter(dDF(:,2),1.1*ones(size(dDF,1),1), scatsize, 'o','Marker', '+');
end
ylim([.8 1.2])
legend({'spikes', 'licks', 'boutStart', 'boutEnd', 'reward'})
hold off
% ylim([-1 1])
% xlim([1637 1648])

%%
subplot(2,2,3)
boutIntv = [boutIntvStart boutIntvEnd];
spikesinbouts = spikeTimes(logical(isExcluded(spikeTimes, boutIntv(1,:))));
[N,~,spbin] = histcounts(spikesinbouts, lickTimes);
relspiketime = spikesinbouts - lickTimes(spbin);
lickdiff = diff(lickTimes);
ili = lickdiff(spbin);
sprads = 2*pi*(relspiketime ./ ili);
meanvec = mean(exp(1i*sprads));
meanMRVmag = abs(meanvec);
vecang = angle(meanvec);
polar([zeros(size(sprads,1),1) sprads]',repmat([0 1],size(sprads,1),1)', 'k')
hold on
polar([0; vecang], [0; meanMRVmag], 'r')
hold off

% out.mmVmags = [];
% out.mmVangs = [];
% for i = 1:size(boutIntv,1)
%     boutlicks = eventTimes(all([eventTimes>=boutIntv(i,1) eventTimes<=boutIntv(i,2)],2));
%     boutspikes = spiketimes(all([spiketimes>boutIntv(i,1) spiketimes<boutIntv(i,2)],2));
% %     if length(boutspikes) < 10
% %         continue
% %     end
%     [N,~,spbin] = histcounts(boutspikes, boutlicks);
%     spikerads = [];
%     for s = 1:length(unique(spbin)) % for each spike, find the proportion of it's ILI
%         ili = boutlicks(spbin(s)+1) - boutlicks(spbin(s));
%         sptime = boutspikes(s) - boutlicks(spbin(s));
%         prpSpL = sptime/ili;
%         spikerads = [spikerads; 2*pi*prpSpL];
%     end
%     meanvec = mean(exp(1i*spikerads));
%     meanMRVmag = abs(meanvec);
%     vecang = angle(meanvec);
%     cla
%     polar([zeros(size(spikerads,1),1) spikerads]',repmat([0 1],size(spikerads,1),1)', 'k')
%     hold on
%     polar([0; vecang], [0; meanMRVmag], 'r')
%     hold off
%     out.mmVmags = [out.mmVmags; meanMRVmag];
%     out.mmVangs = [out.mmVangs; vecang];
% end



end
