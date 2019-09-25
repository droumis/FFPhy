

% Script options.
runscript = 1;
savedata = 1;
savedir = '/datatmp/kkay/ProcessedData/';
savefile = [savedir 'PlaceFieldsTrajMaps_allvelocities'];



if runscript == 1

%Animal selection
%-----------------------------------------------------
%animals = {'Bond','Frank','Ten'};
%animals = {'Dudley','Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Miles','Conley','Bond','Frank','Nine','Ten'};
%animals = {'Ten'};
%animals = {'Dudley','Miles','Conley', 'Ten'};
%animals = {'Bond','Frank','Nine'};
%animals = {'Frank'};
animals = {'Chapati','Egypt'};
%-----------------------------------------------------


%Filter creation
%--------------------------------------------------------
epochfilter = 'isequal($type, ''run'')';
%epochfilter{1} = ['isequal($environment, ''TrackB'') && ($dailyexposure == 1)'];
%epochfilter{1} = ['isequal($environment, ''TrackA'') && ($dailyexposure == 1)'];
%for i = 1:14
    %epochfilter{i} = ['isequal($environment, ''TrackA'') | isequal($environment, ''TrackB'')'];
%   epochfilter{i} = ['isequal($environment, ''TrackA'')'];
%nd

cellfilter = '($meanrate < 7)';
ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100)) && ($meanrate < 20)';
ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100)) && ($meanrate < 20)';
ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100)) && ($meanrate < 20)';

timefilter = { {'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'}, ...     % exclude ripples
               {'getlinstate', '(($traj ~= -1) & (abs($velocity) >= 0))', 6} };                            % velocity cutoff


%timefilter = {};

iterator = 'singlecellanal';
%iterator = 'multicellanal';

ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'iterator', iterator);
ca2f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca2cellfilter,'iterator', iterator);
ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'iterator', iterator);

%ca1f = setfilterfunction(ca1f, 'cellprop', {'spikes'});
%ca3f = setfilterfunction(ca3f, 'cellprop', {'spikes'});
%cf = setfilterfunction(cf, 'filtercalclinfields', {'spikes', 'linpos'});
ca1f = setfilterfunction(ca1f, 'filtercalclinfields', {'spikes', 'linpos'});
ca2f = setfilterfunction(ca2f, 'filtercalclinfields', {'spikes', 'linpos'});
ca3f = setfilterfunction(ca3f, 'filtercalclinfields', {'spikes', 'linpos'});


%'ca1'
ca1f = runfilter(ca1f);

%'ca2'
ca2f = runfilter(ca2f);

%'ca3'
ca3f = runfilter(ca3f);
 
end

%% Save data.

if savedata == 1
     save(savefile);
end


%% Plot CA1.

if 0
for a=1:length(animals)
for i = 1:length(ca1f(a).output{1})
    figure
    r = ca1f(a).output{1}(i);
    for s = 1:length(r.trajdata)   % maybe not all 4 trajectories explored
        subplot(4,1,s);
        plot(r.trajdata{s}(:,1), r.trajdata{s}(:,5),'k','Linewidth',4);
        ylim([0 40])
    end
    subplot(4,1,1);
    title(sprintf('animal %d : %d %d %d %d',a,r.index));
    %pause
end
end
end

%% Plot CA2.

if 1
for a=1:length(animals)
for i = 1:length(ca2f(a).output{1})
    figure
    r = ca2f(a).output{1}(i);
    for s = 1:length(r.trajdata)   % maybe not all 4 trajectories explored
        subplot(4,1,s);
        plot(r.trajdata{s}(:,1), r.trajdata{s}(:,5),'g','Linewidth',4);
        ylim([0 40])
    end
    subplot(4,1,1);
    title(sprintf('animal %d : %d %d %d %d',a,r.index));
    %pause
end
end
end

%% Plot CA3.

if 0
for a=1:length(animals)
for i = 1:length(ca3f(a).output{1})
    figure
    r = ca3f(a).output{1}(i);
    for s = 1:length(r.trajdata)   % maybe not all 4 trajectories explored
        subplot(4,1,s);
        plot(r.trajdata{s}(:,1), r.trajdata{s}(:,5),'r','Linewidth',4);
        ylim([0 40])
    end
    subplot(4,1,1);
    title(sprintf('animal %d : %d %d %d %d',a,r.index));
    %pause
end
end

end