

% Script options.
runscript = 1;

if runscript

%Animal selection
%-----------------------------------------------------
animals = {'Egypt'};
%-----------------------------------------------------

%Filter creation
%--------------------------------------------------------
epochfilter = 'isequal($type, ''run'')';
dayfilter = 1:11;
%epochfilter = 'isequal($type, ''run'')';

ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && isequal($type, ''principal''))';
ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100) && isequal($type, ''principal''))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && isequal($type, ''principal''))';

timefilter = { {'getriptimes', '($nripples == 0)', [], 'cellfilter', '(isequal($area, ''CA1''))'}, ...     % exclude ripples
               {'kk_get2dstate', '($velocity >= 0)'} };                            % velocity cutoff

iterator = 'singlecellanal';

ca1f = createfilter('days',dayfilter,'animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter,'iterator', iterator);
ca2f = createfilter('days',dayfilter,'animal',animals,'epochs',epochfilter,'cells',ca2cellfilter,'excludetime', timefilter,'iterator', iterator);
ca3f = createfilter('days',dayfilter,'animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime', timefilter,'iterator', iterator);

ca1f = setfilterfunction(ca1f, 'dfakk_rawfields', {'spikes', 'pos'});
ca2f = setfilterfunction(ca2f, 'dfakk_rawfields', {'spikes', 'pos'});
ca3f = setfilterfunction(ca3f, 'dfakk_rawfields', {'spikes', 'pos'});

%'ca1'
ca1f = runfilter(ca1f);
%'ca2'
ca2f = runfilter(ca2f);
%'ca3'
ca3f = runfilter(ca3f);
 
end




%% Plot CA1.

if 1
for a=1:length(animals)
    
    % consolidate cells
    indices = []
    for i = 1:length(ca1f(a).output{1})
        indices = [indices ; ca1f(a).output{1}(i).index]
    end
    cellindices = unique(indices(:,[1 3 4]),'rows');
    
for k = 1:size(cellindices,1)
    % find all epochs that cell participates in
        figure
        indices2 = [];
    while rowfind(cellindices(k,:),indices(:,[1 3 4])) ~= 0
        row = rowfind(cellindices(k,:),indices(:,[1 3 4]));
        indices2 = [indices2 ; row];
        indices(row,:) = [nan nan nan nan];
    end
    for kk = 1:size(indices2,1)
        subplot(1,size(indices2,1),kk)
        r = ca1f(a).output{1}(indices2(kk));
        % plot all positions
        plot(r.posdata(:,2),r.posdata(:,3),'Color',[.85 .85 .85],'Linewidth',2);
        hold on
        % plot all spikes
        plot(r.spikes(:,2),r.spikes(:,3),'.k');    
        title(sprintf('animal %d : %d %d %d %d',a,r.index));
    end
end
end
end

%% Plot CA2.

if 1
for a=1:length(animals)
    
    % consolidate cells
    indices = []
    for i = 1:length(ca2f(a).output{1})
        indices = [indices ; ca2f(a).output{1}(i).index]
    end
    cellindices = unique(indices(:,[1 3 4]),'rows');
    
for k = 1:size(cellindices,1)
    % find all epochs that cell participates in
        figure
        indices2 = [];
    while rowfind(cellindices(k,:),indices(:,[1 3 4])) ~= 0
        row = rowfind(cellindices(k,:),indices(:,[1 3 4]));
        indices2 = [indices2 ; row];
        indices(row,:) = [nan nan nan nan];
    end
    for kk = 1:size(indices2,1)
        subplot(1,size(indices2,1),kk)
        r = ca2f(a).output{1}(indices2(kk));
        % plot all positions
        plot(r.posdata(:,2),r.posdata(:,3),'Color',[.85 .85 .85],'Linewidth',2);
        hold on
        % plot all spikes
        plot(r.spikes(:,2),r.spikes(:,3),'.','Color',[.1 .7 .1]);    
        title(sprintf('animal %d : %d %d %d %d',a,r.index));
    end
end
end
end

%% Plot CA3.

if 1
for a=1:length(animals)
    
    % consolidate cells
    indices = []
    for i = 1:length(ca3f(a).output{1})
        indices = [indices ; ca3f(a).output{1}(i).index]
    end
    cellindices = unique(indices(:,[1 3 4]),'rows');
    
for k = 1:size(cellindices,1)
    % find all epochs that cell participates in
        figure
        indices2 = [];
    while rowfind(cellindices(k,:),indices(:,[1 3 4])) ~= 0
        row = rowfind(cellindices(k,:),indices(:,[1 3 4]));
        indices2 = [indices2 ; row];
        indices(row,:) = [nan nan nan nan];
    end
    for kk = 1:size(indices2,1)
        subplot(1,size(indices2,1),kk)
        r = ca3f(a).output{1}(indices2(kk));
        % plot all positions
        plot(r.posdata(:,2),r.posdata(:,3),'Color',[.85 .85 .85],'Linewidth',2);
        hold on
        % plot all spikes
        plot(r.spikes(:,2),r.spikes(:,3),'.r');    
        title(sprintf('animal %d : %d %d %d %d',a,r.index));
    end
end
end
end
