
clear points;
runscript = 1;


if runscript==1
    
for s=1:3
%Animal selection
%-----------------------------------------------------
animals = {'Chapati','Egypt'};
animdir = {'/datatmp/kkay/ChapatiData/','/data12/mari/EgyptData/'}

%Filter creation
%--------------------------------------------------------


% epoch filter
epochfilter{1} = ['isequal($type, ''run'') || isequal($type, ''sleep'')'];

% cell filter
ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100))';
ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100))';
ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100))';

% time filter
if s==1
   str = '((abs($velocity) <= 0.5))';
elseif s==2
   str = '((abs($velocity) >= 8))';
else
   str = '((abs($velocity) > 0.5) & (abs($velocity) < 8))';
end
    
timefilter = { {'get2dstate', str} };

% iterator
iterator = 'singlecellanal';

% filter creation
ca1f{s} = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'cells',ca1cellfilter, 'iterator', iterator);
ca2f{s} = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'cells',ca2cellfilter, 'iterator', iterator);
ca3f{s} = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'cells',ca3cellfilter, 'iterator', iterator);

%ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter, 'iterator', iterator);
%ca2f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca2cellfilter, 'iterator', iterator);
%ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter, 'iterator', iterator);

% set analysis function
ca1f{s} = setfilterfunction(ca1f{s}, 'dfskk_calctotalmeanrate_spikewidth', {'spikes'});
ca2f{s} = setfilterfunction(ca2f{s}, 'dfskk_calctotalmeanrate_spikewidth', {'spikes'});
ca3f{s} = setfilterfunction(ca3f{s}, 'dfskk_calctotalmeanrate_spikewidth', {'spikes'});

% run analysis
ca1f{s} = runfilter(ca1f{s});
ca2f{s} = runfilter(ca2f{s});
ca3f{s} = runfilter(ca3f{s});

end
end


% Plot

%%% Plot all CA1 data from both animals.

caf=ca3f;

alldata=[];
for a=1:length(animals)
    for s=1:length(caf)  % number of states
        for c=1:length(caf{s}(a).output{1})
            % [animal state day epoch tetrode cell meanrate spikewidth] 
            toappend = [a s caf{s}(a).output{1}(c).index caf{s}(a).output{1}(c).rate caf{s}(a).output{1}(c).spikewidth];
            alldata=[alldata ; toappend];       % note appending column corresponding to animal to front
        end
    end
end
fig1 = figure
hold on; dummy1=[]; dummy2=[]; dummy3=[];
for k=1:length(alldata)
    if alldata(k,2)==1
        dummy1=[dummy1 ; alldata(k,7),alldata(k,8)];
    elseif alldata(k,2)==2
        dummy2=[dummy2 ; alldata(k,7),alldata(k,8)];
    elseif alldata(k,2)==3
        dummy3=[dummy3 ; alldata(k,7),alldata(k,8)];
    end
end
        scatter(dummy1(:,1),dummy1(:,2),500,'x','LineWidth',4,'CData',[0 0 0]);
        hold on;
        scatter(dummy2(:,1),dummy2(:,2),500,'x','LineWidth',4,'CData',[1 .5 .5]);
        hold on;
        scatter(dummy3(:,1),dummy3(:,2),500,'x','LineWidth',4,'CData',[.5 .5 1]);        
xlim([0 max(alldata(:,7))])
ylim([0 max(alldata(:,8))])
xlabel('mean firing rate (Hz)')
ylabel('spikewidth')
title('CA3 unit classification','FontSize',24,'FontWeight','bold','FontWeight','bold')

% generate "probelist" variable (nx2) of questionable points for which you want waveforms
        % do this by manually selecting points using Brush/Select Data tool
        % in Matlab figure window

keyboard
points;

% Now plot waveforms for each of these points.

% First find animal-day-epoch-tetrode-cell coordinate of selected points


    fig2=figure
    title('waveforms of selected units')
    allcoord=[];
    
for p=1:size(points,1)
    
    figure(fig2)
    subplot(5,6,p)
    
    ind=rowfind(points(p,:),alldata(:,[7 8]));
    coord=alldata(ind,1:6);
    allcoord=[allcoord ; coord];

    % retrieve spike times in original outputstruct
    for c=1:length(caf(coord(1)).output{coord(2)})
        if sum(caf(coord(1)).output{coord(2)}(c).index==coord(3:6))==4  % unit coordinate matches
            spiketimes=caf(coord(1)).output{coord(2)}(c).spiketimes;
            continue
        end
    end
    % load corresponding spikes .mat file containing waveforms
    dir1=dir(sprintf('%s*%02d',animdir{coord(1)},coord(3)));
    dir2=dir(sprintf('%s%s/%02d-*',animdir{coord(1)},dir1.name,coord(5)));
    dir3=dir(sprintf('%s%s/%s/*%02d-%02d.mat',animdir{coord(1)},dir1.name,dir2(1).name,coord(3),coord(5)));
    load(sprintf('%s%s/%s/%s',animdir{coord(1)},dir1.name,dir2(1).name,dir3.name))
    
    color=[rand rand rand];
    spikewaveforms(spiketimes,waves,timestamps,color);
    title(sprintf('%s state %d // %d %d %d %d',animals{coord(1)},coord(2),coord(3:6)))
    
    % go back and re-plot colored points in original figure
    figure(fig1)
    hold on;
    scatter(alldata(ind,7),alldata(ind,8),250,'filled','o','CData',color)
    
end

    allcoord