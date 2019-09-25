global subplot_count;
subplot_count = 1;

Veqn = '>=0';
minV =  str2num(Veqn(end));
maxstage = 3; % [1 2 3]
minVPF = 3; %cm/sec
minPeakPF = 3;
lessthan=0;
includestates = 6;

%Animal selection
%-----------------------------------------------------
animals = {'I1'};
%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filterionno

% get all spikes during movement

days='[10]';%,'1:10';
%days = '[1:1]';
%days = '[9:9]';

%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

%epochfilter{1} = [''];

epochtype='Run';

%epochfilter{1} = ['isequal($epochtype, ''Run'')'];
epochfilter{1} = ['isequal($epoch, 2)'];

cellfilter = '(isequal($area, ''CA1'') && ($meanrate <10 )  && ($numspikes>100) )'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%timefilter = { {'JY_getlinvelocity', '($velocity <0.07) & ($velocity >-0.07)'} };

timefilter2 = { {'JY_getriptimes','($nripples > 0)', [], 5,'cellfilter', '(isequal($area, ''CA1''))'},{'JY_getlinvelocity', strcat('$velocity < ', num2str(minVPF))}};
timefilter = { {'JY_getriptimes','($nripples == 0)', [], 5,'cellfilter', '(isequal($area, ''CA1''))'},{'JY_getlinvelocity', strcat('$velocity >= ',num2str(minVPF))}};
%timefilter = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}, {'JY_getlinvelocity', '$velocity <0.01'}};
%timefilter2 = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter2 = { {'getriptimes','($nripples > 0)',[],'minthresh',2,'cellfilter',cellfilter,'tetfilter',[1 2 3 4]}};

%timefilter = { {'JY_getlinvelocity', '(($velocity) >= 0))', 6} };
trialf = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);
intertrialf = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter2);
trialtimes = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter);

%only include cells with placefields
%if minPeakPF>0
%    includecells = calcincludecells(minVPF, minPeakPF, animals, epochfPF, lessthan);
%5    f = excludecellsf(f, includecells);
%end
%-----------------------------------------------------------


%run function- single cells
%--------------------------------------------
iterator = 'singlecellanal';
iterator2 = 'singleepochanal';

trialf = setfilteriterator(trialf,iterator);
intertrialf = setfilteriterator(intertrialf,iterator);
%trialtimes = setfilteriterator(trialtimes,iterator2);

trialf=setfilterfunction(trialf, 'JY_filterspikesbytrial', {'spikes','data'});
intertrialf=setfilterfunction(intertrialf, 'JY_filterspikesbytrial', {'spikes','data'});
%trialtimes=setfilterfunction(trialtimes, 'JY_getreroutetrialtime', {'data'});

% spikes during running
trialf=runfilter(trialf);
% spikes after running
intertrialf=runfilter(intertrialf);
% trial times
%trialtimes=runfilter(trialtimes);

% get unique trial spikes
tifraction=[];
itfraction=[];
plottrialspike=[];
plotintertrialspike=[];
for i=1:size(trialf.output{1,1},2)
    trialn=size(trialf.output{1,1}{1,i}.uniquetrialspike(trialf.output{1,1}{1,i}.uniquetrialspike>0),1);
    itrialn=size(intertrialf.output{1,1}{1,i}.uniqueintertrialspike(intertrialf.output{1,1}{1,i}.uniqueintertrialspike>0),1);
    if trialn==0 || itrialn==0
        tifraction(1,i)=0;
        tifraction(2,i)=0;
        tifraction(3,i)=0;
        
        itfraction(1,i)=0;
        itfraction(2,i)=0;
        itfraction(3,i)=0;
        
    else
        % how many are followed by intertrial spikes
        intertrialn=size(intersect(trialf.output{1,1}{1,i}.uniquetrialspike(trialf.output{1,1}{1,i}.uniquetrialspike>0),...
        intertrialf.output{1,1}{1,i}.uniqueintertrialspike(intertrialf.output{1,1}{1,i}.uniqueintertrialspike>0)),1);
        
        tifraction(1,i)=intertrialn./trialn;
        tifraction(2,i)=trialn;
        tifraction(3,i)=intertrialn;
        
        
        
        
        itfraction(1,i)=intertrialn./itrialn;
        itfraction(2,i)=itrialn;
        itfraction(3,i)=intertrialn;
        
        % get positions of all trial spikes
        trialspikepos=trialf.output{1,1}{1,i}.data(trialf.output{1,1}{1,i}.trialspike>0,2:3);
        intertrialspikepos=intertrialf.output{1,1}{1,i}.data(intertrialf.output{1,1}{1,i}.intertrialspike>0,2:3);
        
        plottrialspike=[plottrialspike;trialspikepos];
        plotintertrialspike=[plotintertrialspike;intertrialspikepos];
    end
    i=i+1;
end
figure;
[tihisty]=hist(tifraction(1,:),0:0.1:1);
tihisty=tihisty./sum(tihisty);
bar(0:0.1:1,tihisty,'histc');

h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w');
hold on;


[ithisty]=hist(itfraction(1,:),0:0.1:1);
ithisty=ithisty./sum(ithisty);
bar(0:0.1:1,ithisty,'histc');

figure;

bar(0:0.1:1,[tihisty' ithisty'],'histc');


h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','w');
figure;
plot(plottrialspike(:,1),plottrialspike(:,2),'.r');
hold on;
plot(plotintertrialspike(:,1),plotintertrialspike(:,2),'.b');
%clear all;
