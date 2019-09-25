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



epochtype='Run';

%epochfilter{1} = ['isequal($epochtype, ''Run'')'];
epochfilter{1} = ['isequal($epoch, 4)'];

cellfilter = '(isequal($area, ''CA1'') && ($meanrate <10 )  && ($numspikes>100) )'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%timefilter = { {'JY_getlinvelocity', '($velocity <0.07) & ($velocity >-0.07)'} };

timefilter = { {'JY_getriptimes','($nripples == 0)', [], 5,'cellfilter', '(isequal($area, ''CA1''))'},{'JY_getlinvelocity', strcat('$velocity > ',num2str(minVPF))}};

%timefilter = { {'JY_getlinvelocity', '(($velocity) >= 0))', 6} };
trialf = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);

%run function- single cells
%--------------------------------------------
iterator = 'singlecellanal';

trialf = setfilteriterator(trialf,iterator);

trialf=setfilterfunction(trialf, 'JY_calcinfieldtime', {'spikes','linpos'},minV);

% spikes during running
trialf=runfilter(trialf);

figure;


n=6;

posx=data{1,10}{1,2}.Pos.correcteddata(trialf.output{1,1}{1,n}.inplacefieldtimeindices(:,1),2);

posy=data{1,10}{1,2}.Pos.correcteddata(trialf.output{1,1}{1,n}.inplacefieldtimeindices(:,1),3);


plot(data{1,10}{1,2}.Pos.correcteddata(:,2),data{1,10}{1,2}.Pos.correcteddata(:,3),'k');
hold on;
plot(posx,posy,'.r');

xlim([0 160]);
ylim([0 120]);






% 
% % get unique trial spikes
% tifraction=[];
% itfraction=[];
% plottrialspike=[];
% plotintertrialspike=[];
% for i=1:size(trialf.output{1,1},2)
%     trialn=size(trialf.output{1,1}{1,i}.uniquetrialspike(trialf.output{1,1}{1,i}.uniquetrialspike>0),1);
%     itrialn=size(intertrialf.output{1,1}{1,i}.uniqueintertrialspike(intertrialf.output{1,1}{1,i}.uniqueintertrialspike>0),1);
%     if trialn==0 || itrialn==0
%         tifraction(1,i)=0;
%         tifraction(2,i)=0;
%         tifraction(3,i)=0;
%         
%         itfraction(1,i)=0;
%         itfraction(2,i)=0;
%         itfraction(3,i)=0;
%         
%     else
%         % how many are followed by intertrial spikes
%         intertrialn=size(intersect(trialf.output{1,1}{1,i}.uniquetrialspike(trialf.output{1,1}{1,i}.uniquetrialspike>0),...
%         intertrialf.output{1,1}{1,i}.uniqueintertrialspike(intertrialf.output{1,1}{1,i}.uniqueintertrialspike>0)),1);
%         
%         tifraction(1,i)=intertrialn./trialn;
%         tifraction(2,i)=trialn;
%         tifraction(3,i)=intertrialn;
%         
%         
%         
%         
%         itfraction(1,i)=intertrialn./itrialn;
%         itfraction(2,i)=itrialn;
%         itfraction(3,i)=intertrialn;
%         
%         % get positions of all trial spikes
%         trialspikepos=trialf.output{1,1}{1,i}.data(trialf.output{1,1}{1,i}.trialspike>0,2:3);
%         intertrialspikepos=intertrialf.output{1,1}{1,i}.data(intertrialf.output{1,1}{1,i}.intertrialspike>0,2:3);
%         
%         plottrialspike=[plottrialspike;trialspikepos];
%         plotintertrialspike=[plotintertrialspike;intertrialspikepos];
%     end
%     i=i+1;
% end
% figure;
% [tihisty]=hist(tifraction(1,:),0:0.1:1);
% tihisty=tihisty./sum(tihisty);
% bar(0:0.1:1,tihisty,'histc');
% 
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w');
% hold on;
% 
% 
% [ithisty]=hist(itfraction(1,:),0:0.1:1);
% ithisty=ithisty./sum(ithisty);
% bar(0:0.1:1,ithisty,'histc');
% 
% figure;
% 
% bar(0:0.1:1,[tihisty' ithisty'],'histc');
% 
% 
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','b','EdgeColor','w');
% figure;
% plot(plottrialspike(:,1),plottrialspike(:,2),'.r');
% hold on;
% plot(plotintertrialspike(:,1),plotintertrialspike(:,2),'.b');
% %clear all;
