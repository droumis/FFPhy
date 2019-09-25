% plots the combined normalised place fields of all places cells active in
% a rippple



Veqn = '>=0'
minV =  str2num(Veqn(end))
maxstage = 3% [1 2 3]
minVPF = 2 %cm/sec
minPeakPF = 3
lessthan=0
includestates = 6

%Animal selection
%-----------------------------------------------------
%animals = {'L2'};

%animals = {'M2','M1','M3','K3','L2','L3','N1'};
animals = {'K3','L3','M2', 'N3', 'N1','L2','M3','M1','P2'};

%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filter

days='[1:8]';%,'1:10';



%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

epochfilter{1} = ['isequal($epochtype, ''Run'')'];
%epochfilter{1} = ['isequal($epoch, 2)'];
%cellfilter = '(isequal($area, ''CA1'') && ($meanrate <7))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%timefilter = { {'JY_getlinvelocity', '$velocity <0'} };
%timefilter = { {'JY_getriptimes','($nripples ==0)', [], 3,'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}, {'JY_getlinvelocity', '$velocity <0.01'}};
%timefilter2 = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter2 = { {'getriptimes','($nripples > 0)',[],'minthresh',2,'cellfilter',cellfilter,'tetfilter',[1 2 3 4]}};

%timefilter = { {'JY_getlinvelocity', '(($velocity) >= 0))', 6} };
%f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter);
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter);


%only include cells with placefields
%if minPeakPF>0
%    includecells = calcincludecells(minVPF, minPeakPF, animals, epochfPF, lessthan);
%5    f = excludecellsf(f, includecells);
%end
%-----------------------------------------------------------


%run function- single cells
%--------------------------------------------
iterator = 'JY_singleepochanal';

f = setfilteriterator(f,iterator);

outall=setfilterfunction(f, 'JY_getreroutetrialvelocity', {'data','linpos'});

outall=runfilter(outall);
day_velocity=[];
meanvelocity=[];
stdvelocity=[];
for mm=1:size(outall,2)
    out=outall(1,mm);
    uniquedays=unique(out.epochs{1,1}(:,1));
    outdata=out.output{1,1};
    % plot trajectory distance
    
    
    % col=jet(length(out.output{1,1}));
    %
    % for i=1:length(out.output{1,1});
    %     plot(out.output{1,1}{1,i}.trajectorydistance,'.','Color',col(i,:));
    %     hold on;
    % end
    %
    %
    % figure;
    %
    % for i=1:length(out.output{1,1});
    %     X=[1:1:length(out.output{1,1}{1,i}.normdist)]';
    %     Y=out.output{1,1}{1,i}.normdist;
    %     Z=ones(length(out.output{1,1}{1,i}.normdist),1)*(length(out.output{1,1})-i+1);
    %     plot3(X,Z,Y,'.','Color',col(i,:));
    %     %bar3(Y,,'detached');
    %     %display(max(out.output{1,1}{1,i}.normdist))
    %     grid on;
    %     hold on;
    % end
    
    
    % plot trajectory time
    % setup variables
    winsize=11;
    epoch_len = [];
    day_len = [];
    trialsegments_all = [];
    day_totalseg = [];
    numseg_mean = [];
    numseg_std = [];
    day_totaltime = [];
    totaltime_mean = [];
    totaltime_std = [];
    trials_day=[];
    barrier=[];
    
    
    uniquedays=unique(out.epochs{1,1}(:,1));
    
    % loop though data for each day
    for i = 1:size(uniquedays)
        trials_day = [];
        jind=find(out.epochs{1,1}(:,1)==uniquedays(i,1));
        for j = 1:size(jind,1);
            trials_day=[trials_day; out.output{1,1}{1,jind(j)}.velocity];
            % collecting trial segments across epochs and store which trial each day and
            % epoch ends so it can be used to plot later
            %     day_len = [day_len; out.output{1,1}{1,jind(j)}.ntrial];
            %     barrier = [barrier; out.output{1,1}{1,jind(j)}.barrier];
            
        end
        %     totaltime_mean =[totaltime_mean; moving(trials_day,winsize,@mean)];
        %     totaltime_std = [totaltime_std; moving(trials_day,winsize,@std)];
        
        % get distribution of velocities
        
        day_velocity{mm}{i}=trials_day;
        day_velocityplot(i,:)=histc(trials_day,[0:1:60]);
        % normalise
        %day_velocity{mm}{i}=day_velocity{mm}{i}/sum(day_velocity{mm}{i});
        day_velocityplot(i,:)=day_velocityplot(i,:)/sum(day_velocityplot(i,:));
        
        % get one mean for each day
        
        daymeanvelocity(mm,i)=mean(abs(trials_day));
        
        epochs = [];
        
        
        % Collect the time of each trial across all epochs in the day
        
        %     totaltime_mean = [totaltime_mean movingstat(trial_time',winsize,@mean)'];
        %     totaltime_std = [totaltime_std movingstat(trial_time',winsize,@std)'];
        
        %trialsegments_all = [trialsegments_all trialsegments_day];
        
    end
    % figure;
    %
    
    meanvelocity(mm,:)=mean(day_velocityplot,1);
    stdvelocity(mm,:)=std(day_velocityplot,1);
    %semilogy(day_velocityplot','.');
    
    
    
    
    % epochsummary=[out.epochs{1,1}(:,1) day_len];
    %
    % % Plotting time per trial
    % figure;
    % set(gcf,'position',[0 0 800 250]);
    % set(gcf,'PaperPositionMode','auto');
    % set(gca,'FontSize',12);
    %
    % % make grey/white strips
    % o=1;
    % col=[];
    %
    % while o<size(uniquedays,1)+1;
    %     currcol=[1 1 1;0.8 0.8 0.8];
    %     col=[col;currcol];
    %     o=o+2;
    % end
    %
    % daycumsum=[];
    % for i=1:size(uniquedays,1)
    %     daycumsum=[daycumsum; sum(epochsummary(epochsummary(:,1)==uniquedays(i,1),2))];
    % end
    %
    % daycumsum=[0;cumsum(daycumsum)];
    %
    % n=0;
    %
    % for n=1:size(daycumsum,1)-1;
    %     linen=line([daycumsum(n,1) daycumsum(n+1,1)],[0 0],'Color',col(n,:),'LineWidth',1000);
    %     text(((daycumsum(n+1,1)-daycumsum(n,1))/2+daycumsum(n,1)),18,num2str(uniquedays(n,1)),'HorizontalAlignment','center','FontSize',12);
    %     set(get(get(linen,'Annotation'),'LegendInformation'),...
    %         'IconDisplayStyle','off'); % Exclude line from legend
    % %     text(((daycumsum(n+1,1)-daycumsum(n,1))/2+daycumsum(n,1)),1.15,[num2str(daylength(n+1,2)/100,3) 'm'],'HorizontalAlignment','center','FontSize',10);
    %     n=n+1;
    % end
    % hold on;
    % % add line to show barrier trials
    %
    % barrieron=find(barrier(:,1==1));
    % for i=1:size(barrieron,1);
    %     lineb=line([barrieron(i) barrieron(i)],[8 8],'LineStyle','s', 'Color','g','LineWidth',2);
    %     set(get(get(lineb,'Annotation'),'LegendInformation'),...
    %         'IconDisplayStyle','off'); % Exclude line from legend
    %     hold on;
    % end
    %
    % plot(totaltime_mean,'LineWidth',2);
    %
    % plot(totaltime_mean+totaltime_std,'r:','LineWidth',2);
    % plot(totaltime_mean-totaltime_std,'r:','LineWidth',2);
    %
    % % annotating
    %
    % % for ii = 1:length(day_len)
    % %     plot([sum(day_len(1:ii)) sum(day_len(1:ii))]-(winsize-1)*(ii), [0 200], 'k');
    % %     %text(sum(day_len(1:ii))-(winsize)*ii, 85, sprintf('Day %d\nn=%d\ntime=%d',days(ii),day_len(ii),floor(day_totaltime(ii))), ...
    % %         %'HorizontalAlignment','right','FontSize',10);
    % %
    % % end
    %
    % xlim([0 max(daycumsum(:,1))]);
    % %ylim([0 max(totaltime_mean+totaltime_std)]);
    % ylim([0 50]);
    
    % title(sprintf('Mean trial speed for %s',out.animal{1,3}));
    % xlabel('Trial');
    % ylabel('Mean speed (cm/s)');
    
    
    
    % % saving
    % [s,mess,messid] = mkdir(sprintf('%sPlot/behav/',directoryname));
    % print(sprintf('%sPlot/behav/%s_trialtimeavg_d%d-%d_w%d',directoryname,fileprefix,days(1),days(end),winsize),'-depsc');
end

% % plot distribution of velocity
% 
% errorbar(repmat(1:1:61,7,1)',meanvelocity',stdvelocity')
% set(gca,'yscale','log')
% legend(animals,'location','eastoutside')
% 
% figure;
% hold on;
% errorbar(repmat(1:1:61,4,1)',meanvelocity([2 3 5 7],:)',stdvelocity([2 3 5 7],:)','r')
% errorbar(repmat(1:1:61,3,1)',meanvelocity([1 4 6],:)',stdvelocity([1 4 6],:)','b')
% 
% 
% set(gca,'yscale','log')
% 
% xlabel('Velocity cm/s')
% title('Normalised mean velocity distribution for all days +/-SD')
% xlim([0 65])
% ylim([0.000001 1])
% 
% 
% meanvelocity(7,:)=meanvelocity(7,:)/2;
% 
% figure;
% plot(meanvelocity')
% set(gca,'yscale','log')
% legend(animals,'location','eastoutside')
% 
% figure;
% plot(meanvelocity([1 4 6],:)','b')
% set(gca,'yscale','log')
% hold on;
% plot(meanvelocity([2 3 5 7],:)','r')
% 
% legend(animals,'location','eastoutside')

%% plot distribution of velocities

% day 1:4

control=[];
for ii=1:4
   for jj=1:4
       control(ii,:)=histc(day_velocity{ii}{jj},[0:5:100]);
       control(ii,:)=control(ii,:)/sum(control(ii,:));
       %control(ii,:)=cumsum(control(ii,:));
       
       %plot(control(ii,:),'b');
       %set(gca, 'YScale','log')
       %hold on;
   end
end

controlm=mean(control,1);
controls=std(control,[],1);

inactivation=[];
for ii=5:9
   for jj=1:4
       inactivation(ii,:)=histc(day_velocity{ii}{jj},[0:5:100]);
       inactivation(ii,:)= inactivation(ii,:)/sum( inactivation(ii,:));
       %inactivation(ii,:)=cumsum(inactivation(ii,:));
       
        %plot(inactivation(ii,:),'r');
        %hold on;
   end
end


inactivationm=mean(inactivation(5:9,:),1);
inactivations=std(inactivation(5:9,:),[],1);
% 
% figure;
% plot(controlm(controlm>0),'b');
% hold on;
% plot((controlm(controlm>0)-controls(controlm>0)),'.b');
% hold on;
% plot((controlm(controlm>0)+controls(controlm>0)),'.b');
% 
% plot(inactivationm(inactivationm>0),'r');
% hold on;
% plot((inactivationm(inactivationm>0)-inactivations(inactivationm>0)),'.r');
% hold on;
% plot((inactivationm(inactivationm>0)+inactivations(inactivationm>0)),'.r');
% 
% controllow=controlm-controls;
% controllow(controllow==0)=0.00000001;
% controllows(controllow<=0)=0.00000000000000001;
% inactivationlow=inactivationm-inactivations;
% inactivationlow(inactivationlow==0)=0.0000001;
% inactivationlow(inactivationlow<=0)=0.00000000000000001;

% figure;
% shadedplot([0:5:100],(inactivationm-inactivations),...
%     (inactivationm+inactivations),[1 0.7 0.7], 'r');
% 
% 
% 
% hold on;
% 
% shadedplot([0:5:100], (controlm-controls),...
%     (controlm+controls),[0.7 0.7 1], 'b');


set(gca, 'YScale','log')

bar([controlm;inactivationm]','grouped');

xlabel('Days');
xlim([0 20]);
set(gca,'XTick',[0:5:25]);
set(gca,'XTickLabel',{'0', '20', '40', '60', '80', '100'});
%set(gca,'YTick',[0:5:20]);

ylim([0 1]);

title('Distribution of velocities days 1-4');
xlabel('Velocity (cm/s)');
ylabel('Proportion of time');
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/velocitydistribution1_4'),'-depsc');
close;
% day 5:8

control=[];
for ii=1:4
   for jj=5:8
       control(ii,:)=histc(day_velocity{ii}{jj},[0:5:100]);
       control(ii,:)=control(ii,:)/sum(control(ii,:));
       %control(ii,:)=cumsum(control(ii,:));
       
       %plot(control(ii,:),'b');
       %set(gca, 'YScale','log')
       %hold on;
   end
end

controlm=mean(control,1);
controls=std(control,[],1);





inactivation=[];
for ii=5:9
   for jj=5:8
       inactivation(ii,:)=histc(day_velocity{ii}{jj},[0:5:100]);
       inactivation(ii,:)= inactivation(ii,:)/sum( inactivation(ii,:));
       %inactivation(ii,:)=cumsum(inactivation(ii,:));
       
        %plot(inactivation(ii,:),'r');
        %hold on;
   end
end


inactivationm=mean(inactivation(5:9,:),1);
inactivations=std(inactivation(5:9,:),[],1);

% figure;
% plot(controlm(controlm>0),'b');
% hold on;
% plot((controlm(controlm>0)-controls(controlm>0)),'.b');
% hold on;
% plot((controlm(controlm>0)+controls(controlm>0)),'.b');
% 
% plot(inactivationm(inactivationm>0),'r');
% hold on;
% plot((inactivationm(inactivationm>0)-inactivations(inactivationm>0)),'.r');
% hold on;
% plot((inactivationm(inactivationm>0)+inactivations(inactivationm>0)),'.r');
% 
% controllow=controlm-controls;
% controllow(controllow==0)=0.00000001;
% controllows(controllow<=0)=0.00000000000000001;
% inactivationlow=inactivationm-inactivations;
% inactivationlow(inactivationlow==0)=0.0000001;
% inactivationlow(inactivationlow<=0)=0.00000000000000001;
% 
% figure;
% shadedplot([0:5:100],(inactivationm-inactivations),...
%     (inactivationm+inactivations),[1 0.7 0.7], 'r');
% 
% 
% 
% hold on;
% 
% shadedplot([0:5:100], (controlm-controls),...
%     (controlm+controls),[0.7 0.7 1], 'b');


%set(gca, 'XScale','log')
figure;
bar([controlm;inactivationm]','grouped');

xlabel('Days');
xlim([0 20]);
set(gca,'XTick',[0:5:25]);
set(gca,'XTickLabel',{'0', '20', '40', '60', '80', '100'});
%set(gca,'YTick',[0:5:20]);

ylim([0 1]);

title('Distribution of velocities days 5-8');
xlabel('Velocity (cm/s)');
ylabel('Proportion of time');
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/velocitydistribution5_8'),'-depsc');
close;




% ----------------------
% plot each box in different colours
% based on solution from http://www.mathworks.com/matlabcentral/answers/22

figure;

datacontrol=[reshape(daymeanvelocity(1:4,1:4),1,16)',  reshape(daymeanvelocity(1:4,5:8),1,16)',];
datacontrolposition=[1 3];
datainactivation=[reshape(daymeanvelocity(5:9,1:4),1,20)', reshape(daymeanvelocity(5:9,5:8),1,20)',];
datainactivationposition=[1.5 3.5];

box_control = boxplot(datacontrol,'colors','b','positions',datacontrolposition,'width',0.3); 
set(gca,'XTickLabel',{' '})

hold on;

box_inactivation = boxplot(datainactivation,'colors','r','positions',datainactivationposition,'width',0.3); 
set(gca,'XTickLabel',{' '})

xlabel('Days');
xlim([0.5 4]);
set(gca,'XTick',[1.25 3.25]);
set(gca,'XTickLabel',{'1-4','5-8'})
set(gca,'YTick',[0:5:20]);
xlim([0.5 4]);
ylim([0 22]);

title('Mean daily velocity');
ylabel('Velocity (cm/s)');
xlabel('Days');


% calculate p values between each group for each day period
p14=ranksum(reshape(daymeanvelocity(1:4,1:4),1,16),reshape(daymeanvelocity(5:9,1:4),1,20));
p58=ranksum(reshape(daymeanvelocity(1:4,5:8),1,16),reshape(daymeanvelocity(5:9,5:8),1,20));


text('Position',[1.25,18],'String',sprintf('p=%s',num2str(p14,2)),'HorizontalAlignment','center');
text('Position',[3.25,18],'String',sprintf('p=%s',num2str(p58,2)),'HorizontalAlignment','center');
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/meandayvelocity'),'-depsc');
close;
