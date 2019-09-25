function [] = PlotRerouteSeg_multiday( animaldir,animalname,daylist, barrier)
% Plots the performance for runs with different segments travelled or
% distance travelled over days
% animaldir     - directory where the output from batchDIOextract is no need
%                 to type full path eg.'/F12_/'
% animalname    - name of the animal, eg 'F1'
% daylist       - vector of day or days to be analysed, eg [1], [1:3]
% barrier       - vector containing epoch and path on which barriers are placed
%                 first number indicates the day,second indicates the epoch for that day
%                 and third specify the intersections between which a
%                 barrier is placed
%                 eg. [2 2 1 2]
% needs connectivity matrix located in
% /home/jai/Src/NSpikeProcess/rerouteconnectivity.mat
% assumes partial segment crossings as 40cm (see segment identity 9 of distance table in connectivity matrix)

% set moving window for running average
movwin=10;
% get the current directory to go back to after the function is done
currentdir = pwd;
% get day by day data
% declare variables
% all visits across days
dayvisit=[];
% length of visit of each day across days
daylength=[0 0];
% stats for each epoch of each day [day epoch run_number]
dayepochstat=[];
% stats for barriers for each run of each day
daybarrierruns=[];
% load connectivity matrix and distance table
load('/home/jai/Src/NSpikeProcess/rerouteconnectivity.mat');
rcon=rerouteconnectivity;
clear rerouteconnectivity;
for d=1:size(daylist,2);
    % open files
    dayt = num2str(daylist(1,d));
    % ---- File loading ----
    % See if day number needs 0
    dsz = '';
    if (daylist(1,d) < 10)
        dsz = '0';
    end
    % Load the mat file
    datadir = '/data14/jai/'; % Specify data directory and load the file
    dayt = num2str(daylist(1,d)); % Converts the day and epoch into text
    sfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'data',dsz,dayt,'.mat');
    load(sfilename);
    posfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'linpos',dsz,dayt,'.mat');
    load(posfilename);
    % find no. epochs
    j =  size(data{1,daylist(1,d)},2);
    % data for all epochs of each day
    epochrunmat=[];
    epochmat=[];
    epochstat=[];
    currdaybarrierruns=[];
    for k=1:size(data{1,daylist(1,d)},2);
        if isempty(data{1,daylist(1,d)}{1,k});
            k=k+1;
         elseif strcmp(data{1,daylist(1,d)}{1,k}.Stats.epochtype,'Sleep')==1;
             k=k+1;
            
        else
            if strcmp(data{1,daylist(1,d)}{1,k}.Stats.animal,animalname)==1; %find all epochs for animal
                epoch = k;
                epocht = num2str(epoch);
                %rewardedevents=data{1,day}{1,epoch}.Events.Wellevents((data{1,day}{1,epoch}.Events.Wellevents(:,3)~=-50),:);
                
                rewardedevents=data{1,daylist(1,d)}{1,epoch}.Events.Welltriggers; % reward events as reward release rather than reward trigger
                if ~isempty(rewardedevents);
                    rewardedevents=rewardedevents((rewardedevents(:,4)>100),:);
                    segments=[]; %clean matrix
                    barrierruns=[];
                    % get start and end times for each run
                    %             for i=1:size(data{1,daylist(1,d)}{1,epoch}.Run,1);
                    %                 tstart=data{1,daylist(1,d)}{1,epoch}.Run(i,3);
                    %                 tend=data{1,daylist(1,d)}{1,epoch}.Run(i,4);
                    %
                    for i=1:size(rewardedevents,1)-1;
                        timeind=rewardedevents(i,1)+rewardedevents(i,4);
                        tstartref=min(find(data{1,daylist(1,d)}{1,epoch}.Pos.correcteddata(:,1)*10000>(timeind))); % for some reason +10000 seems to be accurate rather than no offset
                        timeind2=rewardedevents(i+1,1)+10000; % when rat first triggers reward well, location is not within defined well region, so take location at 10000 later
                        tendref=max(find(data{1,daylist(1,d)}{1,epoch}.Pos.correcteddata(:,1)*10000<(timeind2)));
                        
                        
                        % get segment info from linpos file
                        %                 startind=find(linpos{1,day}{1,epoch}.statematrix.time(:,1)*10000>tstart);
                        %                 tstartref=min(startind);
                        %                 endind=find(linpos{1,day}{1,epoch}.statematrix.time(:,1)*10000<tend);
                        %                 tendref=max(endind);
                        % get segment info between times
                        seg=linpos{1,daylist(1,d)}{1,epoch}.statematrix.segmentIndex(tstartref:tendref,1);
                        % get exit times during period
                        
                        estartref=find(linpos{1,daylist(1,d)}{1,epoch}.trajmatrix(:,1)*10000>timeind+10000); %removed +10000
                        if ~isempty(estartref);
                            estartref=estartref(1,1);
                            [val eendref]=min(abs(linpos{1,daylist(1,d)}{1,epoch}.trajmatrix(:,1)*10000-timeind2));
                            if ~isempty(estartref);
                                exitref=linpos{1,daylist(1,d)}{1,epoch}.trajmatrix(estartref(1,1):eendref,:);
                            else i=size(data{1,daylist(1,d)}{1,epoch}.Run,1);
                            end
                            % concatenate columns to generate integer representation of
                            % path
                            for q=1:size(exitref,1);
                                exitref(q,4)=str2num(sprintf('%d%d',exitref(q,2),exitref(q,3)));
                            end
                            % remove segments containing barrier
%                             for b=1:size(barrier,1);
%                                 if barrier(b,1)==daylist(1,d);
%                                     if barrier(b,2)==epoch;
%                                         fdirect=str2num(sprintf('%d%d',barrier(b,3),barrier(b,3)));
%                                         forward=find(exitref(:,4)~=fdirect);
%                                         exitref=exitref(forward,:);
%                                         rdirect=str2num(sprintf('%d%d',barrier(b,4),barrier(b,4)));
%                                         reverse=find(exitref(:,4)~=rdirect);
%                                         exitref=exitref(reverse,:);
%                                     end
%                                 end
%                                 b=b+1;
%                             end
                            
                            % find segment number and distance from connectivity matrix rcon
                            for q=1:size(exitref,1);
                                exitref(q,5)=rcon(exitref(q,2),exitref(q,3));
                                exitref(q,6)=disttb(exitref(q,5),1);
                            end
                            % remove incomplete path crossings
                            exitref=exitref(find(exitref(:,2)-exitref(:,3)~=0),:);
                            segments(i,1)=size(exitref,1); % put values in matrix
                            segments(i,2)=sum(exitref(:,6)); % distance for path
                            
                            % get barrier index for each run
                            % calculate run by run if barrier is on
                            
                            barrierevents=data{1,daylist(1,d)}{1,epoch}.Events.Barrier;
                            if ~isempty(barrierevents);
                                
                                % get run start and run end times
                                TS=data{1,daylist(1,d)}{1,epoch}.Pos.correcteddata(tstartref,1)*10000;
                                TE=data{1,daylist(1,d)}{1,epoch}.Pos.correcteddata(tendref,1)*10000;
                                BS=barrierevents(:,1);
                                BE=barrierevents(:,1)+barrierevents(:,4);
                                barriertest=(TS >= BS & TS < BE) | (TE <= BE & TE > BS) | (TS <= BS & TE >= BE);
                                if sum(barriertest)>0;
                                    segments(i,3)=1;
                                else segments(i,3)=0;
                                end
                            else segments(i,3)=0;
                            end
                            
                            
                            if segments(i,1)==0;
                                display(sprintf('Respecify well location for day %d epoch %d reward index %d start at %d end at %d',daylist(1,d),epoch,i,tstartref,tendref));
                            end
                            i=i+1;
                        end
                        
                    end
                else k=k+1;
                end
                
                %currbarrierruns=data{1,daylist(1,d)}{1,epoch}.Run(:,6);
                %barrierruns=[barrierruns;currbarrierruns];
                epochrunmat=[epochrunmat;size(data{1,daylist(1,d)}{1,epoch}.Run,1)];
                epochmat=[epochmat;segments];
                
                % matrix with [day epoch run_number]
                currepochstat=[];
                currepochstat(1,1)=daylist(1,d);
                currepochstat(1,2)=epoch;
                currepochstat(1,3)=size(data{1,daylist(1,d)}{1,epoch}.Run,1);
                epochstat=[epochstat;currepochstat];
                %currdaybarrierruns=[currdaybarrierruns;barrierruns];
                k=k+1;
            else display(sprintf('No %s data from day %s epoch %s ',animalname,num2str(daylist(1,d)),num2str(k)));
                k=k+1;
            end
            
        end
        
    end
    %currday=size(epochmat,1);
    if isempty(epochrunmat);
        d=d+1;
    else
        currday(1,1)=sum(epochrunmat,1); % number of runs for the day
        currday(1,2)=sum(epochmat(:,2)); % distance for the day
        daylength=[daylength;currday]; % number of runs for each day
        dayvisit=[dayvisit;epochmat]; % number of segments for each run over days
        dayepochstat=[dayepochstat;epochstat];
        daybarrierruns=[daybarrierruns;currdaybarrierruns];
        d=d+1;
    end
end


% mark each day with a different background colour
daycumsum=cumsum(daylength,1);
col=[];
o=1;
while o<size(daycumsum,1);
    currcol=[1 1 1;0.75 0.75 0.75];
    col=[col;currcol];
    o=o+2;
end

% figure;
%
% title(sprintf('Run distance for %s day %s to %s', animalname,num2str(min(daylist)),num2str(max(daylist))));
% hold on;
% bar(dayvisit);
% xlabel('Run');
% ylabel('Segments');
% ylim([0 10]);
% xlim([0 size(dayvisit,1)]);


% new figure for moving average

figure;
o=1;
while o<size(daycumsum,1);
    currcol=[1 1 1;0.75 0.75 0.75];
    col=[col;currcol];
    o=o+2;
end

% add day lines and day names

daynames=unique(dayepochstat(:,1),'first');

for n=1:size(daycumsum,1)-1;
    linen=line([daycumsum(n,1) daycumsum(n+1,1)],[0 0],'Color',col(n,:),'LineWidth',5000);
    text(((daycumsum(n+1,1)-daycumsum(n,1))/2+daycumsum(n,1)),1.2,num2str(daynames(n,1)),'HorizontalAlignment','center','FontSize',10);
    set(get(get(linen,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off'); % Exclude line from legend
    
    n=n+1;
end

% add line to show barrier

% for r=1:size(barrier,1);
%     dayind= find(dayepochstat(:,1)==barrier(r,1));
%     epochind= find(dayepochstat(dayind,2)==barrier(r,2));
%     xind=sum(dayepochstat(1:dayind(epochind,1)-1,3));
%     xind2=sum(dayepochstat(1:dayind(epochind,1),3));
%
%     lineb=line([xind xind],[0 1.1],'LineStyle','--', 'Color','k','LineWidth',2);
%     set(get(get(lineb,'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude line from legend
%     hold on;
%     linec=line([xind2 xind2],[0 1.1],'LineStyle','--', 'Color',[0.4 0.4 0.4],'LineWidth',2);
%     set(get(get(linec,'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude line from legend
%
% end

% add line to show barrier trials

barrieron=find(dayvisit(:,3)==1);
for i=1:size(barrieron,1);
    
    lineb=line([barrieron(i) barrieron(i)],[1 1],'LineStyle','.', 'Color','k','LineWidth',10);
    set(get(get(lineb,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off'); % Exclude line from legend
    hold on;
end

hold on;

daycumsum=cumsum(daylength,1);
% get all segment runs
colors=[];
% make a skewed colormap
col=colormap(hsv(15));

for c=16:max(dayvisit(:,1));
    col(c,:)=col(15,:);
end

for v=min(dayvisit(:,1)):max(dayvisit(:,1));
    dayvisitones=zeros(size(dayvisit(:,1)));
    dayvisitones((dayvisit(:,1)==v),1)=1;
    step=movwin;
    avgrunmat = moving(dayvisitones,step);
    xlabel('runs');
    ylabel('Proportion of moving window');
    plot(avgrunmat,'Color',col(v,:),'LineWidth',3);
    
    % exlude line if segment travelled is greater than 20
    
    if v>20;
        set(get(get(linen,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off'); % Exclude line from legend
    end
    
    currcol={};
    currcol{1,1}=num2str(v);
    title({sprintf('Runs with n segments for %s days %s to %s',animalname,num2str(min(daynames)),num2str(max(daynames)));...
        sprintf('%s run moving average',num2str(step))});
    ylim([0 1.3]);
    set(gca,'ytick',[0:0.2:1]);
    xlim([0 size(dayvisitones,1)]);
    hold on;
    colors=[colors;currcol];
    
    if size(colors,1)>20;
        colors=colors(1:20,:);
    end
    v=v+1;
end
legend(colors','Location','EastOutside');

% ----------------------------------------
% new figure for distance based analysis
%-----------------------------------------

figure;
o=1;
col=[];
while o<size(daycumsum,1);
    currcol=[0.8 0.8 0.8;0.65 0.65 0.65];
    col=[col;currcol];
    o=o+2;
end

% add day lines and day names

daynames=unique(dayepochstat(:,1),'first');

n=0;
for n=1:size(daycumsum,1)-1;
    linen=line([daycumsum(n,1) daycumsum(n+1,1)],[0 0],'Color',col(n,:),'LineWidth',1000);
    text(((daycumsum(n+1,1)-daycumsum(n,1))/2+daycumsum(n,1)),1.2,num2str(daynames(n,1)),'HorizontalAlignment','center','FontSize',10);
    set(get(get(linen,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off'); % Exclude line from legend
    text(((daycumsum(n+1,1)-daycumsum(n,1))/2+daycumsum(n,1)),1.15,[num2str(daylength(n+1,2)/100,3) 'm'],'HorizontalAlignment','center','FontSize',10);
    n=n+1;
end

% add line to show barrier
% r=0;
% for r=1:size(barrier,1);
%     dayind= find(dayepochstat(:,1)==barrier(r,1));
%     epochind= find(dayepochstat(dayind,2)==barrier(r,2));
%     xind=sum(dayepochstat(1:dayind(epochind,1)-1,3));
%     xind2=sum(dayepochstat(1:dayind(epochind,1),3));
%
%     % barrier on
%     lineb=line([xind xind],[0 1.1],'LineStyle','--', 'Color','k','LineWidth',2);
%     set(get(get(lineb,'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude line from legend
%     hold on;
%     % barrier off
%     linec=line([xind2 xind2],[0 1.1],'LineStyle','--', 'Color',[0.4 0.4 0.4],'LineWidth',2);
%     set(get(get(linec,'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off'); % Exclude line from legend
%
% end

% add line to show barrier trials

barrieron=find(dayvisit(:,3)==1);
for i=1:size(barrieron,1);
    
    lineb=line([barrieron(i) barrieron(i)],[1 1],'LineStyle','.', 'Color','k','LineWidth',10);
    set(get(get(lineb,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off'); % Exclude line from legend
    hold on;
end


hold on;
v=0;
daycumsum=cumsum(daylength,1);
% get all segment runs
colors=[];
distlist=unique(dayvisit(:,2));
distkey=[min(distlist):(max(distlist)-min(distlist))/8:max(distlist)];

for v=1:size(distlist,1);
    dayvisitones=zeros(size(dayvisit(:,1)));
    dayvisitones((dayvisit(:,2)==distlist(v,1)),1)=1;
    step=movwin;
    avgrunmat = moving(dayvisitones,step);
    xlabel('runs');
    ylabel('Proportion of moving window');
    col=[];
    %col=colormap(hsv(size(distlist,1)));
    % make a skewed colormap
    col=colormap(hsv(12));
    
    for c=13:size(distlist(:,1));
        col(c,:)=col(12,:);
    end
    
    
    plot(avgrunmat,'Color',col(v,:),'LineWidth',2);
    currcol={};
    currcol{1,1}=num2str(distlist(v,1));
    title({sprintf('Distance (cm) of runs for %s days %s to %s',animalname,num2str(min(daynames)),num2str(max(daynames)));...
        sprintf('%s run moving average',num2str(step))});
    ylim([0 1.3]);
    set(gca,'ytick',[0:0.2:1]);
    xlim([0 size(dayvisitones,1)]);
    hold on;
    colors=[colors;currcol];
end
%legend(colors','Location','EastOutside');
step=ceil(size(distlist,1)/9);
colorbar('YTickLabel', num2str(distlist(2:step:end,1)/100,2), 'Location','EastOutside');


end


