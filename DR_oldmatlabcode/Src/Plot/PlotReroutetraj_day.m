function [ output_args ] = PlotReroutetraj_day( animaldir,animalname,day )
% Plots trajectory of each run
% datadirectory - directory where the output from batchDIOextract is
% located, eg '/data14/jai/F12_/'
% animalname -  name of the animal, eg F1
% days - vector of day or days to be analysed, eg [1], [1:3]
%


% get the current directory to go back to after the function is done
currentdir = pwd;

% ---- File loading ----
% See if day number needs 0
dsz = '';
if (day < 10)
    dsz = '0';
end

% Specify data directory and load the file
%animalname = animaldir(1:end-1);
datadir = '/data14/jai/';
%datadir = '/home/daliu/tmp/';

% Converts the day and epoch into text
dayt = num2str(day);

sfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'data',dsz,dayt,'.mat');

% Load the mat file


load(sfilename);

% load posfile
%posfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'_',dsz,dayt,'linpos.mat');
posfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'linpos',dsz,dayt,'.mat');
load(posfilename);

% plot trajectory between times

for k=1:size(data{1,day},2);
    if isempty(data{1,day}{1,k});
        k=k+1;
        
    elseif strcmp(data{1,day}{1,k}.Stats.epochtype,'Sleep')==1;
        k=k+1;
    else
        
        if strcmp(data{1,day}{1,k}.Stats.animal,animalname)==1;
            figure;
            
            epoch = k;
            epocht = num2str(epoch);
            
            % old rewarded event based on well tigger
            rewardedevents=data{1,day}{1,epoch}.Events.Wellevents((data{1,day}{1,epoch}.Events.Wellevents(:,3)~=-50),:);
            
            % new rewarded events based on reward delivery 12/1/2011
            rewardedevents=data{1,day}{1,epoch}.Events.Welltriggers;
            
            barrierevents=data{1,day}{1,epoch}.Events.Barrier;
            
            dist=0;
            nrows=ceil(size(rewardedevents,1)/10);
            columns=10;
            gap_h=0.01;
            gap_w=0.01;
            marg_h=[0.01 0.1];
            marg_w=[0.01 0.01];
            
            
            ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
            
            distance=[];
            
            %             for i=1:size(data{1,day}{1,epoch}.Run,1);
            %                 tstart=data{1,day}{1,epoch}.Run(i,3);
            %                 tend=data{1,day}{1,epoch}.Run(i,4);
            
            % find times rat leaves the intersection zone after a
            % reward is triggered
            
            for i=1:size(rewardedevents,1)-1;
                timeind=rewardedevents(i,1);
                toffset=0;
                tstartref=max(find(data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000<(timeind+toffset)));
                timeind2=rewardedevents(i+1,1);
                tendref=max(find(data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000<(timeind2+toffset)));
                
                if ~isempty(tstartref);
                    % use for maze with no arms
                    %                 [val tstartref]=min(abs((data{1,day}{1,epoch}.Pos.correcteddata(:,1)-tstart)));
                    %                 findind2=find(data{1,day}{1,epoch}.Pos.correcteddata(:,1)<tend);
                    %                 tendref=findind2(end,1);
                    
                    
                    %                 tstartref=find(data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000==tstart);
                    %                 tendref=find(data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000==tend);
                    
                    % get traj between times
                    
                    traj=[data{1,day}{1,epoch}.Pos.correcteddata(tstartref:tendref,2) data{1,day}{1,epoch}.Pos.correcteddata(tstartref:tendref,3)];
                    
                    % calculate angle between consecutive vectors to determine change
                    % in direction
                    
                    for c=2:size(traj,1)-1;
                        %traj(c,1) traj(c,2) - [traj(c-1,1) traj(c-1,2)]
                        x1 = traj(c-1,1);
                        x2= traj(c,1);
                        x3=traj(c,1);
                        x4=traj(c+1,1);
                        y1=traj(c-1,2);
                        y2=traj(c,2);
                        y3=traj(c,2);
                        y4=traj(c+1,2);
                        %v2 = [traj(c+1,1) traj(c+1,2)] - [traj(c,1) traj(c,2)];
                        %if v1+v2~=0;
                        
                        angle = atan2(abs((x2-x1)*(y4-y3)-(x4-x3)*(y2-y1)), ...
                            (x2-x1)*(x4-x3)+(y2-y1)*(y4-y3));
                        traj(c+1,3)=angle*180/pi;
                        
                        c=c+1;
                        
                        %end
                    end
                    
                    traj=[traj data{1,day}{1,epoch}.Pos.correcteddata(tstartref:tendref,5)];
                    
                    % calculate path travelled
                    
                    dist=0;
                    for j=1:size(traj,1)-1;
                        dt=sqrt((traj(j+1,1)-traj(j,1))^2+(traj(j+1,2)-traj(j,2))^2);
                        dist=dist+dt;
                        j=j+1;
                    end
                    
                    distance(i,1)=dist;
                    
                    %subplot(nrows,5,i);
                    
                    axes(ha(i));
                    
                    plot(traj(:,1), traj(:,2));
                    
                    startpos=[traj(1,1) traj(1,2)];
                    endpos=[traj(end,1) traj(end,2)];
                    
                    hold on;
                    scatter(startpos(1,1),startpos(1,2),'s','MarkerFaceColor','m','MarkerEdgeColor','m'); %start marked by magenta
                    hold on;
                    %scatter(endpos(1,1),endpos(1,2),'MarkerFaceColor','c','MarkerEdgeColor','c'); % end marked by cyan
                    
                    
                    hold on;
                    
                    % plot reversal points by looking at angle between points, greater
                    % than 90 degrees
                    
                    revn=find(traj(:,3)>90);
                    rev=traj(revn,:);
                    
                    % add direction file to traj
                    
                    d=1;
                    
                    revn1=[1;revn];
                    
                    while d<=size(revn1,1)-2;
                        traj(revn1(d,1):revn1(d+1,1)-1,4)=1;
                        
                        traj(revn1(d+1):revn1(d+2,1)-1,4)=0;
                        d=d+2;
                    end
                    
                    % plots points of reversal
                    %plot(rev(:,1), rev(:,2),'o','MarkerEdgeColor','g','MarkerSize',5);
                    
                    hold on;
                    
                    rtrajn=find(traj(:,3)==0);
                    rtraj=traj(rtrajn,:);
                    
                    % plots reversal
                    
                    %plot(rtraj(:,1), rtraj(:,2),'.r','MarkerSize',5);
                    
                    % get points of slow movement
                    %                 slown=find(traj(:,3)<3);
                    %                 slow=traj(slown,:);
                    %                 plot(slow(:,1), slow(:,2),'o','MarkerEdgeColor','g','MarkerSize',2);
                    
                    
                    xlimval=[min(data{1,day}{1,epoch}.Pos.correcteddata(:,2)), max(data{1,day}{1,epoch}.Pos.correcteddata(:,2))];
                    ylimval=[min(data{1,day}{1,epoch}.Pos.correcteddata(:,3)), max(data{1,day}{1,epoch}.Pos.correcteddata(:,3))];
                    
                    fac=0.2;
                    
                    xlimval=[xlimval(1)-fac*(xlimval(2)-xlimval(1)), xlimval(2)+fac*(xlimval(2)-xlimval(1))];
                    ylimval=[ylimval(1)-fac*(ylimval(2)-ylimval(1)), ylimval(2)+fac*(ylimval(2)-ylimval(1))];
                    
                    % lable each plot with run number
                    %title(sprintf('Run %d',i),'FontSize',10);
                    ylim(ylimval);
                    xlim(xlimval);
                    %XTickLabel('FontSize',4);
                    
                    distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.85*(ylimval(2)-ylimval(1))+ylimval(1)];
                    titlelabel=[0.6*(xlimval(2)-xlimval(1))+xlimval(1) 0.08*(ylimval(2)-ylimval(1))+ylimval(1)];
                    % distance label
                    text(distlabel(1),distlabel(2),sprintf('%s m',num2str(dist/100,2)),'FontSize',10);
                    text(titlelabel(1),titlelabel(2),sprintf('Run %d',i),'FontSize',10);
                    
                    % barrier label
                    barrierlabel=[0.6*(xlimval(2)-xlimval(1))+xlimval(1) 0.85*(ylimval(2)-ylimval(1))+ylimval(1)];
                    
                    if ~isempty(barrierevents);
                        
                        % get run start and run end times
                        TS=data{1,day}{1,epoch}.Pos.correcteddata(tstartref,1)*10000;
                        TE=data{1,day}{1,epoch}.Pos.correcteddata(tendref,1)*10000;
                        BS=barrierevents(:,1);
                        BE=barrierevents(:,1)+barrierevents(:,4);
                        barriertest=(TS >= BS & TS < BE) | (TE <= BE & TE > BS) | (TS <= BS & TE >= BE);
                        if sum(barriertest)>0;
                            text(barrierlabel(1),barrierlabel(2),'barrier','FontSize',10);
                        end
                    end
                    
                    set(gca,'xtick',[],'ytick',[]);
                    set(gca,'PlotBoxAspectRatio',[1 1 1]);
                    
                    
                    
                    i=i+1;
                end
            end
            
            % Print title for whole figure
            
            ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
                1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
            
            text(0.5, 0.98,sprintf('%s day %s epoch %s square=start',animalname,dayt,epocht),'HorizontalAlignment'...
                ,'center','VerticalAlignment', 'top');
%             if ~isempty(distance); % skip if the distance is not calculated
%                 figure;
%                 bar(distance/100);
%                 title(sprintf('%s day %s epoch %s',animalname,dayt,epocht));
%                 xlabel('Run');
%                 ylabel('Distance (m)');
%                 ylim([0 10]);
%                 xlim([0 size(distance,1)+1]);
%                 totaldist=sum(distance);
%                 xlimval=get(gca,'xlim');
%                 ylimval=get(gca,'ylim');
%                 distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.90*(ylimval(2)-ylimval(1))+ylimval(1)];
%                 text(distlabel(1),distlabel(2),sprintf('Tot. dist. %s m',num2str(totaldist/100,3)),'FontSize',15);
%             end
            cd(strcat(datadir,animaldir,'/'));
   
    plotdir = dir('Plot');
    if (isempty(plotdir))
        %an a plot folder needs to be created
        !mkdir Plot
    end
    
    % change to that directory and saves the figure with file name
    % animal_day_epoch 
    
    
    cd(strcat(strcat(datadir,animaldir,'/Plot/')));
    figurename = strcat('trajectory_',animalname,'_',dayt,'_',epocht);
    
    saveas(gcf, figurename, 'fig');
    saveas(gcf, figurename, 'pdf');
    
    %Closes the figure
    close;
     
            k=k+1;
        else display(sprintf('No %s data in epoch %s ',animalname,num2str(k)));
            k=k+1;
        end
        
    end
    
      
    
end
end


