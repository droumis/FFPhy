
% This script plots the xy position of the animal between adjacent trials

    % plots -all- trajectories between adjacent timestamps of input stamps
        % vector
    

day=3;
epochs=2;
pos = loaddatastruct('/data12/mari/Egy/','egy','pos',day);


for e=epochs        % run over all epochs

% 1. (most inclusive, check first) 
%     Plot all trajectories from sj_findwellsfromdio1_Egypt:
                
 outputs=[rewardinfo{day}{e}(:,5) rewardinfo{day}{e}(:,3) rewardinfo{day}{e}(:,4)];
 inputs=[rewardinfo{day}{e}(:,2) rewardinfo{day}{e}(:,3) rewardinfo{day}{e}(:,4)];
    stamps=sortrows([outputs ; inputs],1);

%To plot only reward times that have 500 ms delay:
% outputs=rewardtimes{day}{e}(:,1);
% inputs=rewardtimes{day}{e}(:,2);
% times=sort([inputs ; outputs]);

%To plot all error trial input times:
%  inputs=errortimesall{day}{e}(:,3);
%  times=sort(inputs);

% %To plot error trial input times after exclusions:
%inputs=errortimes{day}{e}(:,1);
%times=sort(inputs);

%To plot error trial pseudo output times:
% outputs=pseudotimes{day}{e}(:,1);
% times=sort(outputs);

%To plot inputs or outputs from DIO:
% well8=DIO{day}{e}{9}.pulsetimes;   
% well9=DIO{day}{e}{10}.pulsetimes;  
% well8in=DIO{day}{e}{25}.pulsetimes;
% well9in=DIO{day}{e}{26}.pulsetimes;
% times=sort([well8;well9;well8in;well9in]);

%To plot inputs or outputs from DIO on linear track:
%  well1=DIO{day}{e}{25}.pulsetimes;   
% well2=DIO{day}{e}{26}.pulsetimes;  
% times=sort([well1;well2]);

%Plot times against xy position

counter=0;

for t=1:(length(times)-1) 
        counter=counter+1;
        
    if (counter == 1)
        figure
        set(gcf,'Position',get(0,'Screensize'));
        subplot(4,6,counter)
    elseif (counter > 9)
        pause
        close
        figure
        set(gcf,'Position',get(0,'Screensize'));
        counter=1;
        subplot(3,3,counter)
    end
    
    subplot(3,3,counter)
    
    % plot all traversed paths as light grey
    plot(pos{day}{e}.data(:,2),pos{day}{e}.data(:,3),'LineWidth',2,'Color',[0.8 0.8 0.8]);    
    hold on

    % check if rewarded (0 if not, 1 if yes)
    rewarded = stamps(t+1,2);
    if rewarded
        % change background color if rewarded
        set(gca,'Color',[1 1 0.5]);
    end
    
    % obtain corresponding position indices
    starttime=stamps(t,1)/10000;
        pos_start=lookup(starttime,pos{day}{e}.data(:,1));  
    endtime=stamps(t+1,1)/10000;
        pos_end=lookup(endtime,pos{day}{e}.data(:,1));
    diff = endtime-starttime;
        
    % plot trajectory
    plot(pos{day}{e}.data(pos_start:pos_end,2),pos{day}{e}.data(pos_start:pos_end,3),'k','LineWidth',2);
    plot(pos{day}{e}.data(pos_start,2),pos{day}{e}.data(pos_start,3),'r*','LineWidth',2);                   %red dot=start time of interval

    % draw colored dots at start and end
    plot(pos{day}{e}.data(pos_end,2),pos{day}{e}.data(pos_end,3),'b*','LineWidth',3);                       %blue dot=end time of interval
            string1 = sprintf('start (red) %d',times(t));
            string2 = sprintf('stop (blue) %d',times(t+1));
            title(char(string1,string2,num2str(diff)))
    
    % label inbound or outbound
    if stamps(t+1,3)==10
        text(70,65,'in','FontSize',18,'FontWeight','bold')
    elseif stamps(t+1,3)==11
        text(70,65,'out','FontSize',18,'FontWeight','bold')
    end
            
    ylim([0 100])
    xlim([25 125])
    

end

end

