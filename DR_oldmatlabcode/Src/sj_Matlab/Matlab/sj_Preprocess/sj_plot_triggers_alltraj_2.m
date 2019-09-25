

% Change from Kenny's - Kenny was plotting inbound triggers - outbound
% triggers (Column 2 - Column 5). Focus on 1 column: 2 or 5, and plot
% trajectories


% This script plots the xy position of the animal between adjacent trials

% plots -all- trajectories between adjacent timestamps of input stamps
% vector


day=2;
epochs=2;
%pos = loaddatastruct('/data12/mari/Egy/','egy','pos',day);
pos = loaddatastruct('/data25/sjadhav/HPExpt/HPb_direct/','HPb','pos',day);

plotcolumn = 5; % column 2 is input triggers, in column 5, rewarded is replaced with output triggers 

for e=epochs        % run over all epochs
    
    %  Plot all trajectories from sj_findwellsfromdio1_Egypt:
    %  Note that this means that you'll be plotting parsed events (all the
    %  timestamps reported in rewardinfo) in a ZIGZAG, from 5th column to 2nd
    %  column, then 2nd column to 5th column on next row, and so on
    
    data=rewardinfo{day}{e};
    
    %iterate through rows, though will plot two trajectories per row
    for r=1:size(data,1)
        figure
        % 1 corresponds to previous-to-current traj         
        % 2 corresponds to current-to-next- traj    
        for ii=1:2

            subplot(1,2,ii)
            set(gcf,'Position',get(0,'Screensize'));
            
            % if last event in epoch, then display 'done'
            try
                data(r+1,:);
            catch
                if ii==2
                    figure
                    disp('done with all entries')
                    continue
                end
            end
            
           % plot all traversed paths as light grey
            plot(pos{day}{e}.data(:,2),pos{day}{e}.data(:,3),'LineWidth',2,'Color',[0.8 0.8 0.8]);
            hold on
            
            % obtain corresponding position indices
            
            
            if ii==1                    % previous to current trajectory
                if r==1
                    starttime = data(r,plotcolumn);
                else
                    starttime = data(r-1,plotcolumn);
                end
                pos_start=lookup(starttime/10000,pos{day}{e}.data(:,1));
                endtime=data(r,plotcolumn);
                pos_end=lookup(endtime/10000,pos{day}{e}.data(:,1));
                diff = endtime-starttime;
            elseif ii==2                % current to next trajectory
                starttime=data(r,plotcolumn);
                pos_start=lookup(starttime/10000,pos{day}{e}.data(:,1));
                if r == size(data)
                    endtime = data(r,plotcolumn);
                else
                    endtime = data(r+1,plotcolumn);
                end
                pos_end=lookup(endtime/10000,pos{day}{e}.data(:,1));
                diff = endtime-starttime;
            end
    
       
            % plot trajectory
            plot(pos{day}{e}.data(pos_start:pos_end,2),pos{day}{e}.data(pos_start:pos_end,3),'k','LineWidth',2);
            plot(pos{day}{e}.data(pos_start,2),pos{day}{e}.data(pos_start,3),'r*','LineWidth',2);                   %red dot=start time of interval
            
            % draw colored dots at start and end
            plot(pos{day}{e}.data(pos_end,2),pos{day}{e}.data(pos_end,3),'b*','LineWidth',3);                       %blue dot=end time of interval
            string1 = sprintf('start (red) %d',starttime);
            string2 = sprintf('stop (blue) %d',endtime);   
            
            
            % check if current (previous to current) is rewarded, change bckgd color if so
            if ii==1 % label inbound or outbound
                if data(r,4)==10
                    text(70,65,'in','FontSize',18,'FontWeight','bold')
                elseif data(r,4)==11
                    text(70,65,'out','FontSize',18,'FontWeight','bold')
                end
                
                if data(r,3)==1
                    set(gca,'Color',[1 1 0.3]);   % rewarded is yellow
                    %title(char(string1,string2,'Correct trial - got reward'),'FontSize',17,'FontWeight','bold','Color','y')
                elseif data(r,3)==0
                    set(gca,'Color',[1 .85 .85])  % error is red
                    %title(char(string1,string2,'Error trial - no reward'),'FontSize',17,'FontWeight','bold','Color','r')
                end
                
                title(char('Curr traj',string1,string2,num2str(diff/10000)),'FontSize',17,'FontWeight','bold')
            end
            
            % check if current to next is rewarded, change bckgd color if so
            if ii==2  
                if data(r+1,4)==10
                    text(70,65,'in','FontSize',18,'FontWeight','bold')
                elseif data(r+1,4)==11
                    text(70,65,'out','FontSize',18,'FontWeight','bold')
                end
                
                if r~=size(data,1) % if last well, no next trajectory
                    if data(r+1,3)==1
                        set(gca,'Color',[1 1 0.3]);   % rewarded is yellow
                        %title(char(string1,string2,'Correct trial - got reward'),'FontSize',17,'FontWeight','bold','Color','y')
                    elseif data(r+1,3)==0
                        set(gca,'Color',[1 .85 .85])  % error is red
                        %title(char(string1,string2,'Error trial - no reward'),'FontSize',17,'FontWeight','bold','Color','r')
                    end
                end
                
                title(char('Next traj',string1,string2,num2str(diff/10000)),'FontSize',17,'FontWeight','bold')
            end

            
%             % if no delay, then print in red
%             if diff==0
%                 title(char(string1,string2,'no delay -- should be error'),'FontSize',17,'FontWeight','bold','Color','r')
%             end
            
            
            
            
            ylim([0 120])
            xlim([40 160])
            
            % stop after each trajectory
            if ii==2
                pause
                close
            end
            
        end
        
    end
end
    
