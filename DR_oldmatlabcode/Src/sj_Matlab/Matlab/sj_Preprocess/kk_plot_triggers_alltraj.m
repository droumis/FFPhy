
% This script plots the xy position of the animal between adjacent trials

% plots -all- trajectories between adjacent timestamps of input stamps
% vector


day=3;
epochs=2;
pos = loaddatastruct('/data12/mari/Egy/','egy','pos',day);


for e=epochs        % run over all epochs
    
    %  Plot all trajectories from sj_findwellsfromdio1_Egypt:
    %  Note that this means that you'll be plotting parsed events (all the
    %  timestamps reported in rewardinfo) in a ZIGZAG, from 5th column to 2nd
    %  column, then 2nd column to 5th column on next row, and so on
    
    data=rewardinfo{day}{e};
    
    %iterate through rows, though will plot two trajectories per row
    for r=1:size(data,1)
        figure
        % 1 corresponds to input-to-output traj         (horizontal zig)
        % 2 corresponds to output-to-next-input traj    (diagonal zag)
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
            
            % check if rewarded, change bckgd color if so
            if ii==1
                if data(r,3)==1
                    set(gca,'Color',[1 1 0.3]);   % rewarded is yellow
                elseif data(r,3)==0
                    set(gca,'Color',[1 .85 .85])  % error is red
                end
            end
            
            % obtain corresponding position indices
            if ii==1                    % zig
                starttime=data(r,5);
                pos_start=lookup(starttime/10000,pos{day}{e}.data(:,1));
                endtime=data(r,2);
                pos_end=lookup(endtime/10000,pos{day}{e}.data(:,1));
                diff = endtime-starttime;
            elseif ii==2                % zag
                starttime=data(r,2);
                pos_start=lookup(starttime/10000,pos{day}{e}.data(:,1));
                endtime=data(r+1,5);
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
            title(char(string1,string2,num2str(diff/10000)),'FontSize',17,'FontWeight','bold')
            % if no delay, then print in red
            if diff==0
                title(char(string1,string2,'no delay -- should be error'),'FontSize',17,'FontWeight','bold','Color','r')
            end
                
            
            % label inbound or outbound
            if ii==2
                if data(r+1,4)==10
                    text(70,65,'in','FontSize',18,'FontWeight','bold')
                elseif data(r+1,4)==11
                    text(70,65,'out','FontSize',18,'FontWeight','bold')
                end
            end
            
            ylim([0 100])
            xlim([25 125])
            
            % stop after each trajectory
            if ii==2
                pause
                close
            end
            
        end
        
    end
end
    
