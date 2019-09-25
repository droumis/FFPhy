% plot the xy position of the animal between adjacent trials

% need pos struct in workspace
% need unrewarded and rewarded times
    % vector of timestamps between events

   
%% Define variables
clear output;

if 1
d=3;
epochs=[2];

pos = loaddatastruct('/data12/mari/Egy/','egy','pos',d);

for e=epochs        % run over all epochs

    
%% REWARD    
% 1. (most inclusive, check first) To plot all reward inputs and outputs from sj_findwellsfromdio1_Egypt:
        % notice that we're keeping reward/error flag]
outputs=[rewardinfo{d}{e}(:,2);
inputs=rewardinfo{d}{e}(:,5);
times=sort([outputs ; inputs]);

%To plot only reward times that have 500 ms delay:
% outputs=rewardtimes{d}{e}(:,1);
% inputs=rewardtimes{d}{e}(:,2);
% times=sort([inputs ; outputs]);

%% ERROR
%To plot all error trial input times:
%  inputs=errortimesall{d}{e}(:,3);
%  times=sort(inputs);

% %To plot error trial input times after exclusions:
%inputs=errortimes{d}{e}(:,1);
%times=sort(inputs);

%To plot error trial pseudo output times:
% outputs=pseudotimes{d}{e}(:,1);
% times=sort(outputs);


%% DIO Directly
%To plot inputs or outputs from DIO:
% well8=DIO{d}{e}{9}.pulsetimes;   
% well9=DIO{d}{e}{10}.pulsetimes;  
% well8in=DIO{d}{e}{25}.pulsetimes;
% well9in=DIO{d}{e}{26}.pulsetimes;
% times=sort([well8;well9;well8in;well9in]);

%To plot inputs or outputs from DIO on linear track:
%  well1=DIO{d}{e}{25}.pulsetimes;   
% well2=DIO{d}{e}{26}.pulsetimes;  
% times=sort([well1;well2]);




%Plot times against xy position

counter=0;

for t=1:(length(times)-1)
        counter=counter+1;
        
    if (counter == 1)
        figure
        set(gcf,'Position',get(0,'Screensize'));
        subplot(4,6,counter)
    elseif (counter > 24)
        pause
        close
        figure
        set(gcf,'Position',get(0,'Screensize'));
        counter=1;
        subplot(4,6,counter)
    end
    

    subplot(4,6,counter)
    
    plot(pos{d}{e}.data(:,2),pos{d}{e}.data(:,3),'LineWidth',2,'Color',[0.8 0.8 0.8]);    % all paths in grey
    hold on
    startsample=times(t);        % random start index within the epoch
        starttime=startsample/10000;   % in sec, since open all files
        endsample=times(t+1);
    endtime=endsample/10000;
        pos_start=lookup(starttime,pos{d}{e}.data(:,1));                                                % index of start time in pos vector
        pos_end=lookup(endtime,pos{d}{e}.data(:,1));                                                    %index of end time in pos vector
    plot(pos{d}{e}.data(pos_start:pos_end,2),pos{d}{e}.data(pos_start:pos_end,3),'k','LineWidth',3);
    plot(pos{d}{e}.data(pos_start,2),pos{d}{e}.data(pos_start,3),'r*','LineWidth',3);                   %red dot=start time of interval

    plot(pos{d}{e}.data(pos_end,2),pos{d}{e}.data(pos_end,3),'b*','LineWidth',3);                       %blue dot=end time of interval
            string1 = sprintf('red/start %d',times(t));
            string2 = sprintf('blue/stop %d',times(t+1));
            title(char(string1,string2))

end

end
end

