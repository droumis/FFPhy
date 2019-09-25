%plots position at reward times and audio cue times

%%Three ingredients:
%%  1. <event> struct i.e. DIO struct (after <event>dayprocess and <event>extract)
%%  2. pos struct     (NOT rawpos)
%%  3. eegstruct      (for epoch_starttime) - use loadeegstruct
%% Set these parameters manually

animalname = 'Egypt';                              % to label graphs later
eventtype = 'reward';                              % must be the same name as the events' data structure
reward_bits = 9:11;
days = 8;                                           % days to analyze                                    % tetrodes to analyze
epochs = 2;                                         % epochs to analyze: should be run epochs
Fs = 1500;                                          %eeg sampling rate
Ps = 30;                                            %position (mpeg) sampling rate: 30 frames per second

switch eventtype
    case 'reward'
        events=DIO;
end


%%
%%First, collect reward output times as eventtimes
eventtimes=cell(days(end),epochs(end));

for d=days
    for e=epochs
        for w=reward_bits
            eventtimes{d,e}=[eventtimes{d,e}; events{d}{e}{w}.pulsetimes/10000]; % converts pulsetimes to seconds
        end
    end
end


%%
%%Second, find when position times = eventtimes{d,e}
%Remember pos.data is an array where columns are: 'time x y dir velocity'

postimes=cell(days(end),epochs(end));

for d=days
    for e=epochs
        postimes{d,e}=pos{d}{e}.data(:,1); % assigns position times to the variable postimes
    end
end

postimestamps=postimes{d,e};
eventtimestamps(1:53,1)=eventtimes{d,e}(1:53);
eventtimestamps(54:106,1)=eventtimes{d,e}(54:106);
evtimestamps=sort(eventtimestamps(:,1));

% for i = 1:length(evtimestamps)
%     for ets = evtimestamps(i)
%             I=find(postimestamps==ets);      
%     end
% end
% 
% for ets = evtimestamps(:,1)
%     if ets==postimestamps(:,1)
%         plot(pos{d}{e}.data(:,2), pos{d}{e}.data(:,3),'kd','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize', 5)
%     end
% end
% sort(evtimestamps);
% wells = ismember(postimestamps,evtimestamps);
% 
% timestampfinder=find(wells);        
%         
% 
% for d=days
%     for e=epochs
%         epoch_starttime=eegstruct{d}{e}{14}.starttime;          % time, in seconds, since Open All Files
%         posindex=(eventtimes{d,e}-epoch_starttime)/Fs;
%     end
% end




figure
hold on

plot(pos{d}{e}.data(:,2), pos{d}{e}.data(:,3),'b')

% if postimes{d}{e}=eventtimes{d}{e}
%     plot(pos{d}{e}.data(:,2), pos{d}{e}.data(:,3),'kd','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize', 5)
% end
