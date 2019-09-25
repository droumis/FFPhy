%% plots selected tetrodes' EEG side-by-side for, random start times

function out = plotEVENTS_Bashir(eegstruct,eventstruct,posstruct,eventname, ...
    day,epoch,tetrodes,reftetrode,windowsize,no_events,varargin)

%% eegstruct can be the raw eeg, or a filtered eeg (consolidated struct)
%% eventstruct contains times of ripples, gammal, hightheta, etc.
%% tetrodes is a vector of tetrodes
%% windowsize is in s

animalname='Bashir'
num_tetrodes=length(tetrodes);
Fs=eegstruct{day}{epoch}{14}.samprate;
halfwin=floor(windowsize*Fs/2);            %% window in samples
velocity_state_1=2;
velocity_state_2=8;
CA1_tetrodes=14;
CA3_tetrodes=12;
CA2_tetrodes=9;
nearCA2_tetrodes=13;
DG_tetrodes=10;

state=0;

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'state'
            state = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

if state==0            % plot regardless of state
    total_events=length(eventstruct{day}{epoch}{reftetrode}.midind);
    events=floor((rand(1,no_events))*total_events);                             % chooses random events to plot
else
    velocities=posstruct{day}{epoch}.data(:,5);                                 % all velocities
    position_indices=eventstruct{day}{epoch}{reftetrode}.posind;                % position indices of all events
        valid_event_indices=find(~isnan(position_indices));                            % for non-position data available eventss
        event_velocities=velocities(valid_event_indices);  % instantaneous velocities of events
    if state==1
        valid_events=event_velocities<velocity_state_1;                         % logical
    elseif state==2
        valid_events=event_velocities>velocity_state_2;
    else
        valid_events=event_velocites>velocity_state_1 || event_velocities<velocity_state2;
    end
    total_events=sum(valid_events);
    events=find(valid_events);
    events=events(floor(rand(1,no_events)*total_events)+1);                     % eventstruct indices of events to be plotted
end



for k=1:length(events)
    figure
    hold on
    index=eventstruct{day}{epoch}{reftetrode}.midind(events(k));
    pos_index=eventstruct{day}{epoch}{reftetrode}.posind(events(k));
    halfwindow_possamp=round(windowsize*29.97003/2);
    velocity=posstruct{day}{epoch}.data(pos_index,5);      % inst. v
    
    subplot(4,1,1)
    plot( ... % (1:(halfwindow_possamp*2+1))/29.97003, ...
           posstruct{day}{epoch}.data((pos_index-halfwindow_possamp):(pos_index+halfwindow_possamp),5));
    axis tight
    ylim([0 60])
    
    subplot(4,1,2:4) 
    hold on
    for t=tetrodes
        if sum(t==CA3_tetrodes)         %% CA3
            plot((1:(1+2*halfwin))/Fs-windowsize/2, ...
                eegstruct{day}{epoch}{t}.data((index-halfwin):(index+halfwin))-1000, ...
                'r','LineWidth',2);
        elseif sum(t==DG_tetrodes)     %% DG
            plot((1:(1+2*halfwin))/Fs-windowsize/2, ...
                eegstruct{day}{epoch}{t}.data((index-halfwin):(index+halfwin))-2000, ...
                'm','LineWidth',2);
%         elseif sum(t==nearCA2_tetrodes)     %% near CA2
%             plot((1:(1+2*halfwin))/Fs-windowsize/2, ...
%                 eegstruct{day}{epoch}{t}.data((index-halfwin):(index+halfwin))+1000, ...
%                 'Color',[0 0.5 0],'LineWidth',2);               
        elseif sum(t==CA2_tetrodes)     %% CA2
            plot((1:(1+2*halfwin))/Fs-windowsize/2, ...
                eegstruct{day}{epoch}{t}.data((index-halfwin):(index+halfwin))+1000, ...
                'g','LineWidth',2);            
        elseif sum(t==CA1_tetrodes)                    %% CA1
            plot((1:(1+2*halfwin))/Fs-windowsize/2, ...
                eegstruct{day}{epoch}{t}.data((index-halfwin):(index+halfwin)), ...
                'k','LineWidth',2);
        end
    end
    axis tight
    ylim([-3200 2200]);
    %%% TITLE code.
    h = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
        1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    varname=@(eventname) inputname(1);
    text(0.5, 1,sprintf('single %s event, Day %d, reftet=%d, velocity %f', ...
        eventname,day,reftetrode,velocity),'HorizontalAlignment'...
        ,'center','VerticalAlignment', 'top');


end






