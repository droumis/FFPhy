%% plots selected tetrodes' EEG side-by-side for DIO pulse times

function out = plotEEG_DIO(eegstruct,eventstruct,refbit,day,epoch,tetrodes,windowsize,varargin)

%% eegstruct can be the raw eeg, or a filtered eeg (consolidated struct)
%% eventstruct contains times of ripples, gammal, hightheta, etc.
%% tetrodes is a vector of tetrodes
%% windowsize is in s

Fs=eegstruct{day}{1}{14}.samprate;
halfwin=floor(windowsize*Fs/2);            %% window in samples
random_flag=0;
no_events=10;                            % if want some subset of all events

    pulsetimes=eventstruct{day}{epoch}{refbit}.pulsetimes(:,1);
    
if random_flag==1                                       % plot regardless of state
   events=floor((rand(1,no_events))*pulsetimes);                             % chooses random events to plot
else
   events=pulsetimes;
end


    hold on

for k=1:length(events)
        figure
    epoch_starttime=eegstruct{day}{epoch}{14}.starttime;
    epoch_times=epoch_starttime:1/Fs:(epoch_starttime+(1/Fs)*length(eegstruct{day}{epoch}{14}.data));
    index=lookup(pulsetimes(k)/10000,epoch_times);
    halfwindow_possamp=round(windowsize*29.97003/2);
    for t=tetrodes
            plot(1000*((1:(1+2*halfwin))/Fs-windowsize/2), ...
                eegstruct{day}{epoch}{t}.data((index-halfwin):(index+halfwin))-1000, ...
                'k','LineWidth',2);
    end
end






