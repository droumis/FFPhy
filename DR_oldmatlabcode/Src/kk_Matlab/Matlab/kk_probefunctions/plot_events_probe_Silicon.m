%% plots selected tetrodes' EEG side-by-side for, random start times

function out = plot_events_probe(eegstruct,eventstruct,session,epoch,channels,refchan,windowsize,no_events)

%% eegstruct can be the raw eeg, or a filtered eeg (consolidated struct)
%% eventstruct contains times of ripples, gammal, hightheta, etc.
%% windowsize is in s

animalname='Mitt';
Fs=1000;
windowsizesamp=floor(windowsize*Fs);            %% window in samples
    totalnoevents=length(eventstruct{session}{epoch}{refchan}.midind);
    events=floor((rand(1,no_events))*totalnoevents);                             % chooses random events to plot
    yshift=1000;

for k=1:length(events)
    figure
    hold on
    index=eventstruct{session}{epoch}{refchan}.midind(events(k));
    for c=1:15
        skipflag=0;
        if c>2
            skipflag=1;
        end
        if c==refchan
            plot((0:windowsizesamp)/Fs, ...
                eegstruct{session}{epoch}{c}.data((index-windowsizesamp/2):(index+windowsizesamp/2))-(c+skipflag)*yshift,'r','LineWidth',2)
        else
            plot((0:windowsizesamp)/Fs, ...
                eegstruct{session}{epoch}{c}.data((index-windowsizesamp/2):(index+windowsizesamp/2))-(c+skipflag)*yshift,'k','LineWidth',2)
        end
       
    end
    
    for c=16:31
            plot((0:windowsizesamp)/Fs + windowsize*1.2, ...
                 eegstruct{session}{epoch}{c}.data((index-windowsizesamp/2):(index+windowsizesamp/2))-(c-15)*yshift,'k','LineWidth',2)        
    end
    
    axis tight
    
%     ylim([-3200 2200]);
%     %%% TITLE code.
%     h = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
%         1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%     varname=@(eventname) inputname(1);
%     text(0.5, 1,sprintf('single %s event, Day %d, reftet=%d, velocity %f', ...
%         eventname,day,reftetrode,velocity),'HorizontalAlignment'...
%         ,'center','VerticalAlignment', 'top');


end






