%% plots selected tetrodes' EEG side-by-side for, random start times

function out = plot_events_probe(animalprefix,eegstruct,eventstruct,sitemap,session,epoch,channels,refchan,windowsize,no_events)

%% eegstruct can be the raw eeg, or a filtered eeg (consolidated struct)
%% eventstruct contains times of ripples, gammal, hightheta, etc.
%% windowsize is in s

Fs=1500;
windowsizesamp=floor(windowsize*Fs);            %% window in samples
    totalnoevents=length(eventstruct{session}{epoch}{refchan}.midind);
    events=floor((rand(1,no_events))*totalnoevents);                             % chooses random events to plot
    yshift=350;
    xshift=windowsize*1.2;
    nochanshank=size(sitemap,1);
    noshanks=size(sitemap,2);
    
for k=1:length(events)
    figure
    hold on
    index=eventstruct{session}{epoch}{refchan}.midind(events(k));
    
    for s=1:noshanks
    for z=1:nochanshank
        
        c=sitemap(z,s);
        color='k';
        if c==refchan
            color='r';
        end
        if ~isnan(c)
            plot((0:windowsizesamp)/Fs + xshift*(s-1), ...
                eegstruct{session}{epoch}{c}.data((index-windowsizesamp/2):(index+windowsizesamp/2))-z*yshift,color,'LineWidth',2)
        end
    end
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






