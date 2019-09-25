% given a selected tetrode (e.g. CA1), calculates spectrograms time locked
% to event (startind)

%% THREE ingredients:
  %%  1. eeg struct     (all epochs for a day, use Mattias' loadeegstruct)
  %%  2. <event> struct  (after <event>dayprocess and <event>extract)
  %%  3. pos struct     (NOT rawpos)

%% Set these parameters manually.
    % Chronux params are specified below (fourth block)
animalname = 'Bashir';                              % to label graphs later
eventtype = 'gammal';                              % must be the same name as the events' data structure
days = 4;                                           % days to analyze
ref_tetrode = 10;                                   % reference tetrode
tetrodes = 8:14;                                    % tetrodes to analyze
epochs = 2:3;                                         % epochs to analyze
velocity_state_1 = 2;                               % state 1 (immobile)
velocity_state_2 = 8;                               % state 2 (run)
windowsize_sec = 0.8;                               %                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            % window length
Fs = 1500;                                          % sample rate of EEG
halfwindow_samp = (windowsize_sec*1500)/2;
tetno=tetrodes(end);
epno=epochs(end);
filename='';

switch eventtype
    case 'ripples'
        events=ripples;
    case 'gammal'
        events=gammal;
    case 'gammah'
        events=gammah;
    case 'gammaff'
        events=gammaff;
    case 'slowripples'
        events=slowripples;
    case 'hightheta'
        events=hightheta;
    otherwise
        disp('this event type is not recognized..'); 
end

%% First, collects all event timestamps from chosen ref. tetrode.

eventtimes=cell(days(end),epno);

for d=days
for e=epochs
    eventtimes{d,e}=[eventtimes{d,e}; events{d}{e}{ref_tetrode}.midtime];
end
end

%% Second, copy the raw eeg data into windows aligned to the event starttime.

eegwindows=cell(days(end),epno,tetno,3);

for d=days
for e=epochs
    epoch_starttime=eegstruct{d}{e}{14}.starttime;
    epoch_times=epoch_starttime:1/Fs:(epoch_starttime+(1/Fs)*length(eegstruct{d}{e}{14}.data));
    for t=tetrodes
        for r=(1+3):(length(eventtimes{d,e})-3)                                   % (!!) disregards first and last 5 events
            centerindex=lookup(eventtimes{d,e}(r),epoch_times);
            velocityindex=lookup(eventtimes{d,e}(r),pos{d}{e}.data(:,1));         %  added for state
            velocity=pos{d}{e}.data(velocityindex,5);
            s=0;
            if velocity<velocity_state_1           % state 1, immobile
                s=1;
            elseif velocity>velocity_state_2       % state 2, run
                s=2;
            else                                   % state 3, intermediate
                s=3;
            end
            eegwindows{d,e,t,s}=cat(1,eegwindows{d,e,t,s}, ...
                eegstruct{d}{e}{t}.data((centerindex-halfwindow_samp):(centerindex+halfwindow_samp))');
        end
    end
end
end

%% Third, tabulate the # of events in each of the 3 states.

statecounts=zeros(3,epno);

d=days
for e=epochs
        for r=(1+3):(length(eventtimes{d,e})-3)
            velocityindex=lookup(eventtimes{d,e}(r),pos{d}{e}.data(:,1));            %   added for state
            velocity=pos{d}{e}.data(velocityindex,5);
            if velocity<velocity_state_1           % state 1
                statecounts(1,e)=statecounts(1,e)+1;
            elseif velocity>velocity_state_2       % state 2
                statecounts(2,e)=statecounts(2,e)+1;
            else                                   % state 3
                statecounts(3,e)=statecounts(3,e)+1;
            end
        end
end

statecounts


%% Fourth, calculate all individual spectrograms.

allwindowedspectrograms=cell(days(end),epno,tetno,3);           % holds all events' spectograms

% chronux params
movingwin = [100 10]/1000;          
params.Fs = 1500;
params.err = [2 0.05];
params.fpass = [0 400];
params.tapers = [3 5];

for d=days
for e=epochs
    for t=tetrodes
        for s=1:3
            flag=0;
            for r=1:size(eegwindows{d,e,t,s},1)
                [S,times,freqs,Serr]=mtspecgramc(eegwindows{d,e,t,s}(r,:),movingwin,params);
                if flag==0                              % this if clause initializes the 3D matrix of all spectrograms for a given tetrode
                    allwindowedspectrograms{d,e,t,s}=nan(size(S,1),size(S,2),size(eegwindows{d,e,t,s},1));
                    flag=1;
                end
                allwindowedspectrograms{d,e,t,s}(:,:,r)=S;   % adds the S variable (2D matrix) to the third matrix dimension
            end
        end
    end
end
end




%% Fifth, calculates the continuous spectrogram, day mean, and std spectra for each tetrode.  (long calculation)

continuousspectrograms=cell(days(end),epno,tetno);
meandayspectra=cell(days(end),tetno);
stddayspectra=cell(days(end),tetno);

for d=days
    for t=tetrodes                                
        dummy=[];
        for e=epochs
            [S_full,junkt,junkf,junkserr] = mtspecgramc(eegstruct{d}{e}{t}.data,movingwin,params);
            dummy=[dummy;S_full];        
            continuousspectrograms{d,e,t}=S_full;
        end
    meandayspectra{t}=mean(dummy,1);
    stddayspectra{t}=std(dummy,1);
    end
end

% z-score the continuous spectrogram

for d=days
for e=epochs
    for t=tetrodes
        continuousspectrograms{d,e,t}=bsxfun(@minus,continuousspectrograms{d,e,t},meandayspectra{t});   
        continuousspectrograms{d,e,t}=bsxfun(@rdivide,continuousspectrograms{d,e,t},stddayspectra{t});
    end
end
end


disp('.finished with computing spectrograms.')
save(filename,'-7.3');


%% Sixth, z-scores the individual spectrograms ('zscorespectrograms').

zscorespectrograms=cell(days(end),epno,tetno,3);

for d=days
    for e=epochs
        for t=tetrodes
            for s=1:3
                if ~isempty(allwindowedspectrograms{d,e,t,s})       % in case there are no events
                    for r=1:size(allwindowedspectrograms{d,e,t,s},3)
                        zscorespectrograms{d,e,t,s}(:,:,r)=bsxfun(@minus,allwindowedspectrograms{d,e,t,s}(:,:,r),meandayspectra{t});
                        zscorespectrograms{d,e,t,s}(:,:,r)=bsxfun(@rdivide,allwindowedspectrograms{d,e,t,s}(:,:,r),stddayspectra{t});
                    end
                end
            end
        end
    end
end



%% Seventh, calculates mean spectrograms: (1) state, (2) epochs (3) overall.

meanspectrograms_state=cell(days(end),tetno,3);             %% discards epoch, pools all within state
meanspectrograms_epoch=cell(days(end),epno,tetno);          %% discards states, pools all within epoch
meanspectrograms_allevents=cell(days(end),tetno);           %% pools all events

epochs=[2];

for d=days
    for t=tetrodes
        
        for s=1:3
            dummy=[];
            for e=epochs
                if ~isempty(zscorespectrograms{d,e,t,s})
                    dummy=cat(3,dummy,zscorespectrograms{d,e,t,s});
                end
            end
            meanspectrograms_state{d,t,s}=mean(dummy,3);
            meanspectrograms_allevents{d,t}=cat(3,dummy,meanspectrograms_allevents{d,t});
        end
        
        meanspectrograms_allevents{d,t}=mean(meanspectrograms_allevents{d,t},3);
        
        for e=epochs
            dummy=[];
            for s=1:3
                if ~isempty(zscorespectrograms{d,e,t,s})
                    dummy=cat(3,dummy,zscorespectrograms{d,e,t,s});
                end
            end
            meanspectrograms_epoch{d,e,t}=mean(dummy,3);
        end
        
    end
end


    
%% Final variable list.

eventtimes;
eegwindows;
statecounts;
allwindowedspectrograms;
%continuousspectrograms;
meandayspectra;
stddayspectra;
zscorespectrograms;
meanspectrograms_state;
meanspectrograms_epoch;
meanspectrograms_allevents;


%% Ninth, plot mean spectrogram for each state.

for d=days
    for s=1:3
        figure
        
        for t=tetrodes
            h = subplot(1,8,t-7);
            imagesc(times,freqs,meanspectrograms_state{d,t,s}',[-0.2,4]);
            set(gca,'YDir','normal');
            string = sprintf('%d',t);
            title(string);
        end
        h = subplot(1,8,8);
        imagesc(times,freqs,meanspectrograms_state{d,t,s}',[-0.2,4]);
                    set(gca,'YDir','normal');
        colormap('jet')
        colorbar
        %%% TITLE code.
        eventcount = sum(statecounts(s,:));
        h = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
            1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.5, 1,sprintf('Mean %s spectrogram, %s, Day %d, reftet=%d, state %d (%d events)', ...
            eventtype,animalname,d,ref_tetrode,s,eventcount),'HorizontalAlignment'...
            ,'center','VerticalAlignment', 'top');
        
    end
end

%% (ALTERNATIVE) Plot mean spectrogram for each epoch.

for d=days
    for e=epochs
        figure
        
        for t=tetrodes
            h = subplot(1,11,t-12);
            imagesc(times,freqs,meanspectrograms_epoch{d,e,t}',[-0.2,5]);
            set(gca,'YDir','normal');
            string = sprintf('%d',t);
            title(string);
        end
        h = subplot(1,11,11);
        imagesc(times,freqs,meanspectrograms_epoch{d,e,t}',[-0.2,5]);
                    set(gca,'YDir','normal');
        colormap('jet')
        colorbar
        %%% TITLE code.
        eventcount = sum(statecounts(:,e));
        h = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
            1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.5, 1,sprintf('Mean %s spectrogram, %s, Day %d, reftet=%d, epoch %d (%d events)', ...
            eventtype,animalname,d,ref_tetrode,e,eventcount),'HorizontalAlignment'...
            ,'center','VerticalAlignment', 'top');
        
    end
end


%% Tenth, plot some random individual event spectrograms.
    
    % set these manually to your taste
    d=4
    e=2;                                 % epoch manually
    events_plot=12;                         % no. events to plot
    tetrodes_plot=8:14;                  % tetrodes to plot
    states_plot=1;
    %
    
    for s=states_plot
        for i=1:events_plot
            N=size(zscorespectrograms{d,e,t,s},3);        % # events to choose from
            eventno=ceil(N*rand);
            figure
            title(string);
            for t=tetrodes_plot
                subplot(1,8,t-7)
                imagesc(times,freqs,zscorespectrograms{d,e,t,s}(:,:,eventno)',[-0.2,4]);
                colormap('hot')
                set(gca,'YDir','normal');
                string = sprintf('%d',t);
                title(string);
            end
            subplot(1,8,8)
            imagesc(times,freqs,zscorespectrograms{d,e,t,s}(:,:,eventno)',[-0.2,4]);
            colorbar
            %%% TITLE code.
            h = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
                1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
            text(0.5, 1,sprintf('Single %s event spectrogram, %s, Day %d, reftet=%d, state %d (event #%d)', ...
                eventtype,animalname,d,ref_tetrode,s,eventno),'HorizontalAlignment'...
                ,'center','VerticalAlignment', 'top');
        end
    end


%% Eleventh, plot the continuous spectrogram (random periods).

d=4;
e=2;
plot_duration=60;    % in seconds
    plotlength_spectra=floor(plot_duration/movingwin(2));      % # of spectra (vertical lines)
    plotlength_possamp=floor(plot_duration*29.97003);
noplots=1;

for i=1:5
    figure
    
    eegdata_startindex=lookup(pos{d}{e}.data(1,1),(eegstruct{d}{e}{14}.starttime):1/Fs: ...       % eeg index at start of pos data
                      (eegstruct{d}{e}{14}.starttime+length(eegstruct{d}{e}{14}.data)/Fs));
    nospectra=size(continuousspectrograms{d,e,14},1);    % # of spectra, which should = 1+(length(eeg.data)/Fs-movingwin(1))/movingwin(2)
    no_deadspectra=round((eegdata_startindex/Fs-movingwin(1))/movingwin(2)+1);              % spectra w/o corresponding pos data
    
    startindex_spectra=floor(rand*(nospectra - ...
                                   plotlength_spectra  - ...
                                   no_deadspectra));            % start index (in spectra indices) of first spectrum to plot, but given
                                                                            % that the spectra start at the first available pos window
                                                                   
    true_time=startindex_spectra*movingwin(2)+movingwin(1)/2;            % time elapsed since start of pos data
    
    pos_index=round(true_time*29.97003);
    spectra_index=no_deadspectra+startindex_spectra;
    
        % velocity plot
    
    subplot(8,1,1)
    plot((1:(plotlength_possamp+1))/29.97003, ...
            pos{d}{e}.data(pos_index:(pos_index+plotlength_possamp),5));
    axis tight;
    ylim([0 50]);
    
    
        % spectrogram plot
    
    for t=8:14
        subplot(8,1,t-6)
        imagesc((1:plotlength_spectra)*movingwin(2),junkf(1:100), ...
            continuousspectrograms{d,e,t}(spectra_index:(spectra_index+plotlength_spectra),:)',[-0.5,2]);
        colormap('bone')
        set(gca,'YDir','normal');
        string = sprintf('%d',t);
    end
end