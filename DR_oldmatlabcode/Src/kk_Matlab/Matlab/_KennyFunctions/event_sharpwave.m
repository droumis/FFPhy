% . Computes signed amplitudes of LFP during <event>: "sharpwaves_ref"

% . Preliminarily creates a sharp-wave version of the ripple struct,
%    i.e. simply analyze SPW from the same-tetrode SWR (no referencing):
%    "sharpwave"

%% Three ingredients:
  %%  1. eeg struct     (all epochs for a day, use Mattias' loadeegstruct)
  %%  2. <event> struct (after <event>dayprocess and <event>extract)
  %%  3. pos struct     (NOT rawpos)

%% Set these parameters manually.

Fs = 1500;                                          % sample rate of EEG
baseline = 150;                                      % length of baseline, in ms 
    baseline = round(80*Fs/1000);
animalname = 'Bashir';                              % to label graphs later
eventtype = 'ripples';                              % must be the same name as the events' data structure
days = 4;                                           % days to analyze
ref_tetrode = 14;                                    % reference tetrode
tetrodes = 1:14;                                    % tetrodes to analyze
epochs = 2:3;                                       % epochs to analyze
    tetno=tetrodes(end);
    epno=epochs(end);
velocity_state_1 = 2;                               % state 1 (immobile)
velocity_state_2 = 8;                               % state 2 (run)
windowsize_sec = 2;                               %
    halfwindow_samp = (windowsize_sec*1500)/2;          % window length

switch eventtype
    case 'ripples'
        events=ripples;
    case 'gammal'
        events=gammal;
    case 'gammah'
        events=gammaf;
    case 'gammaff'
        events=gammaff;
    case 'hightheta'
        events=hightheta;
    otherwise
        disp('this event type is not recognized..'); 
end

%% First, compute sharpwaves. (same-tetrode-ripple referencing)

sharpwaves=cell(days(end),epno,tetno); 

for d=days
    for e=epochs
        for t=tetrodes
            for r=1:length(events{d}{e}{t}.startind)
                start_ind=events{d}{e}{t}.startind(r);
                end_ind=events{d}{e}{t}.endind(r);
                baseline_mean = mean(eegstruct{d}{e}{t}.data((start_ind-baseline):(start_ind-1)));
                event_mean = mean(eegstruct{d}{e}{t}.data(start_ind:end_ind));
                sharpwaves{d}{e}{t}.amplitudes(:,r)=event_mean-baseline_mean;
            end
        end
    end
end

% plot

figure
for d=days
    for e=2
          [dummy,bincenters]=hist(sharpwaves{d}{e}{14}.amplitudes,100);   % pick a nice CA1 tetrode (here, tet 14) to get histogram bin centers
        for t=8:14
            %subplot(14,1,t)
                        subplot(7,1,t-7)
            if t==12
                clr='r';
            elseif t==9
                clr='g';
            elseif t==8 || t==13
                clr=[0 0.5 0];
            elseif t==10
                clr='m';
            else
                clr='k';
            end
               hist(sharpwaves{d}{e}{t}.amplitudes,bincenters)
               h=findobj(gca,'Type','patch');
               set(h,'FaceColor',clr,'EdgeColor','w','facealpha',1)
               ylim([0 20])
        end
    end
end


   
%% Second, makes sharpwaves_ref cell array of structs:
      %  i. copies raw eeg windows referenced on <event> startind
      %          also, transcribes event duration
      %  ii. assigns state
      %  iii. calculates signamp
      %  iv. transcribes state1 and state2 velocity thresholds & reftet

      sharpwaves_ref=cell(days(end),epno,tetno);
      
      for d=days
          for e=epochs
              event_start=events{d}{e}{ref_tetrode}.startind;
              event_durations=events{d}{e}{ref_tetrode}.endind-events{d}{e}{ref_tetrode}.startind;
              event_midtimes=events{d}{e}{ref_tetrode}.midtime;
              
              for t=tetrodes
                  
                  sharpwaves_ref{d}{e}{t}.data=[];
                  sharpwaves_ref{d}{e}{t}.eventduration_samp=[];
                  sharpwaves_ref{d}{e}{t}.state=[];
                  sharpwaves_ref{d}{e}{t}.signamp=[];
                  
                  for r=2:length(events{d}{e}{ref_tetrode}.startind)            % (!!) skips the very first event
                      
                      sharpwaves_ref{d}{e}{t}.data(r,:) = ...                                   %  i. copy raw eeg windows
                                eegstruct{d}{e}{t}.data((event_start(r)-halfwindow_samp):(event_start(r)+halfwindow_samp));
                      sharpwaves_ref{d}{e}{t}.eventduration_samp = [sharpwaves_ref{d}{e}{t}.eventduration_samp ; event_durations(r)];           % transcribes duration of (ref_tetrode's) event
                      
                      velocityindex=lookup(event_midtimes(r),pos{d}{e}.data(:,1));              %  ii. classify the event's velocity / state
                      velocity=pos{d}{e}.data(velocityindex,5);
                      if velocity<velocity_state_1           % state 1, immobile
                          s=1;
                      elseif velocity>velocity_state_2       % state 2, run
                          s=2;
                      else                                   % state 3, intermediate
                          s=3;
                      end
                      sharpwaves_ref{d}{e}{t}.state=[sharpwaves_ref{d}{e}{t}.state; s];
                      
                      
                      eventlfp = mean(sharpwaves_ref{d}{e}{t}.data(r,(halfwindow_samp+1):(event_durations(r)+halfwindow_samp+1))); % event duration
                      baselinelfp = mean(sharpwaves_ref{d}{e}{t}.data(r,(halfwindow_samp-baseline):halfwindow_samp));  % baseline
                      sharpwaves_ref{d}{e}{t}.signamp=[sharpwaves_ref{d}{e}{t}.signamp ; eventlfp-baselinelfp];             % iii. calculate signamp of spw
                  end
                  
                  
                  sharpwaves_ref{d}{e}{t}.reftet=ref_tetrode;                                    % iv. transcribing stuff
                  sharpwaves_ref{d}{e}{t}.state1velocity=velocity_state_1;
                  sharpwaves_ref{d}{e}{t}.state2velocity=velocity_state_2;
              end
          end
      end

%% Third, tabulate the # of events in each of the 3 states.

statecounts=zeros(3,epno);

for e=epochs
        eventtimes=events{d}{e}{ref_tetrode}.midind;
        for r=1:length(eventtimes)            %% inst. velocity at the midpoint of the event
            velocityindex=lookup(eventtimes(r),pos{days}{e}.data(:,1));            %   added for state
            velocity=pos{days}{e}.data(velocityindex,5);
            if velocity<velocity_state_1           % state 1
                statecounts(1,e)=statecounts(1,e)+1;
            elseif velocity>velocity_state_2       % state 2
                statecounts(2,e)=statecounts(2,e)+1;
            else                                   % state 3
                statecounts(3,e)=statecounts(3,e)+1;
            end
        end
end



%% Fifth, calculates mean spectrogram of all events over all epochs, within each state.

meaneventspectrograms_epochlump = cell(1,tetno,3);   % mean spectrogram for all events from each state

for t=tetrodes
    for s=1:3
        dummy=[];
             for e=epochs
                 dummy=cat(3,dummy,allwindowedspectrograms{1,e,t,s});   % concatenate events over epoch
             end
        meaneventspectrograms_epochlump{1,t,s}=mean(dummy,3);
    end
end

%% Sixth, calculates day mean and std spectra for each tetrode.  (long calculation)

continuousspectrograms=cell(1,epno,tetno);
meandayspectra=cell(1,tetno);
stddayspectra=cell(1,tetno);

for d=days
    for t=tetrodes                                
        dummy=[];
        for e=epochs
            [S_full,junkt,junkf,junkserr] = mtspecgramc(eegstruct{days}{e}{t}.data,movingwin,params);
            dummy=[dummy;S_full];        
            continuousspectrograms{1,e,t}=S_full;
        end
    meandayspectra{t}=mean(dummy,1);
    stddayspectra{t}=std(dummy,1);
    end
end

save('Bashir_day04_state_gammaltriggeredspectrograms_ref14','-v7.3');

%% Seventh, z-scores both the event-triggered ('zscorespectrograms') and continuous spectrograms ('continuousspectrograms').

zscorespectrograms=cell(1,epno,tetno,3);

for e=epochs
    for t=tetrodes
        for s=1:3
            for r=1:size(allwindowedspectrograms{1,e,t,s},3)
                zscorespectrograms{1,e,t,s}(:,:,r)=bsxfun(@minus,allwindowedspectrograms{1,e,t,s}(:,:,r),meandayspectra{t});   
                zscorespectrograms{1,e,t,s}(:,:,r)=bsxfun(@rdivide,allwindowedspectrograms{1,e,t,s}(:,:,r),stddayspectra{t});
            end
        end
    end
end

for e=epochs
    for t=tetrodes
        continuousspectrograms{1,e,t}=bsxfun(@minus,continuousspectrograms{1,e,t},meandayspectra{t});   
        continuousspectrograms{1,e,t}=bsxfun(@rdivide,continuousspectrograms{1,e,t},stddayspectra{t});
    end
end


%% Eighth, calculates mean spectrogram for each state (and overall mean).

meanspectrograms_state=cell(1,tetno,3);        %% discards epoch, lumps all belonging to a state
meanspectrograms_allevents=cell(1,tetno);      %% discards state, lumps all events

for t=tetrodes       
    for s=1:3
        dummy=[];
        for e=epochs
            dummy=cat(3,dummy,zscorespectrograms{1,e,t,s});      
        end
        meanspectrograms_state{1,t,s}=mean(dummy,3);
        meanspectrograms_allevents{1,t}=cat(3,dummy,meanspectrograms_allevents{1,t});
    end
    meanspectrograms_allevents{1,t}=mean(meanspectrograms_allevents{1,t},3);
end
    
%% Final variable list.

eventtimes;
eegwindows;
statecounts;
allwindowedspectrograms;
meaneventspectrograms_epochlump;
continuousspectrograms;
meandayspectra;
stddayspectra;
zscorespectrograms;
meanspectrograms_state;
meanspectrograms_allevents;


%% Ninth, plot mean spectrogram for events.

for s=1:3
    figure
        
    for t=tetrodes
        h = subplot(1,14,t);
        imagesc(times,freqs,meanspectrograms_state{1,t,s}',[-0.2,3]);
            set(gca,'YDir','normal');
        if t==tetrodes(end)
            colorbar
        end
        string = sprintf('%d',t);
        title(string);

    end
    %%% TITLE code.
    eventcount = sum(statecounts(s,:));
    h = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
             1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 1,sprintf('Mean %s spectrogram, %s, Day %d, reftet=%d, state %d (%d events)', ...
                              eventtype,animalname,days,ref_tetrode,s,eventcount),'HorizontalAlignment'...
                ,'center','VerticalAlignment', 'top');
    
end



%% Tenth, plot some random individual event spectrograms.
    
    % set these manually to your taste
    e=2;                                 % epoch manually
    events_plot=3;                         % no. events to plot
    tetrodes_plot=1:14;                  % tetrodes to plot
    states_plot=1;
    %
    
    for s=states_plot
        for i=1:events_plot
            N=size(zscorespectrograms{1,e,t,s},3);        % # events to choose from
            eventno=1;    % ceil(N*rand);
            figure
            title(string);
            for t=tetrodes_plot
                subplot(1,14,t)
                imagesc(times,freqs,zscorespectrograms{1,e,t,s}(:,:,eventno)');
                set(gca,'YDir','normal');
                string = sprintf('%d',t);
                title(string);
            end
            %%% TITLE code.
            h = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
                1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
            text(0.5, 1,sprintf('Single %s event spectrogram, %s, Day %d, reftet=%d, state %d (event #%d)', ...
                eventtype,animalname,days,ref_tetrode,s,eventno),'HorizontalAlignment'...
                ,'center','VerticalAlignment', 'top');
        end
    end


%% Eleventh, plot the continuous spectrogram (random 1 minute periods).

%% to fix:
%% find random one minute interval
%% concatenate tetrodes into a single plot
%% add velocity plot

figure

for t=tetrodes
      subplot(tetrodes(end)+1,1,t)
        imagesc(times,freqs,[-0.2,5]);
            set(gca,'YDir','normal');
        string = sprintf('%d',t);
        title(string);
end
colorbar









