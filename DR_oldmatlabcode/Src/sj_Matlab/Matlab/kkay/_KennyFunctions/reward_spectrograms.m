% given a selected tetrode (e.g. CA1), calculates spectrograms time locked
% to event (startind)

%% Three ingredients:
  %%  1. eeg struct     (all epochs for a day, use Mattias' loadeegstruct)
  %%  2. dio struct  (after <event>dayprocess and <event>extract)
  %%  3. pos struct     (NOT rawpos)

%% Set these parameters manually.
    % Chronux params are specified below (fourth block)
animalname = 'Coriander';                              % to label graphs later
eventtype = 'reward';                              % must be the same name as the events' data structure
reward_bits = 7:8;
days = 5;                                           % days to analyze
tetrodes = 8:14;                                    % tetrodes to analyze
epochs = 6;                                         % epochs to analyze
windowsize_sec = 10;                                 % window length
velocity_state_1 = 2;
velocity_state_2 = 8;
Fs = 1500;                                          % sample rate of EEG
halfwindow_samp = (windowsize_sec*1500)/2;
tetno=tetrodes(end);
epno=epochs(end);
filename='';

switch eventtype
    case 'reward'
        events=DIO;
end

%% First, collects all event timestamps from chosen ref. tetrode.

eventtimes=cell(days(end),epno);

for d=days
for e=epochs
    for b=reward_bits
        eventtimes{d,e}=[eventtimes{d,e}; events{d}{e}{b}.pulsetimes(:,1)/10000];
    end
    eventtimes{d,e}=sort(eventtimes{d,e});
end
end

%% Second, copy the raw eeg data into windows aligned to the reward time.
%% In parallel, copy velocity position data into another structure.

eegwindows=cell(days(end),epno,tetno,3);
velocitywindows=cell(days(end),epno,tetno,3);

for d=days
for e=epochs
    epoch_starttime=eegstruct{d}{e}{14}.starttime;
    epoch_times=epoch_starttime:1/Fs:(epoch_starttime+(1/Fs)*length(eegstruct{d}{e}{14}.data));
    for t=tetrodes
        for r=2:(length(eventtimes{d,e})-1)                                   % (!!) disregards first and last event
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
        for r=(1+5):(length(eventtimes{d,e})-5)
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

epochs=[6];

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


%% Ninth, plot overall mean spectrogram.


figure
  for d=days
        for t=tetrodes
            
            h = subplot(8,1,t-6);
            imagesc(times,freqs(1:40),meanspectrograms_allevents{d,t}(:,1:40)',[0,2]);
            colormap('jet')
            set(gca,'YDir','normal');
            string = sprintf('%d',t);
            title(string);
        end
        
        %%% TITLE code.
        eventcount = sum(statecounts(s,:));
        h = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
            1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.5, 1,sprintf('Mean %s spectrogram, %s, Day %d, overall', ...
            eventtype,animalname,d,s),'HorizontalAlignment'...
            ,'center','VerticalAlignment', 'top');
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
    states_plot=2;
    %
    
    for s=states_plot
        for i=1:size(zscorespectrograms{d,e,t,s},3)
            N=size(zscorespectrograms{d,e,t,s},3);        % # events to choose from
            eventno=ceil(N*rand);
            figure
            title(string);
            for t=tetrodes_plot
                subplot(8,1,t-7)
                imagesc(times,freqs,zscorespectrograms{d,e,t,s}(:,:,i)',[-0.2,3]);
                colormap('bone')
                set(gca,'YDir','normal');
                string = sprintf('%d',t);
                title(string);
            end
            %%% TITLE code.
            h = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
                1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
            text(0.5, 1,sprintf('Single %s event spectrogram, %s, Day %d, (event #%d)', ...
                eventtype,animalname,d,i),'HorizontalAlignment'...
                ,'center','VerticalAlignment', 'top');
        end
    end


%