% given a selected tetrode (e.g. CA1), calculates spectrograms time locked
% to event (startind)

%% THREE ingredients:
  %%  1. eeg struct     (all epochs for a day, use Mattias' loadeegstruct)
  %%  2. <event> struct  (after <event>dayprocess and <event>extract)
  %%  3. pos struct     (NOT rawpos)

%% Set these parameters manually.
    % Chronux params are specified below (fourth block)
animalname = 'Egypt';                              % to label graphs later
eventtype = 'reward';                              % must be the same name as the events' data structure
days = 5;                                           % days to analyze
ref_tetrode = 0;                                   % reference tetrode
reward_bits = 6:8;
tetrodes=[1 2 3 6 7 10 11 13 15 16 17];
tetrode_pairs = [1 17 ; 2 17 ; 3 17; 6 17 ; 7 17 ; 10 17 ; 13 17 ; 15 17 ; 16 17];                                     % tetrodes to analyze
    tetno=max(max(tetrode_pairs));
    pairno=size(tetrode_pairs,1);
ref_tetrode_gd=11;
dereference=0;    % make 1 if want to add reference eeg back
epochs = [2 4 6];                                         % epochs to analyze
velocity_state_1 = 2;                               % state 1 (immobile)
velocity_state_2 = 8;                               % state 2 (run)
windowsize_sec = 20;                               %                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            % window length

%% calculated for you
halfwindow_samp = (windowsize_sec*1500)/2;
tetpairno=size(tetrode_pairs,1);

epno=epochs(end);
filename='';
rewardflag=0;

switch eventtype
    case 'reward'
        events=DIO;
        rewardflag=1;
    case 'errortrials'
        events=errortimes;
        rewardflag=2;
    case 'rewardaudio'
        events=rewardtimes;
        rewardflag=3;
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

if rewardflag == 1         %%% because reward is different than EEG events
    for d=days
        for e=epochs
            for w=reward_bits
                eventtimes{d,e}=[eventtimes{d,e}; events{d}{e}{w}.pulsetimes(:,1)/10000];
            end
            eventtimes{d,e}=sort(eventtimes{d,e});
        end
    end
elseif rewardflag==2          %error trials, well was triggered but no reward delivered
    for d=days
        for e=epochs
            eventtimes{d,e}=[eventtimes{d,e}; events{d,e}(:,3)/10000];
        end
    end
elseif rewardflag==3
    for d=days
        for e=epochs
            eventtimes{d,e}=[eventtimes{d,e}; events{d,e}(:,3)/10000];
        end
    end
else                    %% eeg events like ripples, gamma, etc.
    for d=days
        for e=epochs
            eventtimes{d,e}=[eventtimes{d,e}; events{d}{e}{ref_tetrode}.midtime];
        end
    end
end

%% Second, copy the raw eeg data into windows aligned to the event starttime.

eegwindows=cell(days(end),epno,tetno,3);

for d=days
for e=epochs    
    for t=tetrodes
    epoch_starttime=eegstruct{d}{e}{t}.starttime;          % time, in seconds, since Open All Files
    Fs=eegstruct{d}{e}{t}.samprate;
    epoch_times=epoch_starttime:1/Fs:(epoch_starttime+(1/Fs)*length(eegstruct{d}{e}{t}.data));

        for r=2:(length(eventtimes{d,e})-1)                                   % (!!) since so few error trials, DO NOT THROW AWAY EVENTS
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
            if dereference==1
            groundeeg=eegstruct{d}{e}{ref_tetrode_gd}.data((centerindex-halfwindow_samp):(centerindex+halfwindow_samp))';
            eegwindows{d,e,t,s}=cat(1,eegwindows{d,e,t,s}, ...
                groundeeg+eegstruct{d}{e}{t}.data((centerindex-halfwindow_samp):(centerindex+halfwindow_samp))');
            else
            eegwindows{d,e,t,s}=cat(1,eegwindows{d,e,t,s}, ...
                eegstruct{d}{e}{t}.data((centerindex-halfwindow_samp):(centerindex+halfwindow_samp))');                
            end
        end
    end
end
end

%% Third, tabulate the # of events in each of the 3 states.

statecounts=zeros(3,epno);

d=days
for e=epochs
        for r=(1):(length(eventtimes{d,e}))
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


%% Fourth, calculate all individual coherograms pairs.

allwindowedcohgrams=cell(days(end),epno,pairno);           % holds all events' spectograms

% chronux params
movingwin = [200 20]/1000;          
params.Fs = 1500;
params.err = [2 0.05];
params.fpass = [0 250];
params.tapers = [3 5];

for d=days
for e=epochs
    for p=1:size(tetrode_pairs,1)
        windows1=[];
        windows2=[];
        for s=1:3           % consolidate all eeg windows from all states for eeg
                windows1=[windows1 ; eegwindows{d,e,tetrode_pairs(p,1),s}];
                windows2=[windows2 ; eegwindows{d,e,tetrode_pairs(p,2),s}];
        end
        windows1=windows1';
        windows2=windows2';
        [C,~,~,~,~,times,freqs]=cohgramc(windows1,windows2,movingwin,params);
        allwindowedcohgrams{d,e,p}=C;   % adds the S variable (2D matrix) to the third matrix dimension
    end
end
end




% %% Fifth, calculates the continuous spectrogram, day mean, and std spectra for each tetrode.  (long calculation)
% 
% continuousspectrograms=cell(days(end),epno,tetno);
% meandayspectra=cell(days(end),tetno);
% stddayspectra=cell(days(end),tetno);
% 
% for d=days
%     for t=tetrodes                                
%         dummy=[];
%         for e=epochs
%             [S_full,junkt,junkf,junkserr] = mtspecgramc(eegstruct{d}{e}{t}.data,movingwin,params);
%             dummy=[dummy;S_full];        
%             continuousspectrograms{d,e,t}=S_full;
%         end
%     meandayspectra{t}=mean(dummy,1);
%     stddayspectra{t}=std(dummy,1);
%     end
% end
% 
% % z-score the continuous spectrogram
% 
% for d=days
% for e=epochs
%     for t=tetrodes
%         continuousspectrograms{d,e,t}=bsxfun(@minus,continuousspectrograms{d,e,t},meandayspectra{t});   
%         continuousspectrograms{d,e,t}=bsxfun(@rdivide,continuousspectrograms{d,e,t},stddayspectra{t});
%     end
% end
% end
% 
% 
% disp('.finished with computing spectrograms.')
% save(filename,'-v7.3');
% 
% 
% %% Sixth, z-scores the individual spectrograms ('zscorespectrograms').
% 
% zscorespectrograms=cell(days(end),epno,tetno,3);
% 
% for d=days
%     for e=epochs
%         for t=tetrodes
%             for s=1:3
%                 if ~isempty(allwindowedcohgrams{d,e,t,s})       % in case there are no events
%                     for r=1:size(allwindowedcohgrams{d,e,t,s},3)
%                         zscorespectrograms{d,e,t,s}(:,:,r)=bsxfun(@minus,allwindowedcohgrams{d,e,t,s}(:,:,r),meandayspectra{t});
%                         zscorespectrograms{d,e,t,s}(:,:,r)=bsxfun(@rdivide,allwindowedcohgrams{d,e,t,s}(:,:,r),stddayspectra{t});
%                     end
%                 end
%             end
%         end
%     end
% end



% %% Seventh, calculates mean coherograms: (1) state, (2) epochs (3) overall.
% 
% meancohgrams_state=cell(days(end),pairno,3);             %% discards epoch, pools all within state
% meancohgrams_epoch=cell(days(end),epno,pairno);          %% discards states, pools all within epoch
% meancohgrams_allevents=cell(days(end),pairno);           %% pools all events
% 
% epochs=[2 4 6];
% 
% for d=days
%     for p=size(tetrode_pair,1)
%         
%         for s=1:3
%             dummy=[];
%             for e=epochs
%                 if ~isempty(allwindowedcohgrams{d,e,p,s})
%                     dummy=cat(3,dummy,allwindowedcohgrams{d,e,p,s});
%                 end
%             end
%             meancohgrams_state{d,t,s}=mean(dummy,3);
%             meancohgrams_allevents{d,t}=cat(3,dummy,meancohgrams_allevents{d,p});
%         end
%         
%         meancohgrams_allevents{d,p}=mean(meancohgrams_allevents{d,p},3);
%         
%         for e=epochs
%             dummy=[];
%             for s=1:3
%                 if ~isempty(allwindowedcohgrams{d,e,p,s})
%                     dummy=cat(3,dummy,allwindowedcohgrams{d,e,p,s});
%                 end
%             end
%             meancohgrams_epoch{d,e,t}=mean(dummy,3);
%         end
%         
%     end
% end


    
%% Final variable list.

% eventtimes;
% eegwindows;
% statecounts;
% allwindowedcohgrams;
% %continuousspectrograms;
% meandayspectra;
% stddayspectra;
% zscorespectrograms;
% meanspectrograms_state;
% meanspectrograms_epoch;
% meanspectrograms_allevents;


%% Ninth, plot.

for p=1:size(tetrode_pairs,1)
    figure
    imagesc(times,freqs,mean(allwindowedcohgrams{d,e,p},3)');
    set(gca,'YDir','normal');
    colorbar
%%% TITLE code.
eventcount = sum(statecounts(s,:));
h = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,sprintf('Mean %s coherogram, %s, Day %d, %d - %d', ...
    eventtype,animalname,d,tetrode_pairs(p,1),tetrode_pairs(p,2)),'HorizontalAlignment'...
    ,'center','VerticalAlignment', 'top');
end



