% takes ratio of POST to PRE

%% THREE ingredients:
  %%  1. eeg struct     (all epochs for a day, use Mattias' loadeegstruct)
  %%  2. dio struct  (after <event>dayprocess and <event>extract)

%% Set these parameters manually.
    % Chronux params are specified below (fourth block)
animalname = 'Bashir';                              % to label graphs later
eventtype = 'reward';                              % must be the same name as the events' data structure
reward_bits = 7:8;
days = 4;                                           % days to analyze
tetrodes = 8:14;                                    % tetrodes to analyze
epochs = 2;                                         % epochs to analyze

windowsize_sec = 0.5;                                 % window length

no_pre_windows = 4;
no_post_windows = 4;
no_skip_post=2;

velocity_state_1 = 2;
velocity_state_2 = 8;

Fs = 1500;                                          % sample rate of EEG
windowsize_samp = round(windowsize_sec*Fs);

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

eegwindows=cell(days(end),epno,tetno,2);                %% let 1 be PRE, let 2 be POST
velocitywindows=cell(days(end),epno,tetno,2);

for d=days
for e=epochs
    epoch_starttime=eegstruct{d}{e}{14}.starttime;
    epoch_times=epoch_starttime:1/Fs:(epoch_starttime+(1/Fs)*length(eegstruct{d}{e}{14}.data));
    for t=tetrodes
        for r=2:(length(eventtimes{d,e})-1)                                   % (!!) disregards first and last event
            centerindex=lookup(eventtimes{d,e}(r),epoch_times);
            
            % PRE windows
            for i=1:no_pre_windows
                eegwindows{d,e,t,1}=cat(1,eegwindows{d,e,t,1}, ...
                    eegstruct{d}{e}{t}.data((centerindex-i*windowsize_samp):(centerindex-(i-1)*windowsize_samp))');
            end
            
            % POST windows
            for i=1:no_post_windows
                eegwindows{d,e,t,2}=cat(1,eegwindows{d,e,t,2}, ...
                    eegstruct{d}{e}{t}.data((1+centerindex+(i-1+no_skip_post)*windowsize_samp): ...
                                            (1+centerindex+(i+no_skip_post)*windowsize_samp))');
            end
            
%             % state (unused for now)
%             velocityindex=lookup(eventtimes{d,e}(r),pos{d}{e}.data(:,1));         %  added for state
%             velocity=pos{d}{e}.data(velocityindex,5);
%             s=0;
%             if velocity<velocity_state_1           % state 1, immobile
%                 s=1;
%             elseif velocity>velocity_state_2       % state 2, run
%                 s=2;
%             else                                   % state 3, intermediate
%                 s=3;
%             end
%             eegwindows{d,e,t,s}=cat(1,eegwindows{d,e,t,s}, ...
%                 eegstruct{d}{e}{t}.data((centerindex-halfwindow_samp):(centerindex+halfwindow_samp))');
            
            
        end
    end
end
end



%% Third, calculate the individual spectrograms.

spectra=cell(days(end),epno,tetno,2);           % all raw spectra

% chronux params       
params.Fs = Fs;
params.err = [0 0.05];
params.fpass = [0 400];
params.tapers = [3 5];

for d=days
for e=epochs
    for t=tetrodes
        for s=1:2
            flag=0;
            for r=1:size(eegwindows{d,e,t,s},1)
                [S,freqs]=mtspectrumc(eegwindows{d,e,t,s}(r,:),params);
%                 if flag==0                              % this if clause initializes the 3D matrix of all spectrograms for a given tetrode
%                     allwindowedspectrograms{d,e,t,s}=nan(size(S,1),size(S,2),size(eegwindows{d,e,t,s},1));
%                     flag=1;
%                 end
                spectra{d,e,t,s}(r,:)=S;   % adds the S variable (2D matrix) to the third matrix dimension
            end
        end
    end
end
end




%% Fourth, calculate means of PRE and POST spectra and then calculate POST:PRE ratio.

meanspectra=cell(days(end),tetno,2);          
ratiospectra=cell(days(end),tetno);


for d=days
    for t=tetrodes
        
        for s=1:2
            dummy=[];
            for e=epochs
                if ~isempty(spectra{d,e,t,s})
                    dummy=[dummy ; spectra{d,e,t,s}];
                end
            end
            meanspectra{d,t,s}=mean(dummy,1);
        end
        ratiospectra{d,t}=meanspectra{d,t,2}./meanspectra{d,t,1};
    end
end





    
%% Final variable list.

eventtimes;
eegwindows;
freqs;
spectra;
meanspectra;
ratiospectra;


%% Seventh, plot for Bashir. 


figure
hold on

for t=8:14
    if t==9 
        semilogx(freqs,ratiospectra{d,t},'g','LineWidth',2)
    elseif t==13 || t==8
        semilogx(freqs,ratiospectra{d,t},'LineWidth',2,'Color',[0 0.5 0])
    elseif t==12
        semilogx(freqs,ratiospectra{d,t},'r','LineWidth',2)
    elseif t==10
        semilogx(freqs,ratiospectra{d,t},'m','LineWidth',2)
    else
        semilogx(freqs,ratiospectra{d,t},'k','LineWidth',2)        
    end
    epoch_string=strtrim(sprintf('%d ',epochs));
    title(sprintf('Ratio spectra, PRE- vs. POST- reward, %s, Day %d, epochs %s, windowsize_sec=%d',...
                              animalname,d,epoch_string,windowsize_sec),'HorizontalAlignment'...
        ,'center','VerticalAlignment', 'top');
end
ylim([0 2])

