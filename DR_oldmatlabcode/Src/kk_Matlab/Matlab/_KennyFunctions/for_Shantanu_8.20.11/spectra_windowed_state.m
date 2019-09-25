
%%  classify states based on velocity
%%  non-overlapping window segmentation, state 1 vs. state 2 spectra

%% used for Bashir analysis on day 4, 6.14.12

%% 2 ingredients 
    % 1. eegstruct
    % 2. pos struct

% parameters for analysis:

animalname = 'Bashir';                              % to label graphs later
windowsize=1;                                     % in seconds
pos_Fs=29.97003;                                    % ETSC video standard
    windowsize_possamp_round=round(windowsize*pos_Fs);
days = 4;                                           % days to analyze
tetrodes = 8:14;                                    % tetrodes to analyze
epochs = [2];                                         % epochs to analyze
velocity_state_1 = 2;                               % state 1 (immobile)
velocity_state_2 = 8;                               % state 2 (run)
Fs=eegstruct{days(end)}{epochs(end)}{tetrodes(end)}.samprate;   % probably just 1500
    windowsize_samp=round(windowsize*Fs);

%% First, classify and tabulate windows into state 1, 2, and 3.

windowstates=cell(days(end),epochs(end));           % all windows in epoch
statecounts=zeros(3,epochs(end));                   % tabulation, just for bookkeeping

for d=days                                          % choose which day
    for e=epochs                                             % 7 epochs / day
        nowindows=floor(length(pos{d}{e}.data(:,5))/windowsize_possamp_round);
        velocities=pos{d}{e}.data(:,5);
        windowstates{d,e}=nan(nowindows,1);
        for w=1:nowindows
            meanvel=mean(velocities(((w-1)*windowsize_possamp_round+1):(w*windowsize_possamp_round+1)));
            if meanvel<velocity_state_1                                             % state 1
                windowstates{d,e}(w)=1;
                statecounts(1,e)=statecounts(1,e)+1;
            elseif meanvel>velocity_state_2                                         % state 2
                windowstates{d,e}(w)=2;    
                statecounts(2,e)=statecounts(2,e)+1;
            else
                windowstates{d,e}(w)=3;                                                % mixed state
                statecounts(3,e)=statecounts(3,e)+1;
            end
        end
    end
end

statecounts

%% Second, collect raw eeg windows by state.

eeg_state=cell(days(end),epochs(end),tetrodes(end),3);

for d=days
    for e=epochs
        eeg_startindex=lookup(pos{d}{e}.data(1,1),(eegstruct{d}{e}{14}.starttime):1/Fs: ...
                                        (eegstruct{d}{e}{14}.starttime+length(eegstruct{d}{e}{14}.data)/Fs));
        for t=tetrodes
            for w=1:(length(windowstates{d,e})-5)                   %% iterate over each nonoverlapping window
                true_time=(w-1)*windowsize_possamp_round/pos_Fs;    %% actual time at beginning of window
                start_index=eeg_startindex+round(true_time*Fs);
                if windowstates{d,e}(w)==1
                    eeg_state{d,e,t,1}=[eeg_state{d,e,t,1} ; eegstruct{d}{e}{t}.data(start_index:  ...  % concatenates rows
                        (start_index+windowsize_samp))'];
                elseif windowstates{d,e}(w)==2
                    eeg_state{d,e,t,2}=[eeg_state{d,e,t,2} ; eegstruct{d}{e}{t}.data(start_index:  ...  % concatenates rows
                        (start_index+windowsize_samp))'];
                else
                    eeg_state{d,e,t,3}=[eeg_state{d,e,t,3} ; eegstruct{d}{e}{t}.data(start_index:  ...  % concatenates rows
                        (start_index+windowsize_samp))'];
                end
            end
        end
    end
end
                      

%% Third, calculate spectra.

params.Fs = Fs;
params.err = [0 0.05];
params.fpass = [0 400];
params.tapers = [3 5];

spectra_state = cell(days(end),epochs(end),tetrodes(end),3); 
                                        
for d=days
    for e=epochs
        for t=tetrodes
            for s=1:3
                for w=1:size(eeg_state{d,e,t,s},1)
                    [spectra_state{d,e,t,s}(w,:),frequencies]=mtspectrumc(eeg_state{d,e,t,s}(w,:),params);
                end
            end
        end
    end
end


%% Fourth, z-score spectra to all of the day's epochs.

spectra_daymeans=cell(days(end),tetrodes(end));
spectra_daystd=cell(days(end),tetrodes(end));

for d=days
    for t=tetrodes
        dummy=[];
        for s=1:3               %% lump spectra temporarily
            for e=epochs
                dummy=vertcat(dummy,spectra_state{d,e,t,s});
            end
        end
        spectra_daymeans{d,t}=mean(dummy,1);
        spectra_daystd{d,t}=std(dummy,1);
    end
end
    
spectra_state_zscored = cell(days(end),epochs(end),tetrodes(end),3);   

for d=days
    for s=1:3
        for e=epochs
            for t=tetrodes
                for w=1:size(eeg_state{d,e,t,s},1)
                    spectra_state_zscored{d,e,t,s}(w,:)=(spectra_state{d,e,t,s}(w,:)-spectra_daymeans{d,t})./spectra_daystd{d,t};
                end
            end                                             
        end
    end
end

%% Variable list.

windowstates;                 % running state classification of all windows
statecounts;                  % matrix tabulation of state 1, 2, 3
eeg_state;                    % raw data
spectra_state;                % raw spectra
spectra_daymeans;                  % day mean spectra
spectra_daystd;                   % day std spectra
spectra_state_zscored;        % z-scored spectra



%% Fifth, calculate mean spectra over all windows for each of 3 states.

spectra_state_zscored_mean=cell(days(end),tetrodes(end),3);

for t=tetrodes
    for s=1:3
        dummy=[];
        for e=epochs              %% lump spectra temporarily
            dummy=vertcat(dummy,spectra_state_zscored{d,e,t,s});
        end       
        spectra_state_zscored_mean{d,t,s}=mean(dummy,1);
    end
end

%% Sixth, calculate mean spectra for each state from RAW spectra.

spectra_state_mean=cell(days(end),tetrodes(end),4);    %% let state "4" be all of 2 and 3

for t=tetrodes
    for s=1:3
        dummy=[];
        for e=epochs         
            dummy=vertcat(dummy,spectra_state{d,e,t,s});
        end       
        spectra_state_mean{d,t,s}=mean(dummy,1);
    end
    for s=2:3
        dummy=[];
        for e=epochs            
            dummy=vertcat(dummy,spectra_state{d,e,t,s});
        end       
        spectra_state_mean{d,t,4}=mean(dummy,1);
    end
end




%% Seventh, plot for Bashir. 


figure
hold on

for t=8:14
    if t==9 
        semilogx(frequencies,spectra_state_mean{d,t,2}./spectra_state_mean{d,t,1},'g','LineWidth',2)
    elseif t==13 || t==8
        semilogx(frequencies,spectra_state_mean{d,t,2}./spectra_state_mean{d,t,1},'LineWidth',2,'Color',[0 0.5 0])
    elseif t==12
        semilogx(frequencies,spectra_state_mean{d,t,2}./spectra_state_mean{d,t,1},'r','LineWidth',2)
    elseif t==10
        semilogx(frequencies,spectra_state_mean{d,t,2}./spectra_state_mean{d,t,1},'m','LineWidth',2)
    else
        semilogx(frequencies,spectra_state_mean{d,t,2}./spectra_state_mean{d,t,1},'k','LineWidth',2)        
    end
    epoch_string=strtrim(sprintf('%d ',epochs));
    title(sprintf('Power spectra ratio, %s, Day %d, epochs %s, windowsize=%d, state1vs2 (%d,%d)',...
                              animalname,d,epoch_string,windowsize,velocity_state_1,velocity_state_2),'HorizontalAlignment'...
        ,'center','VerticalAlignment', 'top');
end
ylim([-0.5 2])


%% Seventh, plot for Coriander.


figure
hold on

for t=14:23
    if t==14 || t==15
        semilogx(frequencies,spectra_state_mean{d,t,2}./spectra_state_mean{d,t,1},'r','LineWidth',2)
    else
        semilogx(frequencies,spectra_state_mean{d,t,2}./spectra_state_mean{d,t,1},'k','LineWidth',2)        
    end
        %%% TITLE code.
    epoch_string=strtrim(sprintf('%d ',epochs));
    title(sprintf('Power spectra ratio, %s, Day %d, epochs %s, windowsize=%d, state1vs2 (%d,%d)',...
                              animalname,d,epoch_string,windowsize,velocity_state_1,velocity_state_2),'HorizontalAlignment'...
        ,'center','VerticalAlignment', 'top');
end
ylim([-0.5 2])


    