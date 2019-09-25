
%% calculates gamma power correlations in 100 ms non-overlapping windows, by state (1,2,3)

%% 2 ingredients 
    % 1. <rhythm>struct   (!!) do not use for theta or delta (downsampled)
    % 2. pos struct

% parameters for analysis:

animalname = 'Bashir';                              % to label graphs later
rhythm = 'ripples';
windowsize=.1;                                    % in seconds
ref_tetrode = 14;                                   % ** pick reference
pos_Fs=29.97003;                                    % ETSC video standard
    windowsize_possamp_round=round(windowsize*pos_Fs);
days = 4;                                           % days to analyze
tetrodes = 8:14;                                    % tetrodes to analyze
epochs = 3;                                         % epochs to analyze
velocity_state_1 = 2;                               % state 1 (immobile)
velocity_state_2 = 8;                               % state 2 (run)


    switch rhythm
        case 'ripples'
            rhythmstruct=ripple;
        case 'gammal'
            rhythmstruct=lowgamma;
        case 'gammah'
            rhythmstruct=highgamma;
        case 'gammaff'
            rhythmstruct=fastfastgamma;
        case 'slowripples'
            rhythmstruct=slowripple;
        case 'hightheta'
            rhythmstruct=hightheta;
        otherwise
            disp('this event type is not recognized..');
    end

Fs=rhythmstruct{days(end)}{epochs(end)}{tetrodes(end)}.samprate;   % probably just 1500
    windowsize_samp=round(windowsize*Fs);
    
    
%% First, classify and tabulate windows into state 1, 2, and 3.

windowstates=cell(days(end),epochs(end));  % all windows in epoch
statecounts=zeros(3,epochs(end));  % tabulation, just for bookkeeping

for d=days        % choose which day
    for e=epochs                                             % 7 epochs / day
        nowindows=floor(length(pos{d}{e}.data(:,5))/windowsize_possamp_round);
        velocities=pos{d}{e}.data(:,5);
        windowstates{d,e}=nan(nowindows,1);
        for w=1:(nowindows-5)      % ignore the last several 25 ms windows
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

%% Second, collect windowed power values (segregated by state)

power_state=cell(days(end),epochs(end),tetrodes(end),3);

for d=days
    for e=epochs
        eeg_startindex=lookup(pos{d}{e}.data(1,1),(rhythmstruct{d}{e}{tetrodes(end)}.starttime):1/Fs: ...
            (rhythmstruct{d}{e}{tetrodes(end)}.starttime+length(rhythmstruct{d}{e}{tetrodes(end)}.data)/Fs));
        for t=tetrodes
            for w=1:length(windowstates{d,e}-5)                       % iterate over each window
                true_time=(w-1)*windowsize_possamp_round/pos_Fs;    % actual time at beginning of window
                start_index=eeg_startindex+round(true_time*Fs);
                datawindow = rhythmstruct{d}{e}{t}.data(start_index:(start_index+windowsize_samp),1);
                rms = (mean(datawindow.^2)).^0.5;                   % root-mean-square power
                if windowstates{d,e}(w)==1
                    power_state{d,e,t,1}=[power_state{d,e,t,1} ; rms];
                elseif windowstates{d,e}(w)==2
                    power_state{d,e,t,2}=[power_state{d,e,t,2} ; rms];
                else
                    power_state{d,e,t,3}=[power_state{d,e,t,3} ; rms];
                end
            end
        end
    end
end
                      

%% Variable list.

windowstates;                 % running state classification of all windows
statecounts;                  % matrix tabulation of state 1, 2, 3
power_state;


%% Third, calculate r values for each of 3 states, given reference.

% lump over all epochs

power_state_epochlump = cell(days(end),tetrodes(end),3);   % mean spectrogram for all events from each state

for d=days
    for t=tetrodes
        for s=1:3
            dummy=[];
            for e=epochs
                dummy=[dummy ; power_state{d,e,t,s}];   % concatenate events over epoch
            end
            power_state_epochlump{d,t,s}=dummy;
        end
    end
end

% calculate r

r_values = nan(tetrodes(end),3);
d=4;

for s=1:3
    for t=tetrodes
        R=corrcoef(power_state_epochlump{d,t,s},power_state_epochlump{d,ref_tetrode,s});
        r_values(t,s)=R(1,2);
    end
end
 
figure
bar(8:14,r_values(8:14,1:2))       
colormap('bone')
        %%% TITLE code.
        h = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
            1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
       
        text(0.5, 1,sprintf('%s power correlations, %s, windowsize=%d, Day %d,epochs: %s, ref_tetrode=%d, state1=%d, state2=%d', ...
            rhythm,animalname,windowsize,d,strtrim(sprintf('%d ',epochs)),ref_tetrode,velocity_state_1,velocity_state_2), ...
            'HorizontalAlignment','center','VerticalAlignment', 'top');











