
%% calculates gamma power correlations in 50 ms non-overlapping windows, by state (1,2,3)

%% 2 ingredients 
    % 1. <rhythm>struct   (!!) do not use for theta or delta (downsampled)

% parameters for analysis:

animalname = 'Bashir';                               % to label graphs later
rhythm = 'gammah';
windowsize=.050;                                     % in seconds
ref_tetrode = 9;
% pos_Fs=29.97003;                                    % ETSC video standard
% windowsize_possamp_round=round(windowsize*pos_Fs);
days=5;                                              % days to analyze
tetrodes = 1:14;                                     % tetrodes to analyze
epochs = 1:6;
run_epochs = 2;                                      % epochs to analyze
%velocity_state_1 = 2;                               % state 1 (immobile)
%velocity_state_2 = 8;                               % state 2 (run)


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
    

%% First, collect windowed power values (segregated by state)

power_sleep_run=cell(days(end),tetrodes(end),2);   % note that this lumps run + sleep epochs together: sleep is 1, run is 2

for d=days
    for t=tetrodes
        dummyrun=[];
        dummysleep=[];
        for e=epochs         % RUN
            if sum(e==run_epochs)
                dummyrun=[dummyrun ; rhythmstruct{d}{e}{t}.data(:,1)];
            else             % SLEEP
                dummysleep=[dummysleep ; rhythmstruct{d}{e}{t}.data(:,1)];
            end
        end
        
        for w=1:floor(length(dummyrun)/windowsize_samp)            % RUN
            start_index=1+(w-1)*windowsize_samp;
            datawindow = dummyrun(start_index:(start_index+windowsize_samp));
            rms = (mean(datawindow.^2)).^0.5;                   % root-mean-square power
            power_sleep_run{d,t,2}=[power_sleep_run{d,t,2} ; rms];
        end
        
        for w=1:floor(length(dummysleep)/windowsize_samp)            % SLEEP
            start_index=1+(w-1)*windowsize_samp;
            datawindow = dummysleep(start_index:(start_index+windowsize_samp));
            rms = (mean(datawindow.^2)).^0.5;                   % root-mean-square power
            power_sleep_run{d,t,1}=[power_sleep_run{d,t,1} ; rms];
        end
    end
end
                      

%% Variable list.

power_sleep_run;


%% Second, calculate r values for RUN vs. SLEEP, given reference.

% Calculate r.

r_values = nan(tetrodes(end),2);

for d=days
for p=1:2           % SLEEP: 1   //   RUN: 2
    for t=tetrodes
        R=corrcoef(power_sleep_run{d,t,p},power_sleep_run{d,ref_tetrode,p});
        r_values(t,p)=R(1,2);
    end
end
end

figure
bar3(r_values)
title(sprintf('Power correlations, sleep vs. run, in %s, %s, Day %d, reftet=%d', ...
    rhythm,animalname,d,ref_tetrode),'HorizontalAlignment'...
    ,'center','VerticalAlignment', 'top');

figure
imagesc(r_values,[0.3 1])
title(sprintf('Power correlations, sleep vs. run, in %s, %s, Day %d, reftet=%d', ...
    rhythm,animalname,d,ref_tetrode),'HorizontalAlignment'...
    ,'center','VerticalAlignment', 'top');
colormap('bone')
colorbar







