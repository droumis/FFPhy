% given a selected tetrode (e.g. CA1), calculates power spectrum time locked
% to event (startind)

%% THREE ingredients:
  %%  1. eeg struct     (all epochs for a day, use Mattias' loadeegstruct)
  %%  2. <event> struct  (after <event>dayprocess and <event>extract)
  %%  3. position struct     (NOT rawposition)

%% Set these parameters manually.
    % Chronux params are specified below (fourth block)
    

    
animalname = 'mit';                              % to label graphs later
eventtype = 'ripples';                              % must be the same name as the events' data structure
days = 13;                                           % days to analyze
ref_tetrode = 16;                                   % reference tetrode
tetrodes = [16 22];                                    % tetrodes to analyze
epochs = 1:2;                                         % epochs to analyze
    eegstruct=loadeegstruct('/data13/anna/Mitt/',animalname,'eeg',days,epochs,tetrodes);
velocity_state_1 = 0.1;                               % state 1 (immobile)
velocity_state_2 = 1;                               % state 2 (run)
windowsize_sec = .4;                               %                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            % window length
Fs = 1000;                                          % sample rate of EEG
halfwindow_samp = (windowsize_sec*Fs)/2;
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

eegs=nestcell(days(end),epno,tetno,3);

for d=days
for e=epochs
    
    epoch_starttime=eegstruct{d}{e}{tetrodes(1)}.starttime;
    epoch_times=epoch_starttime:1/Fs:(epoch_starttime+(1/Fs)*length(eegstruct{d}{e}{tetrodes(1)}.data));
    
    for t=tetrodes
        for r=1:length(eventtimes{d,e})                                   % (!!) disregards first and last 5 events
            centerindex=lookup(eventtimes{d,e}(r),epoch_times);     % for centering eeg window
            
            time_elapsed=(eventtimes{d,e}(r)-epoch_starttime)/1000;   % for probe style position, in SECONDS
            
            velocityindex=lookup(time_elapsed,position{e}.data(:,1)-position{e}.data(1,1));         %  added for state
            velocity=position{e}.data(velocityindex,5);
            s=0;
            if velocity<velocity_state_1           % state 1, immobile
                s=1;
            elseif velocity>velocity_state_2       % state 2, run
                s=2;
            else                                   % state 3, intermediate
                s=3;
            end
            eegs{d}{e}{t}{s}=cat(1,eegs{d}{e}{t}{s}, ...
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
            velocityindex=lookup(eventtimes{d,e}(r),position{e}.data(:,1));            %   added for state
            velocity=position{e}.data(velocityindex,5);
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


%% Fourth, calculate spectra.

mean_spectra=nestcell(days(end),tetno,3);           % holds all events' spectra

% chronux params
params.Fs=1000;
params.err = [1 0.05];   % theoretical error bars
params.trialave = 1;
params.fpass = [0 400];
params.tapers = [3 5];
sessions=days;
channels=tetrodes;

eegs_state=nestcell(sessions(end),channels(end));

%pool eegs
for d=sessions
    for c=channels
                    dummy=[];
        for s=1:3
            for e=epochs
                dummy=cat(1,dummy,eegs{d}{e}{c}{s});
            end

        end
      eegs_state{d}{c}=dummy;
    end
end

mean_spectra=nestcell(sessions(end),channels(end),3);
err_spectra=nestcell(sessions(end),channels(end),3);


% calculate mean spectra w/ error bars
for d=sessions
    for c=channels
            [mean_spectra{d}{c},frequencies,err_spectra{d}{c}]=mtspectrumc(double(eegs_state{d}{c}'),params);
        disp('spectrum')
    end
end

norm_spectra=nestcell(sessions(end),channels(end),3);

mitt_mean_spectra=mean_spectra;
mitt_err_spectra=err_spectra;

% calculate normalized mean spectra
for d=13
    for c=channels
            mitt_norm_spectra{d}{c}=mitt_mean_spectra{d}{c}/sum(mitt_mean_spectra{d}{c});
            mitt_norm_err_spectra{d}{c}=mitt_err_spectra{d}{c}/sum(mitt_mean_spectra{d}{c});
        disp('spectrum')
    end
end




for c=3
for c2=16
    figure
hold on
semilogx(frequencies,ann_norm_spectra{14}{c},'Color',[.7 .7 .7],'LineWidth',2)
semilogx(frequencies,ann_norm_err_spectra{14}{c},'Color',[.7 .7 .7],'LineWidth',2)
semilogx(frequencies,mitt_norm_spectra{13}{c2},'Color',[.9 .7 .7],'LineWidth',2)
semilogx(frequencies,mitt_norm_err_spectra{13}{c2},'Color',[.9 .7 .7],'LineWidth',2)
end
end
set(gca,'YScale','log')


for c=23
for c2=22
    figure
hold on
semilogx(frequencies,ann_norm_spectra{14}{c},'Color',[.7 .7 .7],'LineWidth',2)
semilogx(frequencies,ann_norm_err_spectra{14}{c},'Color',[.7 .7 .7],'LineWidth',2)
semilogx(frequencies,mitt_norm_spectra{13}{c2},'Color',[.9 .7 .7],'LineWidth',2)
semilogx(frequencies,mitt_norm_err_spectra{13}{c2},'Color',[.9 .7 .7],'LineWidth',2)
end
end
set(gca,'YScale','log')















