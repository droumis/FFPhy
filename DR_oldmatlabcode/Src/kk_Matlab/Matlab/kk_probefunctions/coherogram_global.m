%% power spectra: raw and in different states

%%  ingredients:
%    1. eegstruct
%    2. pos 

animalname = 'mit';                              % to label graphs later
sessions = 13;                                           % sessions to analyze
channels = [4 5 20 16];                                    % channels to analyze
epochs = 1:2;                                         % epochs to analyze

eegstruct=loadeegstruct('/data13/anna/Mitt/',animalname,'eeg',sessions,epochs,channels);
position;

velocity_state_1 = 0.1;                               % state 1 (immobile)
velocity_state_2 = 1;                               % state 2 (run)
windowsize=1;                                     % set to taste

%% First, classify windows into state 1 vs. 2 (and mixed states, 3)

windowstates=nestcell(sessions(end),epochs(end));       % classification of each window
statecounts=zeros(3,epochs(end));     % tabulation

for s=sessions      
    for e=epochs                                            
        epochvels=position{e}.data(:,5);
        epochlength=(position{e}.data(end,1)-position{e}.data(1,1));  % in seconds
        nowindows=floor(epochlength/windowsize);        % MANUALLY SET TO 45 MINUTES
        windowstates{s}{e}=nan(nowindows,1);
        postimes=position{e}.data(:,1)-position{e}.data(1,1);  % in seconds

        for w=1:nowindows
            starttime=(w-1)*windowsize;
                startindex=lookup(starttime,postimes);
            endindex=lookup(starttime+windowsize,postimes);
            meanvel=mean(epochvels(startindex:endindex));
            if meanvel<velocity_state_1                                               % state 1
                windowstates{s}{e}(w)=1;
                statecounts(1,e)=statecounts(1,e)+1;
            elseif meanvel>velocity_state_2                                           % state 2
                windowstates{s}{e}(w)=2;    
                statecounts(2,e)=statecounts(2,e)+1;
            else
                windowstates{s}{e}(w)=3;                                                % mixed state
                statecounts(3,e)=statecounts(3,e)+1;
            end
        end
    end
end
statecounts

%% Second, collect raw eegs into windows classified by state.

eegs=nestcell(sessions(end),epochs(end),channels(end),3);

for d=sessions
    for e=epochs
        for t=channels
            Fs=eegstruct{d}{e}{t}.samprate;
            data=eegstruct{d}{e}{t}.data;
            numwindows=length(windowstates{d}{e})-5;
               eegs{d}{e}{t}{1}=[];
               eegs{d}{e}{t}{2}=[];
               eegs{d}{e}{t}{3}=[];
            for w=1:numwindows                   %% iterate over each nonoverlapping window
                startindex=(w-1)*windowsize*Fs+1;
                datawin=data(startindex:(startindex+windowsize*Fs))';
                %to save precious time, initialize eegs to number of state
                %windows
                if windowstates{d}{e}(w)==1
                    eegs{d}{e}{t}{1}=[eegs{d}{e}{t}{1} ; datawin];
                elseif windowstates{d}{e}(w)==2
                    eegs{d}{e}{t}{2}=[eegs{d}{e}{t}{2} ; datawin];
                else
                    eegs{d}{e}{t}{3}=[eegs{d}{e}{t}{3} ; datawin];
                end
            end
            disp('eeg')
        end
    end
end


%% Third, calculate individual spectra.

% params.Fs=1000;
% params.err = 0;
% params.fpass = [0 400];
% params.tapers = [3 5];
% 
% spectra=nestcell(sessions(end),epochs(end),channels(end),3);
% 
% for d=sessions
%     for e=epochs
%         for c=channels
%             for s=1:3
%                 for w=1:size(eegs{d}{e}{c}{s},1)
%                     [spectra{d}{e}{c}{s}(w,:),frequencies]=mtspectrumc(double(eegs{d}{e}{c}{s}(w,:)),params);
%                 end
%             end
%             disp('spectrum')
%         end
%     end
% end

%% Pool eegs and then calculate coherograms.

params.Fs=1000;
params.err = [1 0.05];   % theoretical error bars
params.trialave = 1;
params.fpass = [0 400];
params.tapers = [3 5];
channelpairs = [16 4; 16 5 ; 16 20];     %%%%%%%%%%%

eegs_state=nestcell(sessions(end),channels(end),3);

%pool eegs
for d=sessions
    for c=channels
        for s=1:3
            dummy=[];
            for e=epochs
                dummy=cat(1,dummy,eegs{d}{e}{c}{s});
            end
            eegs_state{d}{c}{s}=dummy;
        end
    end
end

mean_cohgram=nestcell(sessions(end),size(channelpairs,1),3);

% calculate coherograms w/ error bars

for d=sessions
        for s=1:3
            for p=1:size(channelpairs,1)
            [mean_cohgram{d}{p}{s},~,~,~,~,frequencies,~,~] = ...
                    coherencyc(double(eegs_state{d}{channelpairs(p,1)}{s}'), ...
                               double(eegs_state{d}{channelpairs(p,2)}{s}'), ...
                               params);
            disp('coherogram')
            end
        end
end


%% Plot coherograms between animal.


    figure
    hold on
    
for p=[1]
        %semilogx(frequencies,ann_mean_cohgram{14}{p}{1},'Color',[.8 .8 .8],'LineWidth',2)
        semilogx(frequencies,ann_mean_cohgram{14}{p}{1},'Color',[0 0 0],'LineWidth',2)
       %semilogx(frequencies,ann_mean_cohgram{14}{p}{3},'Color',[.5 .5 .5],'LineWidth',2)

        %semilogx(frequencies,mitt_mean_cohgram{13}{p}{1},'Color',[.8 0 0],'LineWidth',2)

        %semilogx(frequencies,mitt_mean_cohgram{13}{p}{3},'Color',[.5 0 0],'LineWidth',2)
                set(gca,'XScale','log')
end

for p=1
        semilogx(frequencies,mitt_mean_cohgram{13}{p}{1},'Color',[1 0 0],'LineWidth',2)
end

    %title(num2str(channelpairs(p,:)));

% find difference peaks
    figure
        semilogx(frequencies,mitt_mean_cohgram{13}{1}{2}./ann_mean_cohgram{14}{1}{2},'Color',[1 0 0],'LineWidth',2)
        ylim([-.2 2])











