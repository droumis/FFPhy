%%  non-overlapping window segmentation, state 1 vs. state 2 spectra

%% used for RO1 grant analysis for Grendel, Maggie's animal

% load pos file (from a given day), which should contain velocity data for all epochs for a
% given day

windowsize=2;           % set to taste
possamplefreq=29.97;    % ETSC video standard
windowsize_samp=floor(2*29.97);

% classify windows into state 1 vs. 2 (and mixed states, 0)

windowstates=cell(1,7);  % epochs
statecounts=zeros(3,7);     % 3 states over 7 epochs -- this keeps track of number of windows for each

for d=3        % choose which day
    for e=1:7                                             % 7 epochs / day
        epochvels=pos{d}{e}.data(:,5);
        nowindows=floor(size(epochvels,1)/windowsize_samp);
        windowstates{e}=nan(nowindows,1);
        for w=1:nowindows
            meanvel=mean(epochvels(((w-1)*windowsize_samp+1):(w*windowsize_samp)));
            if meanvel<1                                                            % state 1
                windowstates{e}(w)=1;
                statecounts(1,e)=statecounts(1,e)+1;
            elseif meanvel>3                                                          % state 2
                windowstates{e}(w)=2;    
                statecounts(2,e)=statecounts(2,e)+1;
            else
                windowstates{e}(w)=0;                                                % mixed state
                statecounts(3,e)=statecounts(3,e)+1;
            end
        end
    end
end

% collect raw eegs  (recopying windows into classified states)

state1eeg=cell(1,7,14);
state2eeg=cell(1,7,14);
state3eeg=cell(1,7,14);

eegsamprate=1500;
eegsampwin=windowsize*eegsamprate;

for d=3
for e=1:7
    for t=[5 6 11 12 13]
        for w=1:(size(windowstates{e},1)-15)   %% iterate over each nonoverlapping window
            if windowstates{e}(w)==1
                state1eeg{1,e,t}=[state1eeg{1,e,t} ; eeg{d}{e}{t}.data((eegsampwin*(w-1)+1):  ...  % concatenates rows
                                                                     (eegsampwin*w))'];
            elseif windowstates{e}(w)==2
                state2eeg{1,e,t}=[state2eeg{1,e,t} ; eeg{d}{e}{t}.data((eegsampwin*(w-1)+1):  ...  % concatenates rows
                                                                     (eegsampwin*w))'];
            else
                state3eeg{1,e,t}=[state3eeg{1,e,t} ; eeg{d}{e}{t}.data((eegsampwin*(w-1)+1):  ...  % concatenates rows
                                                                     (eegsampwin*w))'];
            end
        end
    end
end
end
                      
alleeg(:,:,:,1)=state1eeg;         % consolidate into a single cell array, with state added
alleeg(:,:,:,2)=state2eeg;
alleeg(:,:,:,3)=state3eeg;

%% calculate spectra

params.Fs = 1500;
params.err = [2 0.05];
params.fpass = [0 400];
params.tapers = [3 5];

allspectra = cell(1,7,14,3);    % initialize
                                                         
    for s=1:3
        for e=1:7
            for t=1:14
                for w=1:size(alleeg{1,e,t,s},1)
                    [allspectra{1,e,t,s}(w,:),f,Serr]=mtspectrumc(alleeg{1,e,t,s}(w,:),params);
                end
            end                                             
        end
    end


%% z-score (and ratio) spectra to all of the day's epochs (for each of 14 tetrodes)


meanspectra=cell(1,14);
stdspectra=cell(1,14);
for t=1:14
    dummy=[];
    for s=1:3               %% lump spectra temporarily
        for e=1:7
            dummy=vertcat(dummy,allspectra{1,e,t,s});
        end
    end           
    meanspectra{t}=mean(dummy,1);
    stdspectra{t}=std(dummy,1);
end

    
    
zscorespectra = cell(1,7,14,3);    % initialize
ratiospectra = cell(1,7,14,3);    % initialize
 

    for s=1:3
        for e=1:7
            for t=1:14
                for w=1:size(alleeg{1,e,t,s},1)
                    zscorespectra{1,e,t,s}(w,:)=(allspectra{1,e,t,s}(w,:)-meanspectra{t})./stdspectra{t};
                    ratiospectra{1,e,t,s}(w,:)=(allspectra{1,e,t,s}(w,:)./meanspectra{t});
                end
            end                                             
        end
    end


% now take means over all windows for each of 3 states

mean_zscorespectra=cell(1,14,3);
for t=1:14
    for s=1:3
        dummy=[];
        for e=1:7              %% lump spectra temporarily
            dummy=vertcat(dummy,zscorespectra{1,e,t,s});
        end       
        mean_zscorespectra{1,t,s}=mean(dummy,1);
    end
end


for t=[6 11 12]
figure
hold on
plot(f,mean_zscorespectra{1,t,1})
plot(f,mean_zscorespectra{1,t,2},'r')
plot(f,mean_zscorespectra{1,t,3},'g')
s=sprintf('mean z-scored spectra, tetrode %d',t);
title(s);
set(gca,'xscale','log')
end


% take means
mean_ratiospectra=cell(1,14,3);

for t=1:14
    for s=1:3
        dummy=[];
        for e=1:7              %% lump spectra temporarily
            dummy=vertcat(dummy,ratiospectra{1,e,t,s});
        end       
        mean_ratiospectra{1,t,s}=mean(dummy,1);
    end
end

% plot
for t=[6 11 12]
figure
hold on
semilogx(f,mean_ratiospectra{1,t,1})
semilogx(f,mean_ratiospectra{1,t,2},'r')
semilogx(f,mean_ratiospectra{1,t,3},'g')
string = sprintf('Grendel day 4, <1 vs. >3, ratio spectra to day mean, tetrode %d',t);
title(string);
set(gca,'xscale','log')
end

figure
ratio_moving_stopped=(mean_ratiospectra{1,t,2}-mean_ratiospectra{1,t,1})./mean_ratiospectra{1,t,1};
plot(f,ratio_moving_stopped*100);








     
     
     
     
     
