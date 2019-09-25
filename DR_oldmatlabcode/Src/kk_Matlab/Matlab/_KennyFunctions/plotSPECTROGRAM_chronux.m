
%% Plots continuous spectrogram (z-scored to entire day).
   % ingredients
        % 1. eegstruct
        % 2. pos (to overlay velocity on spectrogram)

%% Set these parameters manually.

days = 4;
epochs = 4;
duration = 180;  % in seconds
if 1
eegstruct = loadeegstruct('/data12/mari/Egy/','egy','eeg',days);
pos = loaddatastruct('/data12/mari/Egy/','egy','pos',days);
end
tetrodes = [2 17];
Fs = 1500;                                      

% chronux params
movingwin = [1500 10]/1000;
params.Fs = Fs;
params.err = [2 0.05];
params.fpass = [0 20];
params.tapers = [3 5];

tetno=tetrodes(end);
epno=epochs(end);


%% First, calculate the continuous spectrogram for all epochs for the day.
%% At the same time, concatenate the spectrograms from each epoch and
%% calculate the mean and std.

spectrograms={};
meandayspectrum={};
stddayspectrum={};

for d=days
    for t=tetrodes                                
        dummy=[];
        for e=epochs
            [S,timevec,freqs,junkserr] = mtspecgramc(eegstruct{days}{e}{t}.data(1:(duration*Fs)),movingwin,params);
            spectrograms{d}{e}{t}=S;
            dummy=[dummy;S];
        end
    meandayspectrum{d}{t}=mean(dummy,1);
    stddayspectrum{d}{t}=std(dummy,1);
    end
end

% Second, z-score the spectrograms.

for d=days
for e=epochs
    for t=tetrodes
        spectrograms{d}{e}{t}=bsxfun(@minus,spectrograms{d}{e}{t},meandayspectrum{d}{t});   
        spectrograms{d}{e}{t}=bsxfun(@rdivide,spectrograms{d}{e}{t},stddayspectrum{d}{t});
    end
end
end


    
% Third, plot z-score spectra, all

for d = days
    for e=epochs
for t=2
        figure
        
        % pos
        subplot(6,1,1)
        postimes = pos{d}{e}.data(:,1)-pos{d}{e}.data(1,1);
        lastindex = lookup(duration,postimes);
        plot(postimes(1:lastindex),pos{d}{e}.data(1:lastindex,5));
        
        % spectrogram
        subplot(6,1,2:6)
        imagesc(timevec,freqs,mean(spectrograms{d}{e}{t},3)',[0,3]);
            set(gca,'YDir','normal');
        string = sprintf('%d',t);
        title(string);
end
    end
end







        
        
        
        
        
        
        
        
        
        
        
