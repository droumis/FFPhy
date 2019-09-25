
%% Plots continuous spectrogram (z-scored to entire day).
   % ingredients
        % 1. eegstruct
        % 2. pos (to overlay velocity on spectrogram)

%% Set these parameters manually.

days = 4;
epochs = 2:3;
epoch_to_plot = 2;
tetrodes = 1:14;
Fs = 1500;                                      

% chronux params
movingwin = [100 10]/1000;
params.Fs = Fs;
params.err = [2 0.05];
params.fpass = [0 400];
params.tapers = [3 5];

tetno=tetrodes(end);
epno=epochs(end);


%% First, calculate the continuous spectrogram for all epochs for the day.
%% At the same time, concatenate the spectrograms from each epoch and
%% calculate the mean and std.

spectrograms=cell(1,epno,tetno);
meandayspectrum=cell(1,tetno);
stddayspectrum=cell(1,tetno);

for d=days
    for t=tetrodes                                
        dummy=[];
        for e=epochs
            [S,junkt,junkf,junkserr] = mtspecgramc(eegstruct{days}{e}{t}.data,movingwin,params);
            spectrograms{1,e,t}=S;
            dummy=[dummy;S];
        end
    meandayspectrum{t}=mean(dummy,1);
    stddayspectrum{t}=std(dummy,1);
    end
end

% Second, z-score the spectrograms.

for e=epochs
    for t=tetrodes
        spectrograms{1,e,t}=bsxfun(@minus,spectrograms{1,e,t},meandayspectrum{t});   
        spectrograms{1,e,t}=bsxfun(@rdivide,spectrograms{1,e,t},stddayspectrum{t});
    end
end


    
% Third, plot z-score spectra, all
    figure
for t=8:14
      subplot(1,14,t)
        imagesc(times,freqs,mean(zscorespectra_lumpmean2{t},3)',[-0.2,5]);
            set(gca,'YDir','normal');
        string = sprintf('%d',t);
        title(string);
end
colorbar







        
        
        
        
        
        
        
        
        
        
        
        
    
    % plot some random individual ripple-triggered z-scored spectrograms
    % from an epoch
    
    e=2;  % set epoch manually
    for i=1:5
    figure
    N=size(zscorespectra_epochs{1,e,t},3);        % # events to choose from
    eventno=ceil(N*rand);
    string=sprintf('individual event no %d',eventno);
    title(string);
    for t=8:14
        subplot(1,14,t)
        imagesc(times,freqs,zscorespectra_epochs{1,e,t}(:,:,eventno)');
        set(gca,'YDir','normal');
        string = sprintf('%d',t);
        title(string);
        end
    end



