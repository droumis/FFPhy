
% load eeg for run sessions on two tetrodes
    load('/data13/mcarr/Fiv/EEG/Fiveeg06-2-11.mat')
    mec6211 = eeg{1,6}{1,2}{1,11}.data;
    load('/data13/mcarr/Fiv/EEG/Fiveeg06-2-16.mat')
    ca16216 = eeg{1,6}{1,2}{1,16}.data;

% set sampling frequency and frequency band
    params = {};
    params.Fs = 1500;
    params.fpass = [40 140];

% generate spectograms for each
    [Smec,tmec,fmec]=mtspecgramc(mec6211,[1 0.5],params);
    [Sca1,tca1,fca1]=mtspecgramc(ca16216,[1 0.5],params);

% plot spectograms
    figure;
    subplot(2,1,1); plot_matrix(Smec,tmec,fmec, 'l');
    subplot(2,1,2); plot_matrix(Sca1,tca1,fca1, 'l');
    
    %____________________________________________________________
    
    % compute coherency between two tetrodes
    
    
    % load eeg for run sessions on two tetrodes
    load('/data13/mcarr/Fiv/EEG/Fiveeg06-2-07.mat')
    mec627 = eeg{1,6}{1,2}{1,7}.data(1:30000);    
    load('/data13/mcarr/Fiv/EEG/Fiveeg06-2-16.mat')
    ca16216 = eeg{1,6}{1,2}{1,16}.data(1:30000);

    % set sampling frequency and frequency band
    params = [];
    params.Fs = 1500;
    params.fpass = [0 250];
    params.err = [1 0.05];
    
    % compute coherency
    [C,phi,S12,S1,S2,f,confC,phistd]=coherencyc(mec627,ca16216,params);
    
    % look at it backwards (to decouple in time, but not in frequency)
    
    n = length(mec627);
    mecback = [];
    for i = 1:n
        mecback(i,1) = mec627(n+1-i);
    end
    
    % compute coherency
    [Cback,phiback,S12back,S1back,S2back,fback,confCback,phistdback]=coherencyc(mecback,ca16216,params);
    
    % plot both
    figure
    plot(f, C)
    hold on
    plot(fback, Cback, 'r')