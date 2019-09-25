function [] = stimspectra(eegstruct, stimdio, windowsize)

%% here eegstruct is my own structure, where first level is epoch (1-7)
%% same for pulsesstruct
% windowsize is given in ms

for e=2
    eeg = eegstruct{e}{5}{e}{11};           %% day is 5, tetrode is 11
    pulses = stimdio{e}{1}.pulsetimes(:,1);
    %% for some odd reason, eeg.samprate is not an integer:
    eeg.samprate = round(eeg.samprate);


%% note that the eeg argument here is a structure containing these three fields:
%% .samprate, .data, and .starttime


%% currently, eeg.data is sampled at 2 kHz, originally 30 kHz

downsampled = 15;                             %% downsampling factor that was imposed in dayprocess
Fs = round(eeg.samprate);                        %% original Fs
t = eeg.starttime:Fs:(eeg.starttime+(length(eeg.data)-1)*Fs);

pulses = pulses/10000;
windowsize = windowsize/1000;
numsamples = round(windowsize/Fs) + 1;
sampleoffset = (numsamples-1)/2;                      %% to accommodate matlab's vector indexing
eegwindows = nan([numsamples length(pulses)]);       %% initializes output vector

%% all arguments are now in timestamps of seconds


%% copy eeg data into windows

for p=1:length(pulses)
    pulseindex = lookup(pulses(p),t);                                 %% looks up the eeg vector index of the specific pulse
    for s = -sampleoffset:1:sampleoffset                     %% iterates over all bins
        eegwindows(s+sampleoffset+1,p) = eeg.data(pulseindex+s);     %% copies data into window
    end
end

%raw EEG plot
%seplot(-sampleoffset:1:sampleoffset,eegwindows,'varType','std');



% chronux params setting
params.Fs = 1/Fs;
params.err = [2 0.05];
params.fpass = [0 250];
params.tapers = [3 5];

% pre- vs. post- spectra plot
%hold on;
%figure;
%title('pre vs. post- spectra, given epoch');
%[spectrum,freqs,Serr] = mtspectrumc(eegwindows(1:10000),params);
%plot(1:300,spectrum(1:300));



if 1
%%% chronux fft, averaged non-overlapping windows
    startsamp = 0;
    nowindows = 500;
    % test run to extract f
    for w=1:nowindows
        onewindow = eeg.data(startsamp+(w-1)*round(windowsize/Fs) + (1:round(windowsize/Fs)));
        [spectra(w,:),f,Serr] = mtspectrumc(onewindow,params);     %% Hanning window
    end
    meanspectrum = mean(spectra,1);
    stdspectrum = std(spectra,[],1);
    meansubspectra = bsxfun(@minus,spectra,meanspectrum);
    zscorespectra = bsxfun(@rdivide,spectra,stdspectrum);
    figure
    subplot(10,1,[1:3])
    plot((startsamp + (1:nowindows*windowsize/Fs))*Fs,eeg.data(startsamp + (1:nowindows*windowsize/Fs)))
    subplot(10,1,[4:10])
    pcolor([1:nowindows],f,zscorespectra');
    shading flat
    caxis([0 4])
    divspectra = bsxfun(@rdivide,spectra,meanspectrum);
    figure
    subplot(10,1,[1:3])
    plot((startsamp + (1:nowindows*windowsize/Fs))*Fs,eeg.data(startsamp + (1:nowindows*windowsize/Fs)))
    subplot(10,1,[4:10])
    pcolor([1:nowindows],f,divspectra');
    shading flat
    caxis([0 5])
    %set(gca,'xlim',[0 300])
    
end


if 0
%%% matlab fft, averaged non-overlapping windows
    hold on;
    title('averaged spectrum, each of 7 epochs');
    startsamp = eeg.samprate*10;                   %% time of first window
    nowindows = 10;
    eegwindows = [];      %% initialize
    for w=1:nowindows
        eegwindows(w,:) = (hanning(windowsize/Fs+1).*             ...       
                             (eeg.data((startsamp+w*windowsize/Fs): ...
                             (startsamp+(w+1)*windowsize/Fs))))';
    end
    spectra = fftn(eegwindows);                     %% Hanning window applied

end
    




end

















