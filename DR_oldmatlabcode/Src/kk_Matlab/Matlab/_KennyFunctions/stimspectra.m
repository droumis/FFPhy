function [] = stimspectra(eegstruct, stimdio, windowsize, velepoch)

%% here eegstruct is my own structure, where first level is epoch (1-7)
%% same for pulsesstruct
% windowsize is given in ms

for e=5 %% specify epoch
    eeg = eegstruct{e}{5}{e}{11};           %% day is 5, tetrode is 11
    pulses = stimdio{e}{1}.pulsetimes(:,1);
    %% for some odd reason, eeg.samprate is not an integer:
    eeg.samprate = round(eeg.samprate);


%% note that the eeg argument here is a structure containing these three fields:
%% .samprate, .data, and .starttime


%% currently, eeg.data is sampled at 2 kHz, originally 30 kHz

downsampled = 15;                             %% downsampling factor that was imposed in dayprocess
Fs = round(eeg.samprate);                        %% original Fs
t = eeg.starttime:(1/Fs):(eeg.starttime+(length(eeg.data)-1)*(1/Fs));

pulses = pulses/10000;
windowsize = windowsize/1000;
numsamples = round(windowsize*Fs) + 1;
sampleoffset = (numsamples-1)/2;                      %% to accommodate matlab's vector indexing
eegwindows = nan([numsamples length(pulses)]);       %% initializes output vector

% chronux params
params.Fs = Fs;
params.err = [2 0.05];
params.fpass = [0 400];
params.tapers = [3 5];

%% all arguments are now in timestamps of seconds


%% copy eeg data into windows

for p=1:length(pulses)
    pulseindex = lookup(pulses(p),t);                                 %% looks up the eeg vector index of the specific pulse
    for s = -sampleoffset:1:sampleoffset                     %% iterates over all bins
       eegwindows(s+sampleoffset+1,p) = eeg.data(pulseindex+s);     %% copies data into window
    end
end


if 0
%% raw EEG plot
figure
seplot(-sampleoffset:1:sampleoffset,eegwindows','varType','std');
plot(-sampleoffset:1:sampleoffset,mean(eegwindows,2));
end



if 0
%% pre- vs. post- spectra plot  +  pulse-by-pulse power change

% pre-pulse spectra
for p=1:length(pulses)         
    [prespectra(:,p),freqs,Serr] = mtspectrumc(eegwindows(1:(sampleoffset+1),p),params);
end
    
% post-pulse spectra
for p=1:length(pulses)         
    [postspectra(:,p),freqs,Serr] = mtspectrumc(eegwindows((sampleoffset+1):numsamples,p),params);
end

% pre + post plot
figure
semilogy(freqs,mean(prespectra,2));
hold on;
semilogy(freqs,mean(postspectra,2),'r');

% z-score plot            %% keep in mind pcolor's color n'lzn range
figure
meanspectrum = mean(prespectra,1);
stdspectrum = std(prespectra,[],1);
diffpostspectra = bsxfun(@minus,postspectra,meanspectrum);
zscorepostspectra = bsxfun(@rdivide,diffpostspectra,stdspectrum);
%pcolor([1:length(pulses)],freqs,zscorepostspectra);

pcolor([1:length(pulses)],freqs(1:20),zscorepostspectra(1:20,:));      %% theta plot
figure
pcolor([1:length(pulses)],freqs(20:100),zscorepostspectra(20:100,:));      %% slow gamma plot
figure
pcolor([1:length(pulses)],freqs(75:200),zscorepostspectra(75:200,:));      %% fast gamma plot

%plot(freqs,mean(zscorepostspectra,2))

end


if 1
%% spectrogram
params.fpass = [0 350];            
movingwin = [500 20]/1000;

for p=1:length(pulses)            % stim window
      [spectrograms2(:,:,p),t2,freqs2,Serr2] = mtspecgramc(eegwindows(:,p),movingwin,params);
end
[fullspectrograms2,junkt,junkf,junkserr] = mtspecgramc(eeg.data,movingwin,params);   %% entire epoch run

stdspectrogram2 = std(fullspectrograms2,1);       %%  std calculated from entire epoch
meanfullspectrograms2 = mean(fullspectrograms2,1);      %% mean for entire epoch  
meanspectrogram2 = mean(spectrograms2,3);
diffspectrogram2 = bsxfun(@minus,meanspectrogram2,meanfullspectrograms2);
zscorespectrogram2 = bsxfun(@rdivide,diffspectrogram2,stdspectrogram2);
t2=(t2-windowsize/2)*1000;


    if 1  % plots individual pulse-triggered spectrograms and velocity during window
        for g=1:5        %% do first 5 pulses    
            diffspectrogram = bsxfun(@minus,spectrograms2(:,:,g),meanfullspectrograms2);
            zscorespectrogram = bsxfun(@rdivide,diffspectrogram,stdspectrogram2);
            figure
            subplot(10,1,[1:8])
            imagesc(t2,freqs2,zscorespectrogram');
            set(gca,'YDir','normal');
            subplot(10,1,[1:8])
            
            velocitywindow = nan(1,numsamples);    %% initialize
            velocityindex = lookup(pulses(p),velepoch(:,1));   %% retrieve velocity alignment
            for s = -sampleoffset:1:sampleoffset                     %% iterates over all bins
                velocitywindow(s+sampleoffset+1) = velepoch(pulseindex+s);     %% copies data into window
            end
            
            plot(-windowsize/2:1/Fs:windowsize/2,velocitywindow)
        end
    end
    

if 0   %% plot mean over all pulses in epoch
    figure
    imagesc(t2,freqs2,zscorespectrogram2');            %% plot of mean of all pulses
    set(gca,'YDir','normal')
end



    if 0    % full epoch spectrogram (need to narrow window)
        difffullspectrogram2 = bsxfun(@minus,fullspectrograms2,meanfullspectrograms2);
        zscorefullspectrogram2 = bsxfun(@rdivide,difffullspectrogram2,stdspectrogram2);    % z-score each time window
        figure
        imagesc(junkt(9000:16000),junkf(30:180),zscorefullspectrogram2(9000:16000,30:180)');   %% alter your choice of time window.. if too large imagesc can't resolve!
        set(gca,'YDir','normal');
    end    
    
    jay=1;
    


        
    
    if 0  %% for state-analysis, just sticking this in -- use debug to stop here w/ velocity matrix
        nostate1pulses = 0;             %% initialize count
        nostate2pulses = 0;
        nostate3pulses = 0;
        state1spectrograms = [];         %% initialize arrays
        state2spectrograms = [];
        state3spectrograms = [];   
        for p=1:length(pulses)
            velocityindex = lookup(pulses(p),velepoch(:,1));
            if (velepoch(velocityindex,2) < .1)           % state 1 
                nostate1pulses = nostate1pulses + 1;
                state1spectrograms = cat(3,state1spectrograms,spectrograms2(:,:,p));
            elseif (velepoch(velocityindex,2) > 8)       % state 2
                nostate2pulses = nostate2pulses + 1;
                state2spectrograms = cat(3,state2spectrograms,spectrograms2(:,:,p));
            else                                         % state 3
                nostate3pulses = nostate3pulses + 1;
                state3spectrograms = cat(3,state3spectrograms,spectrograms2(:,:,p));
            end
        end
            nostate1pulses
            nostate2pulses
            nostate3pulses
            
            
                if 0   %% tests your spectrogram for harmonic noise
                    fundamental_f = 100;
                    testsignal = 4*sin(2*pi*(fundamental_f)*t);
                    [testspec,t3,f3,err3] = mtspecgramc(testsignal,movingwin,params);
                    testspec = bsxfun(@minus,testspec,meanfullspectrograms2);
                    testspec = bsxfun(@rdivide,testspec,stdspectrogram2);
                    figure
                    imagesc(t3,f3,testspec');            %% plot of mean of all pulses
                    set(gca,'YDir','normal')
                end
    
            %% take average, then z-score
         state1 = mean(state1spectrograms,3);
         state2 = mean(state2spectrograms,3);
         state3 = mean(state3spectrograms,3);
         
         state1 = bsxfun(@minus,state1,meanfullspectrograms2);
         state1 = bsxfun(@rdivide,state1,stdspectrogram2);
         
         state2 = bsxfun(@minus,state2,meanfullspectrograms2);
         state2 = bsxfun(@rdivide,state2,stdspectrogram2);

         state3 = bsxfun(@minus,state3,meanfullspectrograms2);
         state3 = bsxfun(@rdivide,state3,stdspectrogram2);
         
         figure
         imagesc(t2,freqs2,state1');            %% state 1
         set(gca,'YDir','normal')

         figure
         imagesc(t2,freqs2,state2');            %% state 2
         set(gca,'YDir','normal')
         
         figure
         imagesc(t2,freqs2,state3');            %% state 3
         set(gca,'YDir','normal')
         
         figure                                  %% state 1 and 2
         imagesc(t2*2,freqs2,cat(2,state1',state2'));
         set(gca,'YDir','normal')
         
         figure                                  %% state 1 and 2 and3
         imagesc(t2*3,freqs2,cat(2,state1',state3',state2'));
         set(gca,'YDir','normal')
    end
    
    
    


end

if 0
%% non-scored spectrograms
params.fpass = [0 20];            %% theta

figure
for p=1:length(pulses)   
      [spectrograms3(:,:,p),t3,freqs3,Serr3] = mtspecgramc(eegwindows(:,p),movingwin,params);
end
meanspectrogram3 = mean(spectrograms3,3);
t3=(t3-windowsize/2)*1000;
pcolor(t3,freqs3,meanspectrogram3');
end


if 0
params.fpass = [50 200];            %% high

figure
for p=1:length(pulses)   
      [spectrograms4(:,:,p),t4,freqs4,Serr4] = mtspecgramc(eegwindows(:,p),movingwin,params);
end
meanspectrogram4 = mean(spectrograms4,3);
t4=(t4-windowsize/2)*1000;
pcolor(t4,freqs4,meanspectrogram4');

end










if 0
%%% chronux fft, averaged non-overlapping windows
    startsamp = 0;
    nowindows = 500;sanfran = 1;
    % test run to extract f
    for w=1:nowindows      
        onewindow = eeg.data(startsamp+(w-1)*round(windowsize*Fs) + (1:round(windowsize*Fs)));
        [spectra(w,:),f,Serr] = mtspectrumc(onewindow,params);                                     
    end
    meanspectrum = mean(spectra,1);
    diffspectra = bsxfun(@minus,spectra,meanspectrum);          %% spectra - meanspectrum -- absolute subtraction
    stdspectrum = std(diffspectra,[],1);
    zscorespectra = bsxfun(@rdivide,diffspectra,stdspectrum);          %% z-score

    figure
    
    
    subplot(10,1,[1:3])
    plot((startsamp + (1:nowindows*windowsize/Fs))*(1/Fs), eeg.data(startsamp + (1:nowindows*windowsize*Fs)))
    subplot(10,1,[4:10])
    pcolor([1:nowindows],f,zscorespectra');
    shading flat
    caxis([0 4])
    
    divspectra = bsxfun(@rdivide,spectra,meanspectrum);
    
    figure
    subplot(10,1,[1:3])
    plot((startsamp + (1:nowindows*windowsize*Fs))*(1/Fs),eeg.data(startsamp + (1:nowindows*windowsize*Fs)))
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





















