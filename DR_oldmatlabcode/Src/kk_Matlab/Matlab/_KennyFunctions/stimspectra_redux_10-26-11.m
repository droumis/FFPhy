
%% here eegstruct is my own ad hoc struct, where first level is epoch (1-7)
%% same for pulsesstruct
% windowsize is given in ms

%%%%%%% COPY EEG DATA INTO EEGWINDOWS

    eegwindows=cell(1,7);
    
for e=1:4 %% specify epoch
    eeg = eegstruct{e}{5}{e}{11};           %% day is 5, tetrode is 11        %% used loadeegstruct
    pulses = stimdio{e}{1}.pulsetimes(2:(end-1),1);
    %% for some odd reason, eeg.samprate is not an integer:
    eeg.samprate = round(eeg.samprate);

%% note that the eeg argument here is a structure containing these three fields:
%% .samprate, .data, and .starttime

Fs = round(eeg.samprate);          %% currently, eeg.data is sampled at 2 kHz, originally 30 kHz             
t = eeg.starttime:(1/Fs):(eeg.starttime+(length(eeg.data)-1)*(1/Fs));

pulses = pulses/10000;
windowsize = 800/1000;
numsamples = round(windowsize*Fs) + 1;
sampleoffset = (numsamples-1)/2;                      %% to accommodate matlab's vector indexing
eegwindows{e} = nan(numsamples,length(pulses));       %% initializes output vector

% chronux params
params.Fs = Fs;
params.err = [2 0.05];
params.fpass = [0 350];
params.tapers = [3 5];

%% all arguments are now in timestamps of seconds

for p=1:length(pulses)
    pulseindex = lookup(pulses(p),t);                                 %% looks up the eeg vector index of the specific pulse
    for s = -sampleoffset:1:sampleoffset                     %% iterates over all bins
       eegwindows{e}(s+sampleoffset+1,p) = eeg.data(pulseindex+s);     %% copies data into window
    end
end
end

%%%%%%%%%%%% CONCATENATE EEGWINDOWS

alleegwindows=[];
for e=[1 2 3 4]
    alleegwindows=cat(2,alleegwindows,eegwindows{e});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% raw EEG plot for each epoch

for e=1:7
figure
%seplot(-sampleoffset:1:sampleoffset,(eegwindows{e})','varType','std');
plot(-sampleoffset:1:sampleoffset,mean(eegwindows{e},2));
string1=sprintf('epoch %d,',e);
string2=sprintf(' %d pulses',size(eegwindows{e},2));
title([string1 string2]);
end

%% raw EEG plot for lumped epochs
figure
hold on
for e=[1 2 3 4 6 7]
%seplot(-sampleoffset:1:sampleoffset,(eegwindows{e})','varType','std');
plot(-sampleoffset:1:sampleoffset,mean(alleegwindows,2));
end



%%%%%%%%%%%%%%% pre- vs. post- spectra, EPOCH BY EPOCH

figure
hold on
for e=[1 2 3 4]
% pre-pulse spectra
for p=1:size(eegwindows{e},2)         
    [prespectra(:,p),freqs,Serr] = mtspectrumc(eegwindows{e}(1:(sampleoffset+1),p),params);
end
    
% post-pulse spectra
for p=1:size(eegwindows{e},2)          
    [postspectra(:,p),freqs,Serr] = mtspectrumc(eegwindows{e}((sampleoffset+1):numsamples,p),params);
end
    if mod(e,2)==0
        semilogy(freqs,mean(postspectra,2)./mean(prespectra,2),'g','LineWidth',2);
    else
        semilogy(freqs,mean(postspectra,2)./mean(prespectra,2),'b','LineWidth',2);
    end
end
set(gca,'xscale','log')

%% pre- vs. post- spectra, ALLEPOCHS

% pre-pulse spectra
for p=1:size(alleegwindows,2)         
    [prespectra(:,p),freqs,Serr] = mtspectrumc(alleegwindows(1:(sampleoffset+1),p),params);
end
    
% post-pulse spectra
for p=1:1:size(alleegwindows,2)         
    [postspectra(:,p),freqs,Serr] = mtspectrumc(alleegwindows((sampleoffset+1):numsamples,p),params);
end

semilogy(freqs,mean(postspectra,2)./mean(prespectra,2),'r','LineWidth',4);
set(gca,'xscale','log')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% mean z-score spectra plot           
figure
meanspectrum = mean(prespectra,2);
stdspectrum = std(prespectra,0,2);
diffpostspectra = bsxfun(@minus,postspectra,meanspectrum);
zscorepostspectra = bsxfun(@rdivide,diffpostspectra,stdspectrum);

plot(freqs,mean(zscorepostspectra,2))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555



if 1
%% spectrogram

params.fpass = [0 400];            
movingwin = [100 10]/1000;

%%%%%%%%%%% DAY MEAN SPECTROGRAM
dummy=[];
        for e=[1 2 3 4]
            [S_epoch,junkt,junkf,junkserr] = mtspecgramc(eegstruct{e}{5}{e}{11}.data,movingwin,params);   %% day 5, tetrode 11
            dummy=[dummy;S_epoch];        
        end
meandayspectra=mean(dummy,1);
stddayspectra=std(dummy,0,1);

%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%% CALCULATE RAW SPECTROGRAMS

spectrograms=cell(1,7);

[dummy,times,freqs,Serr] = mtspecgramc(eegwindows{1}(:,1),movingwin,params);  %% for obtaining dimensions

for e=[1 2 3 4]
    spectrograms{e}=nan(size(dummy,1),size(dummy,2),size(eegwindows{e},2));   %% initialize
    for p=1:size(eegwindows{e},2)
      [spectrograms{e}(:,:,p),times,freqs,Serr] = mtspecgramc(eegwindows{e}(:,p),movingwin,params);
    end
end

%%%%%%%%%%% Z-SCORE SPECTRA

zscorespectrograms=spectrograms;

for e=[1 2 3 4]
    for p=1:size(eegwindows{e},2)
        zscorespectrograms{e}(:,:,p) = bsxfun(@minus,spectrograms{e}(:,:,p),meandayspectra);
        zscorespectrograms{e}(:,:,p) = bsxfun(@rdivide,zscorespectrograms{e}(:,:,p),stddayspectra);
    end
end

%plot (subplot)
figure
for e=[1 2 3 4]
    subplot(1,4,e)
    imagesc(times,freqs,mean(zscorespectrograms{e},3)',[0 2]);
    set(gca,'YDir','normal');
end

%plot (individual)

for e=[1 2 3 4]
    figure
    imagesc(times,freqs,mean(zscorespectrograms{e},3)',[0 1]);
    set(gca,'YDir','normal');
    string=sprintf('epoch %d',e);
    title(string);
end








%     if 1  % plots individual pulse-triggered spectrograms and velocity during window
%         for g=1:5        %% do first 5 pulses    
%             diffspectrogram = bsxfun(@minus,spectrograms2(:,:,g),meanfullspectrograms2);
%             zscorespectrogram = bsxfun(@rdivide,diffspectrogram,stdspectrogram2);
%             figure
%             subplot(10,1,[1:8])
%             imagesc(t2,freqs2,zscorespectrogram');
%             set(gca,'YDir','normal');
%             subplot(10,1,[1:8])
%             
%             velocitywindow = nan(1,numsamples);    %% initialize
%             velocityindex = lookup(pulses(p),velepoch(:,1));   %% retrieve velocity alignment
%             for s = -sampleoffset:1:sampleoffset                     %% iterates over all bins
%                 velocitywindow(s+sampleoffset+1) = velepoch(pulseindex+s);     %% copies data into window
%             end
%             
%             plot(-windowsize/2:1/Fs:windowsize/2,velocitywindow)
%         end
%     end
%     
% 
% if 0   %% plot mean over all pulses in epoch
%     figure
%     imagesc(t2,freqs2,zscorespectrogram2');            %% plot of mean of all pulses
%     set(gca,'YDir','normal')
% end
% 
% 
% 
%     if 0    % full epoch spectrogram (need to narrow window)
%         difffullspectrogram2 = bsxfun(@minus,fullspectrograms2,meanfullspectrograms2);
%         zscorefullspectrogram2 = bsxfun(@rdivide,difffullspectrogram2,stdspectrogram2);    % z-score each time window
%         figure
%         imagesc(junkt(9000:16000),junkf(30:180),zscorefullspectrogram2(9000:16000,30:180)');   %% alter your choice of time window.. if too large imagesc can't resolve!
%         set(gca,'YDir','normal');
%     end    
%     
%     jay=1;
%     


        %%%%%%%%%%%%%%%%%%%%%%%%% VELOCITY CLASSIFICATION OF STATES,
        %%%%%%%%%%%%%%%%%%%%%%%%% SPECTROGRAMS, SPECTRA
    
    if 0  %% for state-analysis, just sticking this in -- use debug to stop here w/ velocity matrix
        nostate1pulses = 0;             %% initialize count
        nostate2pulses = 0;
        nostate3pulses = 0;
        state1spectrograms=[];
        state2spectrograms=[];
        state3spectrograms=[];
        prespectra_state=cell(1,3);   %% state 1, 2, 3
        postspectra_state=cell(1,3);
        camera_Fs=60/(1.001*2);

 for e=[2 4]
    if e==1 || e==3 || e==7
        col=11;
    else
        col=9;
    end
    pulses = stimdio{e}{1}.pulsetimes(2:(end-1),1);
    for p=1:length(pulses)
            velocityindex = lookup(pulses(p),pos{5}{e}.data(:,1)*10000);
            if (pos{5}{e}.data(velocityindex,col) <2)       % state 1
                nostate1pulses = nostate1pulses + 1;
                state1spectrograms = cat(3,state1spectrograms,zscorespectrograms{e}(:,:,p));
                [prespectra_state{1}(:,p),freqs,Serr] = mtspectrumc(eegwindows{e}(1:(sampleoffset+1),p),params);
                [postspectra_state{1}(:,p),freqs,Serr] = mtspectrumc(eegwindows{e}((sampleoffset+1):numsamples,p),params);
            elseif (pos{5}{e}.data(velocityindex,col) > 8)       % state 2
                nostate2pulses = nostate2pulses + 1;
                state2spectrograms = cat(3,state2spectrograms,zscorespectrograms{e}(:,:,p));
                [prespectra_state{2}(:,p),freqs,Serr] = mtspectrumc(eegwindows{e}(1:(sampleoffset+1),p),params);
                [postspectra_state{2}(:,p),freqs,Serr] = mtspectrumc(eegwindows{e}((sampleoffset+1):numsamples,p),params);
            else                                         % state 3
                nostate3pulses = nostate3pulses + 1;
                state3spectrograms = cat(3,state3spectrograms,zscorespectrograms{e}(:,:,p));
                [prespectra_state{3}(:,p),freqs,Serr] = mtspectrumc(eegwindows{e}(1:(sampleoffset+1),p),params);
                [postspectra_state{3}(:,p),freqs,Serr] = mtspectrumc(eegwindows{e}((sampleoffset+1):numsamples,p),params);
 
            end
    end
 end

             nostate1pulses
            nostate2pulses
            nostate3pulses
         
            
            figure
            hold on
        semilogy(freqs,mean(postspectra_state{1},2)./mean(prespectra_state{1},2),'b','LineWidth',2);
        semilogy(freqs,mean(postspectra_state{2},2)./mean(prespectra_state{2},2),'r','LineWidth',2);
        semilogy(freqs,mean(postspectra_state{3},2)./mean(prespectra_state{3},2),'g','LineWidth',2);
set(gca,'xscale','log')
    
    
            
 
           

    figure  
    imagesc(times,freqs,mean(state1spectrograms,3)',[-0.2 3])
    set(gca,'YDir','normal');
    s=sprintf('state1, %d pulses',nostate1pulses);
    title(s);
    colorbar
    figure
    imagesc(times,freqs,mean(state2spectrograms,3)',[-0.2 3])
    set(gca,'YDir','normal');
    s=sprintf('state2, %d pulses',nostate2pulses);
    title(s);
    colorbar
    
    
    
    figure
    imagesc(times,freqs,mean(state3spectrograms,3)',[-0.5 3.2])
    set(gca,'YDir','normal');
    s=sprintf('state3, %d pulses',nostate3pulses);
    title(s);
    
    
    %% individual stim plots
    
    for p=5:10
    figure
    imagesc(times,freqs,state1spectrograms(:,:,p)',[0 6])
    set(gca,'YDir','normal');
    s=sprintf('state1, %d pulses',nostate1pulses);
    title(s);
    end
    
    for p=10:20
    figure
    imagesc(times,freqs,state2spectrograms(:,:,p)',[0 6])
    set(gca,'YDir','normal');
    s=sprintf('state2, %d pulses',nostate2pulses);
    title(s);
    end
    
    
    %%%%%%%% subplot
        figure
    subplot(3,1,1)
    imagesc(times,freqs,mean(state0spectrograms,3)',[0 0.8])
    set(gca,'YDir','normal');
    s=sprintf('state0, %d pulses',nostate0pulses);
    title(s);
    subplot(3,1,2)
    imagesc(times,freqs,mean(state1spectrograms,3)',[0 0.8])
    set(gca,'YDir','normal');
    s=sprintf('state1, %d pulses',nostate1pulses);
    title(s);
   subplot(3,1,3)
    imagesc(times,freqs,mean(state2spectrograms(:,:,[1:13 15:end]),3)',[0 0.8])      %% 14 is cabled
    set(gca,'YDir','normal');
    s=sprintf('state2, %d pulses',nostate2pulses);
    title(s);
    
    

    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
            
            
            
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





















