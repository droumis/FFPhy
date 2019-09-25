function out = wavelettest(eegstruct,days,epochs,tetrodes,pos,scale,mother,windowsize,notests)

% Use this function to validate the wavelet spectrogram on some test data.
% z-scores entire epoch and plots some sample 10 s epochs of wavelet + eeg

% choose, say, a run and sleep epoch

% inputs:
% eegstruct = % put some seconds-minutes longs epochs of your EEG data
        % if an eegstruct -- contains ONLY the d,e,t that you want to analyze
% pos = 
% scale = choose a value between 0 and 1
% mother = choose either 'MORLET', 'DOG', or 'PAUL' -- experiment
% windowsize = duration of test window, in sec
% notests = number of random windows to look at

Fs=eegstruct{days(1)}{epochs(1)}{tetrodes(1)}.samprate;
windowsamples=windowsize*Fs;

%chronux mt params
params.tapers=[3 5];
params.Fs=Fs;
params.fpass=[0 400];
params.Serr=0;
movingwin=[100 25]/1000;

% wavelet calculation
 [wavezspectrogram,~,wavefreq] = waveletterz(eegstruct,Fs,scale,mother,days,epochs,tetrodes);
% dpss calculation
[mtzspectrogram,mttimes,mtfreq] = mtspecgramcz(eegstruct,movingwin,params,days,epochs,tetrodes);

% plot
for k=1:notests

    figure
    
    % pick random day, epoch, and start index + make times vector
    d=days(ceil(rand*length(days)));               
    e=epochs(ceil(rand*length(epochs)));
    N=length(eegstruct{d}{e}{tetrodes(1)}.data);
    startindex=int32(floor(rand*(N-windowsamples)));        % random start index within the epoch
        clock_time=double((startindex-1)/Fs+eegstruct{d}{e}{tetrodes(1)}.starttime);   % in sec, since open all files
    endindex=int32(startindex+windowsamples-1);
    times=((1:windowsamples+1)-1)/Fs;
    
    % wavelet spectrogram plot
    subplot(4,2,1)
    flip=wavezspectrogram{d}{e}{tetrodes}';
    imagesc(times,wavefreq,flip(startindex:endindex,:)',[-1,3]);
    colormap('jet')
    set(gca,'YDir','normal');
    
    % multitaper spectrogram plot
    subplot(4,2,3)   
    startwin=lookup(double((startindex-1)/Fs),mttimes{d}{e}{tetrodes});    % index of first mt window
    %endwin=lookup(double((endindex-1)/Fs),mttimes{d}{e}{tetrodes});        % index of last mt window
    endwin=startwin+windowsamples/(Fs*movingwin(2));
    imagesc(mttimes{d}{e}{tetrodes}(startwin:endwin)-mttimes{d}{e}{tetrodes}(startwin), ...
        mtfreq, ...
        mtzspectrogram{d}{e}{tetrodes}(startwin:endwin,:)', ...
        [-1,3]);
    colormap('jet')
    set(gca,'YDir','normal');
    
    % eeg plot
    subplot(4,2,5)
    plot(times,eegstruct{d}{e}{tetrodes}.data(startindex:endindex));
    ylim([-500 500]);
    
    % velocity plot
    subplot(4,2,7)
    pos_start=lookup(clock_time,pos{d}{e}.data(:,1));                   % index in pos vector
        pos_Fs=29.97003;
        pos_end=floor(pos_start+windowsize*pos_Fs);
    plot(0:1/pos_Fs:((pos_end-pos_start)/pos_Fs),pos{d}{e}.data(pos_start:pos_end,5),'k','LineWidth',3);
    ylim([0 30])
    xlim([0 windowsize])
    
    % position plot
    subplot(4,2,[2 4 6 8])
    hold on
    plot(pos{d}{e}.data(:,2),pos{d}{e}.data(:,3),'LineWidth',2,'Color',[0.8 0.8 0.8]);    % all paths in grey
    plot(pos{d}{e}.data(pos_start:pos_end,2),pos{d}{e}.data(pos_start:pos_end,3),'k','LineWidth',3);  % during snippet, black
    
    
end

end

