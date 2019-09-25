function out = plot_eeg_position(eegstruct,days,epochs,tetrodes,pos,windowsize,notests)

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

Fs=1000;
windowsamples=windowsize*Fs;


% plot
for k=1:notests

    figure
    
    % pick random day, epoch, and start index + make times vector
    d=days(ceil(rand*length(days)));               
    e=epochs(ceil(rand*length(epochs)));
    N=length(eegstruct{d}{e}{tetrodes(1)}.data);
    startindex=int32(floor(rand*(N-windowsamples)));        % random start index within the epoch
    
        clock_time=double((startindex-1)/Fs);   % in sec, since open all files
    endindex=int32(startindex+windowsamples-1);
    times=((1:windowsamples)-1)/Fs;
        
    % eeg plot
    subplot(4,2,[1 3 5])
    yspacing=750;
    hold on
    for t=16:31
        plot(times,eegstruct{d}{e}{t}.data(startindex:endindex)-yspacing*(t-16));
    end
    axis tight;

    
    % velocity plot
    subplot(4,2,7)
    pos_start=lookup(clock_time,(pos{e}.data(:,1)-pos{e}.data(1,1)));                   % index in pos vector
    pos_end=lookup(clock_time+windowsize,(pos{e}.data(:,1)-pos{e}.data(1,1)));
        fps=30;
    plot(pos{e}.data(pos_start:pos_end,5),'k','LineWidth',3);

        axis tight
            ylim([0 5])
    
    % position plot
    subplot(4,2,[2 4 6 8])
    hold on
    plot(pos{e}.data(:,2),pos{e}.data(:,3),'LineWidth',2,'Color',[0.8 0.8 0.8]);    % all paths in grey
    plot(pos{e}.data(pos_start:pos_end,2),pos{e}.data(pos_start:pos_end,3),'k','LineWidth',3);  % during snippet, black
    ylim([15 50])
    xlim([20 80])
    
    
end

end

