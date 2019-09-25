
%% 6.14.12  -- calculates ratio of spectra from one epoch from another

animalname = 'Bashir';                              % to label graphs later
windowsize=0.5;           % in seconds
possamplefreq=29.97;    % ETSC video standard
    windowsize_samp=floor(2*29.97);
days = 4;                                           % days to analyze
tetrodes = 8:14;                                    % tetrodes to analyze
epochs = 2:3;                                       % epochs to analyze

Fs=eegstruct{days(end)}{epochs(end)}{tetrodes(end)}.samprate;   % probably just 1500
    eegsampwin=floor(windowsize*Fs);

params.Fs = 1500;
params.err = [2 0.05];
params.fpass = [0 400];
params.tapers = [3 5];

spectra=cell(size(eegstruct));            
winsize = params.Fs*windowsize ;             %% approach: calculate spectra in 1 s windows

for d=days
    for e=epochs                                            
        nowindows = floor(length(eegstruct{d}{e}{t}.data)/winsize)-1;
        for t=tetrodes
            for w=1:nowindows
                [spectra{d}{e}{t}(w,:),freqs,Serr] = ...
                    mtspectrumc(double(eegstruct{d}{e}{t}.data((1+(w-1)*winsize):(w+1)*winsize)),params);   %% only analyzing 1 minute / epoch
            end
        end
    end
end


%% take average

sleep_epochs = [3];
run_epochs = [2];
sleepmean=cell(tetrodes(end));
runmean=cell(tetrodes(end));

for t=tetrodes
    dummy=[];
    for e=sleep_epochs
        dummy=[dummy; spectra{d}{e}{t}];
    end
    sleepmean{t}=mean(dummy,1);
end

for t=tetrodes
    dummy=[];
    for e=run_epochs
        dummy=[dummy; spectra{d}{e}{t}];
    end
    runmean{t}=mean(dummy,1);
end




%% ratio of run to sleep spectra 

figure
hold on

for t=8:14
    if t==9 
        semilogx(freqs,runmean{t}./sleepmean{t},'g','LineWidth',2)
    elseif t==13 || t==8
        semilogx(freqs,runmean{t}./sleepmean{t},'LineWidth',2,'Color',[0 0.5 0])
    elseif t==12
        semilogx(freqs,runmean{t}./sleepmean{t},'r','LineWidth',2)
    elseif t==10
        semilogx(freqs,runmean{t}./sleepmean{t},'m','LineWidth',2)
    else
        semilogx(freqs,runmean{t}./sleepmean{t},'k','LineWidth',2)        
    end
end
ylim([-0.5 2])








