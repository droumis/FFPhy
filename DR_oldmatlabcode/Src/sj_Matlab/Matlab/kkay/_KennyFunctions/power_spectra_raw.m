
%% first analysis: power spectrum of each tetrode, each epoch (one day)

% input is a struct containing tetrode fields
% eegstruct{something}{day}{tetrode}.<field>

params.Fs = 1500;
params.err = [2 0.05];
params.fpass = [0 400];
params.tapers = [3 5];

spectra=cell(size(eegstruct));            
winsize = params.Fs*2 ;             %% approach: calculate spectra in 1 s windows

for d=5
    for e=1:6                                            
        nowindows = floor(length(eegstruct{d}{e}{t}.data)/winsize)-1;
        for t=8:14
            for w=1:nowindows
                [spectra{d}{e}{t}(w,:),freqs,Serr] = ...
                    mtspectrumc(double(eegstruct{d}{e}{t}.data((1+(w-1)*winsize):(w+1)*winsize)),params);   %% only analyzing 1 minute / epoch
            end
        end
    end
end


 %%%%%%%% average all sleep epoch spectra

for d=5
    for t=8:14                
        sleepdummy{t}=[];
        for e=1:5           %% sleep
            sleepdummy{t} = [sleepdummy{t} ; spectra{d}{e}{t}];   %% only analyzing 1 minute / epoch
        end
    end
end

for d=5
    for t=8:14
        rundummy{t}=[];
        for e=6          %% run
            rundummy{t} = [rundummy{t} ; spectra{d}{e}{t}];   %% only analyzing 1 minute / epoch
        end
    end
end



%% take average

sleepmean=cell(1,14)
runmean=cell(1,14)

for t=8:14
    sleepmean{t}=mean(sleepdummy{t},1);
end

for t=8:14
    runmean{t}=mean(rundummy{t},1);
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








