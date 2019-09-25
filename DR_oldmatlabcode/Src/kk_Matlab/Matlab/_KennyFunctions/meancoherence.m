function out = meancoherence(eegdata1, eegdata2, freqrange)

% Calculates mean coherence between two pairs of eegdata (one epoch) for freq. in range [lowfreq highfreq]
    % eegdata1 and 2 are in the eeg data structure of our lab

    % uses chronux's coherencyc 
    
    params.tapers = [3 5];
    params.Fs = 1500;
    params.err = 0;
    params.trialave = 1;

    windowsize = .5; % default is 1 sec windows
        windowsize_samp = windowsize * params.Fs;
    
    if round(eegdata1.samprate) ~= round(eegdata2.samprate)
        disp('difference in samprate..')
        keyboard
    end
        
    eegtimes1 = geteegtimes(eegdata1);
    eegtimes2 = geteegtimes(eegdata2); 
        % truncate ends
        starttime = eegtimes1(1) + 1;   
        endtime = eegtimes1(end) - 1;
        
        startind1 = lookup(starttime,eegtimes1);
        startind2 = lookup(starttime,eegtimes2);
        endind1 = lookup(endtime,eegtimes1);
            datalength = endind1 - startind1 + 1;        
                numwindows = floor(datalength/windowsize_samp);
    
    % split up data into windows ("trials")
    
    data1 = [];
    data2 = [];
    
    for w = 1:numwindows
        data1(w,:) = double( eegdata1.data((startind1 + (w-1)*windowsize_samp):...
                               (startind1 + w * windowsize_samp - 1) , 1));
        data2(w,:) = double( eegdata2.data((startind2 + (w-1)*windowsize_samp):...
                               (startind2 + w * windowsize_samp - 1) , 1 ));
    end
            
    % calculate coherogram
    [C,~,~,~,~,freqs]=coherencyc(data1,data2,params);
    
    lowindex = lookup(freqrange(1),freqs);
    highindex = lookup(freqrange(2),freqs);
    
    out = mean(C(lowindex:highindex));

    disp(sprintf('coherence: %d',out));


end