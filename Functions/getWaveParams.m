
function wavepms = getWaveParams(waveSet)
wavepms = setParams(waveSet);
end

%add custom wavelet param sets as cases below

function out = setParams(waveSet)
allNTDataCat = [];
switch waveSet
    case 'default'
        srate = 1500;
        win = [-1.5 1.5]; %in seconds. The morlet wavelet needs to be centered at time zero.. so always make this window symmetric about zero
        timeWin = win(1):1/srate:win(2);
        plotwin = [-.5 .5];
        basewin = plotwin(2); %baseline window to use for normalization
        plottimeWin = plotwin(1):1/srate:plotwin(2);
        %baseline indices. length of 1/2 plot window, centered on first plot ind
        % i.e. 0.5 sec basewindow starting at .25 sec before the start of the
        % plotwindow
        baseind = ([(win(2)-plotwin(2))-basewin/2 (win(2)-plotwin(2))+basewin/2]).*srate;
        % baseind = [1 abs(win(1))*srate-(abs(plotwin(1))*srate)]; %the entire duration before the plotwin
        % basedur = abs(win(1))+abs(win(2));
        % baselineind = ([1 (srate*basedur)+1]);
        pval = 0.05;
        zval = abs(norminv(pval)); % convert p-value to Z value
        n_permutes = 500; % number of permutations %1000 takes 24 hours for JZ1 wtrack days 1:6
        min_freq =  4;
        max_freq = 60;
        frexres = 2; %frequency resolution default 4
        numfrex = floor((max_freq - min_freq)/frexres);
        % set range for variable number of wavelet cycles
        % range_cycles = [min_freq*win(2) max_freq*win(2)]; % is there a principled way to do this?
        % more wavelets for the same frequency will reduce the temporal precision.
        % more wavelets for the same frequency will increase the frequency precision
        % range_cycles = [4 floor(max_freq/10)];
        range_cycles = (([min_freq max_freq])/min_freq); % sliding window

        frex = logspace(log10(min_freq),log10(max_freq),numfrex);
        % frex = linspace(min_freq,max_freq,numfrex+1);
        % nWavecycles = linspace(range_cycles(1),range_cycles(end),numfrex); %number of wavelet cycles per freq
        nWavecycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),numfrex); %number of wavelet cycles per freq
        % wp.nWavecycles = repmat(range_cycles,1,num_frex);
        half_wave_size = (length(timeWin)-1)/2;
        
        nevents = size(allNTDataCat,3);
        nWave = length(timeWin);
        nTimeSeries = 1; %keep as 1 right now even if theres a bunch of tetrodes.. used to be called 'nNTrodes'
        nsamps = size(allNTDataCat,2);
%         nData = prod(size(allNTDataCat)); %all ntrodes method
        nData = nsamps*nevents; % per ntrode method
        nConv = nWave+nData-1; % length of the result of convolution.
        nConv2pow = 2^nextpow2(nConv); %find the next pwr of 2 (for FFT optimization purposes)
        zpad2pow = nConv2pow - nConv;
        nNTrodes = size(allNTDataCat,1);
        clear allNTDataCat
        iout = whos;
        outvars = {iout(:).name};
        eval(['outvals = {', strjoin(outvars, ','), '};']);
        out = cell2struct(outvals', outvars,1);
%     case '4-30HzJustfreqs'
%         min_freq =  4;
%         max_freq = 30;
%         frexres = 2;
%         numfrex = floor((max_freq - min_freq)/frexres);
%         frex = logspace(log10(min_freq),log10(max_freq),numfrex);
    case '4-30Hz'
        srate = 1500;
        win = [-2 2];
        plotwin = [-.5 .5];
        basewin = plotwin(2); %baseline window to use for 1/f normalization
        plottimeWin = plotwin(1):1/srate:plotwin(2);
        % baseline indices. length of 1/2 plot window, centered on first plot ind
        % i.e. 0.5 sec basewindow starting at .25 sec before the start of the
        % plotwindow
        
%         baseind = ([(win(2)-plotwin(2))-basewin/2 (win(2)-plotwin(2))+basewin/2]).*srate;
        % no.. instead have the baseline be 1 second, from -1.5 : -.5 if
        % win is 2 s and plotwin is .5. 
        baseind = [.5 1] * srate;
        pval = 0.05;
        zval = abs(norminv(pval)); % convert p-value to Z value
        n_permutes = 1000;
        min_freq =  4;
        max_freq = 30;
        frexres = 2; %frequency resolution default 4
        numfrex = floor((max_freq - min_freq)/frexres);
        % set range for variable number of wavelet cycles
        range_cycles = [min_freq max_freq];
        % more wavelets for the same frequency will reduce the temporal precision.
        % more wavelets for the same frequency will increase the frequency precision
%         range_cycles = [2 floor(max_freq/10)];
%         range_cycles = (([min_freq max_freq])/min_freq); % sliding window
        timeWin = win(1):1/srate:win(2);
        frex = logspace(log10(min_freq),log10(max_freq),numfrex);
        % frex = linspace(min_freq,max_freq,numfrex+1);
        % nWavecycles = linspace(range_cycles(1),range_cycles(end),numfrex); %number of wavelet cycles per freq
        nWavecycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),numfrex); %number of wavelet cycles per freq
        % wp.nWavecycles = repmat(range_cycles,1,num_frex);
        half_wave_size = (length(timeWin)-1)/2;
        
        nevents = size(allNTDataCat,3);
        nsamps = size(allNTDataCat,2);
        nNTrodes = size(allNTDataCat,1);
        
        nWave = length(timeWin);
        nTimeSeries = 1; %keep as 1 right now even if theres a bunch of tetrodes.. used to be called 'nNTrodes'
        % nData = prod(size(allNTDataCat)); %all ntrodes method
        
        nData = nsamps*nevents; % per ntrode method
        nConv = nWave+nData-1; % length of the result of convolution.
        nConv2pow = 2^nextpow2(nConv); %find the next pwr of 2 (for FFT optimization purposes)
        zpad2pow = nConv2pow - nConv;
%     case '4-300HzJustfreqs'
%         srate = 1500;
%         min_freq =  4;
%         max_freq = 300;
%         numfrex = 25;
%         frex = logspace(log10(min_freq),log10(max_freq),numfrex);
%         baseind = [.75 1] * srate;
%         voxel_pval = .01;
%         
    case '4-300Hz'
        srate = 1500;
        win = [-2 2];
        dsamp = 2;
%         timeWin = win(1):1/srate:win(2);
%         plotwin = [-1 1];
%         plottimeWin = plotwin(1):1/srate:plotwin(2);

        basewin = [-1.5 -1]; % in seconds, period before event start to use as baseline
        prewin = [-1 -.5];
        postwin = [.5 1];
%         baseind(1,1) = dsearchn(timeWin',basetime(1));
%         baseind(1,2) = dsearchn(timeWin',basetime(2));
%         
%         baseind = [.75 1] * srate;
        pval = 0.05;
        voxel_pval = 0.01;
        mcc_voxel_pval = 0.05;
        mcc_cluster_pval = 0.05;
        zval = abs(norminv(pval)); % convert p-value to Z value
        n_permutes = 200;
        
        min_freq =  4;
        max_freq = 300;
        numfrex = 25;
        range_cycles = [3 14]; % more decreases temp, increases freq precision        
        
        frex = logspace(log10(min_freq),log10(max_freq),numfrex);
        nWavecycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),numfrex); %number of wavelet cycles per freq
%         hws = (length(timeWin)-1)/2; % half_wave_size 
        
%         nevents = size(allNTDataCat,3);
%         nsamps = size(allNTDataCat,2);
%         nNTrodes = size(allNTDataCat,1);
        
%         nWave = length(timeWin);
%         nTimeSeries = nNTrodes; %1; %keep as 1 right now even if theres a bunch of tetrodes.. used to be called 'nNTrodes'
%         nData = nsamps*nevents*nNTrodes; % per ntrode method is just nsamps*nevents
%         nConv = nWave+nData-1; % length of the result of convolution.
%         nConv2pow = 2^nextpow2(nConv); %find the next pwr of 2 (for FFT optimization purposes)
%         zpad2pow = nConv2pow - nConv; %zero padding added by nextpow2
        
    case 'riptrigMU_instantFR'
        srate = 1000;
        win = [-.5 .5]; %in seconds. The morlet wavelet needs to be centered at time zero.. so always make this window symmetric about zero
        plotwin = [-.25 .25];
        basewin = plotwin(2); %baseline window to use for 1/f normalization
        plottimeWin = plotwin(1):1/srate:plotwin(2);
        %baseline indices. length of 1/2 plot window, centered on first plot ind
        % i.e. 0.5 sec basewindow starting at .25 sec before the start of the
        % plotwindow
        baseind = ([(win(2)-plotwin(2))-basewin/2 (win(2)-plotwin(2))+basewin/2]).*srate;
        % baseind = [1 abs(win(1))*srate-(abs(plotwin(1))*srate)]; %the entire duration before the plotwin
        % basedur = abs(win(1))+abs(win(2));
        % baselineind = ([1 (srate*basedur)+1]);
        pval = 0.05;
        zval = abs(norminv(pval)); % convert p-value to Z value
        n_permutes = 100; % number of permutations %1000 takes 24 hours for JZ1 wtrack days 1:6
        min_freq =  5;
        max_freq = 8;
        frexres = .5; %frequency resolution default 4
        numfrex = floor((max_freq - min_freq)/frexres);
        % set range for variable number of wavelet cycles
        % range_cycles = [min_freq*win(2) max_freq*win(2)]; % is there a principled way to do this?
        % more wavelets for the same frequency will reduce the temporal precision.
        % more wavelets for the same frequency will increase the frequency precision
        % range_cycles = [4 floor(max_freq/10)];
        range_cycles = (([min_freq max_freq])/min_freq); % sliding window
        timeWin = win(1):1/srate:win(2);
        frex = logspace(log10(min_freq),log10(max_freq),numfrex);
        % frex = linspace(min_freq,max_freq,numfrex+1);
        % nWavecycles = linspace(range_cycles(1),range_cycles(end),numfrex); %number of wavelet cycles per freq
        nWavecycles = logspace(log10(range_cycles(1)),log10(range_cycles(end)),numfrex); %number of wavelet cycles per freq
        % wp.nWavecycles = repmat(range_cycles,1,num_frex);
        half_wave_size = (length(timeWin)-1)/2;
        
        nevents = size(allNTDataCat,3);
        nWave = length(timeWin);
        nTimeSeries = 1; %keep as 1 right now even if theres a bunch of tetrodes.. used to be called 'nNTrodes'
        nsamps = size(allNTDataCat,2);
%         nData = prod(size(allNTDataCat)); %all ntrodes method
        nData = nsamps*nevents; % per ntrode method
        nConv = nWave+nData-1; % length of the result of convolution.
        nConv2pow = 2^nextpow2(nConv); %find the next pwr of 2 (for FFT optimization purposes)
        zpad2pow = nConv2pow - nConv;
        nNTrodes = size(allNTDataCat,1);
%         clear allNTDataCat
%         iout = whos;
%         outvars = {iout(:).name};
%         eval(['outvals = {', strjoin(outvars, ','), '};']);
%         out = cell2struct(outvals', outvars,1);
    otherwise
        error('pick an existing wave param set or make a new one')
end
clear allNTDataCat
iout = whos;
outvars = {iout(:).name};
eval(['outvals = {', strjoin(outvars, ','), '};']);
out = cell2struct(outvals', outvars,1);
end