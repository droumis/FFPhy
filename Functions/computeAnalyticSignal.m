
function [powerout, phaseout] = computeAnalyticSignal(lfpstack, varargin)
% compute the analytic signal for the input rip triggered LFP

% lfpstack is a struct array per animal: data, ntrodes, lfptypes, animal
% lfpstack(ian).data = {eegtype}(ntrode x sample x ripple mat)

% power and phase : (ntrode x sample x ripple x frequency)

% Demetris Roumis 2019

saveOutput = 1;
uselfptype = 'eeg';
waveSet = '4-300Hz';

if ~isempty(varargin)
    assign(varargin{:});
end

for ian = 1:length(lfpstack)
    animal = lfpstack(ian).animal;
    fprintf('computing Analytic Signal for %s %s\n', animal, uselfptype);
    fprintf('waveSet: %s\n', waveSet);
    wp = getWaveParams(waveSet);    
    itypeidx = find(strcmp(lfpstack(ian).lfptypes, uselfptype));
    dsampdata = lfpstack(ian).data{itypeidx}(:,1:wp.dsamp:end,:);
    nevents = size(dsampdata,3);
    nsamps = size(dsampdata,2);
    nNtrodes = size(dsampdata,1);
    nData = nsamps*nevents*nNtrodes;
    timeWin = wp.win(1):1/(wp.srate/wp.dsamp):wp.win(2);
    baseind(1,1) = dsearchn(timeWin',wp.basewin(1));
    baseind(1,2) = dsearchn(timeWin',wp.basewin(2));
    lenWave = length(timeWin);
    hws = (length(timeWin)-1)/2; % half_wave_size    
    nConv = lenWave+nData-1; % length of the result of convolution
    nConv2pow = 2^nextpow2(nConv); % next pwr of 2 for FFT speedup
    zpad2pow = nConv2pow - nConv; % zero padding added by nextpow2

    waveletFFT = createCMWfft(wp, nConv2pow);
    
    tic
    dataFFT = zeros(nConv2pow, 1);
    astmp = zeros(nConv2pow, wp.numfrex);
    as = zeros(nsamps, nNtrodes, nevents, wp.numfrex);
     % flatten all the data, then fft. i make time as dim 1 because rows stack first
    dataFFT = fft(reshape(permute(dsampdata, [2 1 3]),[],1),nConv2pow);
    astmp = bsxfun(@times,dataFFT,waveletFFT);% time convolution in freq domain
    astmp = ifft(astmp,nConv2pow);
    % trim length wavlet, padding and reshape to ntrode x sample x ripple x frequency
    as = reshape(astmp(hws+1:end-(hws+zpad2pow),:), nsamps, nNtrodes, nevents, wp.numfrex);
    as = permute(as, [2 1 3 4]);
    
    % Phase
    fprintf('getting phase \n');
    phaseout.ph =  single(angle(as)); % convert to single to save diskspace
    phaseout.wp = wp;
    phaseout.animal = animal;
    phaseout.day = lfpstack(ian).day;
    phaseout.epoch = lfpstack(ian).epoch;
    phaseout.ripStartTime = lfpstack(ian).ripStartTime;
    phaseout.ripEndTime = lfpstack(ian).ripEndTime;
    phaseout.ntrode = lfpstack(ian).ntrodes;
    phaseout.frequency = wp.frex;
    phaseout.lfptype = uselfptype;
    phaseout.dsampsrate = wp.win(1):1/(wp.srate/wp.dsamp):wp.win(2);
    phaseout.dsamp = wp.dsamp;
    phaseout.time = timeWin;
    phaseout.dims = {'ntrode', 'sample', 'ripple', 'frequency'};
    
    % Power
    fprintf('getting power \n');
    powerout.pwr =  single(abs(as).^2);
    powerout.wp = wp;
    powerout.animal = animal;
    powerout.day = lfpstack(ian).day;
    powerout.epoch = lfpstack(ian).epoch;
    powerout.ripStartTime = lfpstack(ian).ripStartTime;
    powerout.ripEndTime = lfpstack(ian).ripEndTime;
    powerout.ntrode = lfpstack(ian).ntrodes;
    powerout.frequency = wp.frex;
    powerout.lfptype = uselfptype;
    powerout.dsampsrate = wp.win(1):1/(wp.srate/wp.dsamp):wp.win(2);
    powerout.dsamp = wp.dsamp;
    powerout.time = timeWin;
    powerout.dims = {'ntrode', 'sample', 'ripple', 'frequency'};
    fprintf('%.02f seconds to compute pwr,ph AS\n', toc);
    
    %% ---------------- Save Output ---------------------------------------------------
    if saveOutput == 1
        fprintf('saving AS\n')
        andef = animaldef(lower('Demetris'));
        save_data(phaseout, sprintf('%s/analyticSignal/', andef{2}), ...
            sprintf('AS_waveSet-%s_%s_phase', waveSet, uselfptype))
        save_data(powerout, sprintf('%s/analyticSignal/', andef{2}), ...
            sprintf('AS_waveSet-%s_%s_power', waveSet, uselfptype))
    end
end
end

function waveletFFT=createCMWfft(wp, nConv2pow)
% create complex morlet wavelets and return FFT
tic
waveletFFT = zeros(nConv2pow, wp.numfrex);
% fwhm = linspace(.8, .7, wp.numfrex);
timeWin = wp.win(1):1/(wp.srate/wp.dsamp):wp.win(2);
parfor fi=1:wp.numfrex % can use parfor
    sine_wave = exp(2*1i*pi*wp.frex(fi).*timeWin); % complex sine wave
    bn = wp.nWavecycles(fi)/(2*pi*wp.frex(fi)); % std of gaussian, dependent on freq and #cycles
    gaus_win = exp(-timeWin.^2./(2*bn^2)); % gaussian
    wavelet = sine_wave .* gaus_win;
%     gaus_fwhm = exp(-(4*log(2)*timeWin).^2/fwhm(fi).^2);
    wavefft = fft(wavelet,nConv2pow); % fft of the wavelet, pad to next power of 2 
    % normalize wavelet to a maximum of 1 to ensure convolution units are same as data
    waveletFFT(:,fi) = (wavefft ./ max(wavefft))';
    fprintf('creating wavelet for freq %d of %d \n',fi,wp.numfrex);
end
fprintf('%.02f seconds to make CMWFFT\n', toc);
end