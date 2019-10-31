
function [powerout, phaseout] = computeAnalyticSignal(lfpstack, varargin)
% compute the analytic signal for the input rip triggered LFP

% lfpstack is a struct array per animal: data, ntrodes, lfptypes, animal
% lfpstack(ian).data = {eegtype}(ntrode x sample x ripple mat)

% power and phase : (ntrode x sample x ripple x frequency)

% Demetris Roumis 2019

saveOutput = 1;
lfptype = 'eeg';
waveSet = '4-300Hz';
env = 'ep';
if ~isempty(varargin)
    assign(varargin{:});
end

for ian = 1:length(lfpstack)
    animal = lfpstack(ian).animal;
    fprintf('computing Analytic Signal for %s %s\n', animal, lfptype);
    fprintf('waveSet: %s\n', waveSet);
    wp = getWaveParams(waveSet);
    t = find(strcmp(lfpstack(ian).lfptypes, lfptype));
    
    dsampdata = lfpstack(ian).data{t}(:,1:wp.dsamp:end,:);
    nNtrodes = size(dsampdata,1);
    nevents = size(dsampdata,3);
    nsamps = size(dsampdata,2);
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
    dataFFT = fft(reshape(permute(dsampdata, [2 1 3]),nData,1),nConv2pow);
    astmp = bsxfun(@times,dataFFT,waveletFFT);% time convolution in freq domain
    astmp = ifft(astmp,nConv2pow);
    % trim length wavlet, padding and reshape to ntrode x sample x ripple x frequency
    as = reshape(astmp(hws+1:end-(hws+zpad2pow),:), nsamps, nNtrodes, nevents, wp.numfrex);
    as = permute(as, [2 1 3 4]);
    
    % Phase
    fprintf('getting phase \n');
    phaseout(ian).ph =  single(angle(as)); % convert to single to save diskspace
    phaseout(ian).wp = wp;
    phaseout(ian).animal = animal;
    phaseout(ian).day = lfpstack(ian).day{t};
    phaseout(ian).epoch = lfpstack(ian).epoch{t};
    phaseout(ian).evStart = lfpstack(ian).evStart{t};
    phaseout(ian).evEnd = lfpstack(ian).evEnd{t};
    phaseout(ian).ntrode = lfpstack(ian).ntrodes{t};
    phaseout(ian).frequency = wp.frex;
    phaseout(ian).lfptype = lfptype;
    phaseout(ian).dsampsrate = wp.win(1):1/(wp.srate/wp.dsamp):wp.win(2);
    phaseout(ian).dsamp = wp.dsamp;
    phaseout(ian).time = timeWin;
    phaseout(ian).dims = {'ntrode', 'sample', 'event', 'frequency'};
    
    % Power
    fprintf('getting power \n');
    powerout(ian).pwr =  single(abs(as).^2);
    powerout(ian).wp = wp;
    powerout(ian).animal = animal;
    powerout(ian).day = lfpstack(ian).day{t};
    powerout(ian).epoch = lfpstack(ian).epoch{t};
    powerout(ian).evStart = lfpstack(ian).evStart{t};
    powerout(ian).evEnd = lfpstack(ian).evEnd{t};
    powerout(ian).ntrode = lfpstack(ian).ntrodes{t};
    powerout(ian).frequency = wp.frex;
    powerout(ian).lfptype = lfptype;
    powerout(ian).dsampsrate = wp.win(1):1/(wp.srate/wp.dsamp):wp.win(2);
    powerout(ian).dsamp = wp.dsamp;
    powerout(ian).time = timeWin;
    powerout(ian).dims = {'ntrode', 'sample', 'event', 'frequency'};
    fprintf('%.02f seconds to compute AS, pwr, ph \n', toc);
    
    %% ---------------- Save Output ---------------------------------------------------
    if saveOutput == 1
        fprintf('saving AS\n')
        andef = animaldef(lower('Demetris'));
        save_data(phaseout(ian), sprintf('%s/analyticSignal/', andef{2}), ...
            sprintf('AS_waveSet-%s_%s_%s_phase', waveSet, lfptype, env))
        save_data(powerout(ian), sprintf('%s/analyticSignal/', andef{2}), ...
            sprintf('AS_waveSet-%s_%s_%s_power', waveSet, lfptype, env))
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
% i think this will ensure that the wavelett and data fft have same
    % freq resolution
    wavefft = fft(wavelet,nConv2pow); % fft of the wavelet, pad to next power of 2 of data
    
    % normalize wavelet to a maximum of 1 to ensure convolution units are same as data
    waveletFFT(:,fi) = (wavefft ./ max(wavefft))';
    fprintf('creating wavelet for freq %d of %d \n',fi,wp.numfrex);
end
fprintf('%.02f seconds to make CMWFFT\n', toc);
end