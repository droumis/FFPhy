
function [powerout, phaseout] = computeAnalyticSignal(lfpstack, varargin)
% compute the analytic signal for a stack of LFP clips, e.g. swr centered.
% lfpstack is a struct array per animal: data, ntrodes, lfptypes, animal
% lfpstack(ian).data = {eegtype}(ntrode x sample x ripple mat)
% power and phase : (ntrode x sample x ripple x frequency)
%
%            Mendacious Mushroom
%                 __.....__
%              .'" _  o    "`.
%            .' O (_)     () o`.
%           .           O       .
%          . ()   o__...__    O  .
%         . _.--"""       """--._ .
%         :"                     ";
%          `-.__    :   :    __.-'
%               """-:   :-"""
%                  J     L
%                  :     :
%                 J       L
%                 :       :
%                 `._____.' 
%
%{ 
Notes:
- forest:bear:cactus:mushroom:beer:leaf

FFPhy V0.1
@DR
%}


saveOutput = 1;
lfptype = 'eeg';
waveSet = '4-300Hz';
env = 'ep';
eventType = 'swr';
if ~isempty(varargin)
    assign(varargin{:});
end

for ian = 1:length(lfpstack)
    tic
    try
        animal = lfpstack(ian).animal{3};
    catch 
        animal = lfpstack(ian).animal;
    end
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
    %     baseind(1,1) = dsearchn(timeWin',wp.basewin(1));
    %     baseind(1,2) = dsearchn(timeWin',wp.basewin(2));
    lenWave = length(timeWin);
    hws = (length(timeWin)-1)/2; % half_wave_size
    nConv = lenWave+nData-1; % length of the result of convolution
    nConv2pow = 2^nextpow2(nConv); % next pwr of 2 for FFT speedup
    zpad2pow = nConv2pow - nConv; % zero padding added by nextpow2
    
%     waveletFFT = zeros(nConv2pow, wp.numfrex);
    as = zeros(nNtrodes, nsamps, nevents, wp.numfrex);
    for fi=1:wp.numfrex % can use parfor (but be currfull)
        waveletFFT = createCMWfft(wp.frex(fi), wp.nWavecycles(fi), wp.win, wp.srate,...
            wp.dsamp, nConv2pow);
    %     dataFFT = zeros(nConv2pow, 1);
    %     astmp = zeros(nConv2pow, wp.numfrex);
    %     as = zeros(nsamps, nNtrodes, nevents, wp.numfrex);
        % flatten all the data, then fft. i make time as dim 1 because rows stack first
        astmp = fft(reshape(permute(dsampdata, [2 1 3]),nData,1),nConv2pow);
        astmp = bsxfun(@times,astmp,waveletFFT);% time convolution in freq domain
        tic
        astmp = ifft(astmp,nConv2pow);
        fprintf('running ifft took %.02f s\n', toc);
        % trim length wavlet, padding and reshape to ntrode x sample x ripple x frequency
        astmp = reshape(astmp(hws+1:end-(hws+zpad2pow),:), nsamps, nNtrodes, nevents, 1);
        as(:,:,:,fi) = permute(astmp, [2 1 3 4]);
    end
    % Phase
    fprintf('getting phase \n');
    phaseout(ian).ph =  single(angle(as)); % convert to single to save diskspace
    phaseout(ian).wp = wp;
    phaseout(ian).animal = animal;
    phaseout(ian).day = lfpstack(ian).day{t};
    phaseout(ian).epoch = lfpstack(ian).epoch{t};
    phaseout(ian).evStart = lfpstack(ian).evStart{t};
    %     phaseout(ian).evEnd = lfpstack(ian).evEnd{t};
    phaseout(ian).ntrode = lfpstack(ian).ntrodes{t};
    phaseout(ian).frequency = wp.frex;
    phaseout(ian).lfptype = lfptype;
    phaseout(ian).dsampsrate = wp.win(1):1/(wp.srate/wp.dsamp):wp.win(2);
    phaseout(ian).dsamp = wp.dsamp;
    phaseout(ian).time = timeWin;
    phaseout(ian).dims = {'ntrode', 'sample', 'event', 'frequency'};
    
    % Power
    fprintf('getting power \n');
    powerout(ian).pwr =  single(abs(as).^2); % convert to single to save diskspace
    powerout(ian).wp = wp;
    powerout(ian).animal = animal;
    powerout(ian).day = lfpstack(ian).day{t};
    powerout(ian).epoch = lfpstack(ian).epoch{t};
    powerout(ian).evStart = lfpstack(ian).evStart{t};
    %     powerout(ian).evEnd = lfpstack(ian).evEnd{t};
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
        fprintf('saving\n')
        andef = animaldef(lower('Demetris'));
        save_data(phaseout(ian), sprintf('%s/analyticSignal/', andef{2}), ...
            sprintf('LFPphase_%s_%s_%s_%s', waveSet, lfptype, env, eventType))
        save_data(powerout(ian), sprintf('%s/analyticSignal/', andef{2}), ...
            sprintf('LFPpower_%s_%s_%s_%s', waveSet, lfptype, env, eventType))
    end
end
end

function wavefft = createCMWfft(freq, nWaves, win, srate, dsamp, nConv2pow)
% create complex morlet wavelets and return FFT
tic
% waveletFFT = zeros(nConv2pow, 1);
% gaus_fwhm = exp(-(4*log(2)*timeWin).^2/fwhm(fi).^2); % ignore this
timeWin = win(1):1/(srate/dsamp):win(2);
sine_wave = exp(2*1i*pi*freq.*timeWin); % complex sine wave
bn = nWaves/(2*pi*freq); % std of gaussian, dependent on freq and #cycles
gaus_win = exp(-timeWin.^2./(2*bn^2)); % gaussian
wavelet = sine_wave .* gaus_win;
% ensure that the wavelett and data fft have same freq resolution
wavefft = fft(wavelet,nConv2pow); % fft of the wavelet, pad to next power of 2 of data
% normalize wavelet to a maximum of 1 to ensure convolution units are same as data
wavefft = (wavefft ./ max(wavefft))';
fprintf('creating waveletFFT for freq %.02f took %.02f s\n', freq, toc);
end