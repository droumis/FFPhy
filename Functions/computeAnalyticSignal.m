
function as = computeAnalyticSignal(lfpstack, varargin)
% compute the analytic signal for the input rip triggered LFP

% lfpstack is a struct array, row per an with fields data, ntrodes, lfptypes, animal
% lfpstack(ian).data is a cell array of eeg types with ntrode x sample x trial

% Demetris Roumis 2019

saveAnalyticSignal = 1;
overwrite = 1;
uselfptypes = {'eeg'};
waveSet = '4-300Hz';

if ~isempty(varargin)
    assign(varargin{:});
end

for ian = 1:length(lfpstack)
    wp = getWaveParams(waveSet);
    for itype = 1:length(uselfptypes)
        phaseout = struct;
        powerout = struct;
        lfptype = uselfptypes{itype};
        itypeidx = find(strcmp(lfpstack(ian).lfptypes, lfptype));
        animal = lfpstack(ian).animal;
        ntrodes = lfpstack(ian).ntrodes;
        fprintf('waveSet: %s\n', waveSet);

        srate = wp.srate/wp.dsamp;
        dsampdata = lfpstack(ian).data{itypeidx}(:,1:wp.dsamp:end,:);
        nevents = size(dsampdata,3);
        nsamps = size(dsampdata,2);
        nNTrodes = size(dsampdata,1);
        nData = nsamps*nevents*nNTrodes;
        
        timeWin = wp.win(1):1/srate:wp.win(2);
        baseind(1,1) = dsearchn(timeWin',wp.basewin(1));
        baseind(1,2) = dsearchn(timeWin',wp.basewin(2));
        lenWave = length(timeWin);
        hws = (length(timeWin)-1)/2; % half_wave_size
        
        nConv = lenWave+nData-1; % length of the result of convolution
        nConv2pow = 2^nextpow2(nConv); % next pwr of 2 for FFT speedup
        zpad2pow = nConv2pow - nConv; % zero padding added by nextpow2

        tic
        waveletFFT = createCMWfft(wp, nConv2pow, timeWin);
        fprintf('%d seconds to make CMWFFT\n', toc);
        
        fprintf('computing Analytic Signal for %s %s\n', animal, lfptype);        
        dataFFT = zeros(nConv2pow, 1);
        astmp = zeros(nConv2pow, wp.numfrex);
        as = zeros(nNTrodes, nsamps, nevents, wp.numfrex);
        
        dataFFT = fft(reshape(dsampdata,nData,1),nConv2pow,1);
        % Time-domain convolution in the frequency domain
        astmp = bsxfun(@times,dataFFT,waveletFFT);
        astmp = ifft(astmp,nConv2pow,1);
        % trim and reshape to 'source', 'sample', 'trial', 'frequency'
        as = reshape(astmp(hws+1:end-(hws+zpad2pow),:), ...
            nNTrodes, nsamps, nevents, wp.numfrex);
        % Phase
        phaseout.ph =  angle(as);
        phaseout.wp = wp;
        phaseout.animal = animal;
        phaseout.lfptype = lfptype;
        phaseout.dsampsrate = srate;
        phaseout.dsamp = wp.dsamp;
        phaseout.dims = {'source', 'sample', 'trial', 'frequency'};
        % Power
        powerout.pwr =  abs(as).^2;
        powerout.wp = wp;
        powerout.animal = animal;
        powerout.lfptype = lfptype;
        powerout.dsampsrate = srate;
        powerout.dsamp = wp.dsamp;
        powerout.dims = {'source', 'sample', 'trial', 'frequency'};
        %% ---------------- Save Analytic Signal Output ---------------------------------------------------
        if saveAnalyticSignal == 1
            andef = animaldef(lower('Demetris'));
            save_data(phaseout, sprintf('%s/analyticSignal/', andef{2}), ...
                sprintf('AS_waveSet-%s_%s_phase', waveSet, lfptype))
            save_data(powerout, sprintf('%s/analyticSignal/', andef{2}), ...
                sprintf('AS_waveSet-%s_%s_power', waveSet, lfptype))
        end
    end
end
end
