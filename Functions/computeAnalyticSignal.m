
function as = computeAnalyticSignal(lfpstack, varargin)
% compute the analytic signal for the input rip triggered LFP

% lfpstack is a struct array, row per an with fields data, ntrodes, lfptypes, animal
% lfpstack(ian).data is a cell array of eeg types with ntrode x sample x trial

% Demetris Roumis

saveAnalyticSignal = 1;
overwrite = 0;
useeeg = 'eeg';
waveSet = '4-300Hz';
if ~isempty(varargin)
    assign(varargin{:});
end

for ian = 1:length(lfpstack)
    eegidx = find(cellfun(@(x) strcmp(x,useeeg),lfpstack(ian).lfptypes, 'un', 1));
    animal = lfpstack(ian).animal;
    ntrodes = lfpstack(ian).ntrodes;
    fprintf('waveSet: %s\n', waveSet);
    wp = getWaveParams(waveSet, lfpstack(ian).data{eegidx});
    
    % create complex morlet wavelets and get FFT
    waveletFFT = createCMWfft(wp);
    pathdef = animaldef(lower('Demetris'));
    dirstr = sprintf('%s/analyticSignal/%s', pathdef{2}, animal);
    
    parfor nt = 1:length(ntrodes) %use parfor
        ASsavestr = sprintf('%s/nt%02d_waveSet-%s_AS.mat',dirstr, nt, wp.waveSet);
        if exist(ASsavestr, 'file') && ~overwrite
            fprintf('AS results detected, skipping ++++++++++ %s \n',ASsavestr)
            continue
        end
        
        fprintf('computing Analytic Signal for ntrode %d of %d \n',nt,wp.nNTrodes);
        
        dataY = zeros(wp.nConv2pow, wp.nTimeSeries);
        astmp = zeros(wp.nConv2pow, wp.numfrex);
        astmp2 = zeros(wp.numfrex, wp.nConv2pow);
        astmp3 = zeros(wp.numfrex, wp.nData);
        as = zeros(wp.nsamps, wp.nevents, wp.numfrex);
        
        % get pow2 padded FFT of nt data
        tmp1 = squeeze(datastack(nt,:,:));
        tmp2 = reshape(tmp1,wp.nData,wp.nTimeSeries);
        dataY = fft(tmp2,wp.nConv2pow,1); %fft works on columns
        
        % run convolution (filtering) : Time-domain convolution in the frequency domain...
        astmp = bsxfun(@times,dataY,waveletFFT); %multiply the pwr spectrum (FFT) of the data and the wavelet
        %  ifft returns the inverse DFT of each column of the AS matrix
        astmp2 = ifft(astmp,wp.nConv2pow,1); % take the inverse transform of convolution
        % trim half length wavelet and zero padding at the tail
        astmp3 = astmp2(wp.half_wave_size+1:end-(wp.half_wave_size+wp.zpad2pow),:);
        % reshape back to nt x samples x events x frequencies
        as = reshape(astmp3,wp.nsamps,wp.nevents,wp.numfrex);
        %     ph =  angle(as);
        
        %% ---------------- Save Analytic Signal Output ---------------------------------------------------
        if saveAnalyticSignal == 1
            saveAS(as,ntrodes(nt),nt,wp,animal)
        end
    end
end
end

function waveletFFT=createCMWfft(wp)
% create wavelet FFT
waveletFFT = zeros(wp.nConv2pow, wp.numfrex);
parfor fi=1:wp.numfrex % can use parfor
    % creating the morlet wavelet by combining the complex sine wave and the gaussian
    sine_wave = exp(2*1i*pi*wp.frex(fi).*wp.timeWin); %make a complex sine wave
    bn = wp.nWavecycles(fi)/(2*pi*wp.frex(fi)); % std of gaussian of wavelet (dependent on the freq and #cycles)
    gaus_win = exp(-wp.timeWin.^2./(2*bn^2)); % gaussian
    wavelet = sine_wave .* gaus_win;
    wavefft = fft(wavelet,wp.nConv2pow); % take the fft of the wavelet pad with zeros to the next power of 2 for speed
    % normalize wavelet to a maximum of 1 to ensure convolution units are same as data
    waveletFFT(:,fi) = (wavefft ./ max(wavefft))';
    fprintf('creating wavelet for freq %d of %d \n',fi,wp.numfrex);
end
end
