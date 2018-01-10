
function computeAnalyticSignal(allNTDataCat, wp, animal, resultfilename, varargin)

saveAnalyticSignal = 1;
if ~isempty(varargin)
    assign(varargin{:});
end
waveletFFT = zeros(wp.nConv2pow, wp.numfrex);
parfor fi=1:wp.numfrex % can use parfor
    % create wavelet and get its FFT
    % creating the morlet wavelet by combingin the complex sine wave and the gaussian
    sine_wave = exp(2*1i*pi*wp.frex(fi).*wp.timeWin); %make a complex sine wave
    bn = wp.nWavecycles(fi)/(2*pi*wp.frex(fi)); % std of the gaussian of the morlet wavelet. dependent on the freq and #cycles
    gaus_win = exp(-wp.timeWin.^2./(2*bn^2)); %the gaussian
    wavelet = sine_wave .* gaus_win;
    wavefft = fft(wavelet,wp.nConv2pow); % take the fft of the wavelet pad with zeroes to the next pwer of 2 for speed
    %                 waveletFFT{introde}{fi} = fft(exp(2*1i*pi*frex(fi).*wp.timeWin) .* exp(-wp.timeWin.^2./(2*wp.nWavecycles(fi)/(2*pi*frex(fi))^2)),nConv2pow);
    %normalize wavelet to a maximum of 1. this will ensure that the units of convolution are the same as in the original data.
    waveletFFT(:,fi) = (wavefft ./ max(wavefft))';
    disp(sprintf('creating wavelet for freq %d of %d',fi,wp.numfrex));
end

animalinfo = animaldef(lower(animal));
animalID = animalinfo{1,3}; %use anim prefix for name
FFanimdir =  sprintf('%s',animalinfo{1,2});
dtypedir = 'analyticSignal';
dirstr = sprintf('%s%s/', FFanimdir, dtypedir);
if ~isdir(dirstr);
    mkdir(dirstr);
end
parfor nt = 1:wp.nNTrodes
    dataY = zeros(wp.nConv2pow, wp.nTimeSeries);
    astmp = zeros(wp.nConv2pow, wp.numfrex);
    astmp2 = zeros(wp.numfrex, wp.nConv2pow);
    astmp3 = zeros(wp.numfrex, wp.nData);
    
    as = zeros(wp.nsamps, wp.nevents, wp.numfrex);
    ph = zeros(wp.nsamps, wp.nevents, wp.numfrex);
    disp(sprintf('computing Analytic Signal for ntrode %d of %d',nt,wp.nNTrodes));
    tmp1 = squeeze(allNTDataCat(nt,:,:));
    tmp2 = reshape(tmp1,wp.nData,wp.nTimeSeries);
    dataY = fft(tmp2,wp.nConv2pow,1); %fft works on columns
    %                 waveletFFT{introde}{fi} = waveletFFT{introde}{fi} ./ max(waveletFFT{introde}{fi});
    % run convolution (filtering) : Time-domain convolution in the frequency domain... because itz so much faster
    astmp = bsxfun(@times,dataY,waveletFFT); %multiply the pwr spectrum (FFT) of the data and the wavelet
    %                 astmp{introde}{fi} = bsxfun(@times,dataY{introde},waveletFFT); %multiply the pwr spectrum of the data and the wavelet
    %  If X is a matrix, ifft returns the inverse DFT of each column of the matrix.
    astmp2 = ifft(astmp,wp.nConv2pow,1); % take the inverse transform
    %trim off the length of half the wavelet at the beginning
    %and at the end. also trim the zero padding at the tail
%     astmp3 = astmp2(wp.half_wave_size+1:end-wp.half_wave_size-wp.zpad2pow,:);
    astmp3 = astmp2(wp.half_wave_size+1:end-(wp.half_wave_size+wp.zpad2pow),:);
    
    as = reshape(astmp3,wp.nsamps,wp.nevents,wp.numfrex); %reshape it back to ntrodes by nsamples by nevents by frequencies
    ph =  angle(as); %this takes sooooooooo long to do 'as' all at once and eats memory.. do inside the loop instead
    
    %% ---------------- Save Analytic Signal Output ---------------------------------------------------
    if saveAnalyticSignal == 1;
        saveAS(as, ph, dirstr, nt, wp)
%         ASsavestr = sprintf('%s/nt%02d_waveSet-%s_AS.mat',dirstr, nt, wp.waveSet);
%         save(ASsavestr, 'as', '-v7.3');
%         disp(sprintf('SAVED ANALYTIC SIGNAL RESULTS ++++++++++ %s',ASsavestr))
%         PHsavestr = sprintf('%s/nt%02d_waveSet-%s_PH.mat',dirstr, nt, wp.waveSet);
%         save(PHsavestr, 'ph', '-v7.3');
%         disp(sprintf('SAVED PHASE RESULTS ++++++++++ %s',PHsavestr))
    end
end

%         for intr = 1:ixpc.wp.nNTrodes %  par %make this run on all ntrodes at once, then loop ntrodes to do the perm tests
%             %             intDataCat = squeeze(ixpc.allNTDataCat{ian}(int,:,:)); %squeeze will make nsampes x wp.nevents..
%             intDataCat = squeeze(allNTDataCat(intr,:,:)); %squeeze will make nsampes x wp.nevents..
%             % concat reshape all the events (peri-rip snips) for speed and to reduce edge artifacts, then take the fft.
%             % fft with the arg of wp.nConv ensures that matlab will zero-pad the half-wavelet duration on either side of the data series
%             % dataY = fft(reshape(iEpTetDataCat,wp.nNTrodes,wp.nData),wp.nConv,2);
%             dataY = fft(reshape(intDataCat,ixpc.wp.nTimeSeries,ixpc.wp.nData),ixpc.wp.nConv2pow,2); % reshape reshapes row-wise into a 1 x d vector
%             [ixpc.as{ian}{int}, ixpc.ph{ian}{int}] = computeAnalyticSignal(dataY, ixpc.wp); %ph, as ~ {int}
%         end

end
