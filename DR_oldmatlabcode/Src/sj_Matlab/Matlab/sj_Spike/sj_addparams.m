function [addn_params, new_waves] = sj_addparams(waves, triggers, dopc)
%[addn_params] = sj_addparams (waves, triggers,dopc)

% 20 May 2011. Shantanu
% Caleb's Peak ALign by C-M. Using Over-sampled waveform to get additional
% parameters. Then Peak Aligning and getting PC if asked for.
% function pp = peak_align(waves)

if nargin<2,
    triggers = 40*ones(size(waves,2),1);
end
if nargin<3
    dopc = 0;
end

dopc;

% FIXED PARAMS
Fs = 31.2; % kHz
OVERSAMPLING = 4;
SNIPPET_LENGTH = 28;
SNIPPET_PRE = 7;

THRESHOLD_IND = 9;
WINDOW = [5,15];
CM_THRESH = triggers; % Can use also
CM_THRESH = 20;
WINDOW_VALLEY_END = 20;

N = size(waves,1);
K = N;
L = length(((WINDOW(1))*OVERSAMPLING) : (WINDOW(2)*OVERSAMPLING));

% generate oversampling matrix
delta = zeros(K,1);
for i = 1:K
    dd = delta;
    dd(i) = 1;
    fw = fft(dd);
    fo = zeros(K*OVERSAMPLING,1);
    fo(1:(K/2)) = fw(1:(K/2));
    fo((end-(K/2-2)):end) = fw(((K/2)+2):end);
    OM(:,i) = OVERSAMPLING*ifft(fo);
end

% Initialize before going into loop
wU = nan(N*OVERSAMPLING, size(waves,2));
new_waves = nan(SNIPPET_LENGTH,4,size(waves,3));

% Additional Parameters:1)Peak-Val RatioxNCb, 2)Peak-Val Latency, 3)Width-FWHM, 4)Valley, 5)Valley a, 6) 3xNch PCs, usually 12
if dopc==1
    addn_params = nan(size(waves,3),(5+3)*size(waves,2)); % 3 PCs for each channel
else
    addn_params = nan(size(waves,3),5*size(waves,2));
end
%center_of_mass = nan(1,size(waves,3));
%ind_startret = nan(1,size(waves,3));

threshold = min(max(waves(THRESHOLD_IND,:,:),[],2)); % this tells us the approximate threshold

%tic;
for i = 1:size(waves,3)
    
    if mod(i,100000)==0
        disp(['      Addn Params: ' num2str(i) ' spikes done']);
        %toc;
        %tic;
    end
    
    w = waves(:,:,i);
    wU = OM*double(w);
    
    % Oversampled Waveform - Get Width, Peak-Val Ratio, Peak-Val Latency, Valley, Valley a
    % Take Max Window 9-15,  For Valley, Find 1st Min Peak in Window: PeakIdx-20
    %[m,im]=max(wU(THRESHOLD_IND*OVERSAMPLING : WINDOW(2)*OVERSAMPLING),[],1);
    
    for nch=1:size(waves,2)
        
        % Peak
        [m,im]=max(wU( (THRESHOLD_IND*OVERSAMPLING) : (WINDOW(2)*OVERSAMPLING),nch ));
        peak(nch) = m(1);
        peakidx(nch) = im(1)+(THRESHOLD_IND*OVERSAMPLING)-1;
        
        % Valley (Parameter b) - After Peakidx
        [m,im] = sj_peakpick(wU( peakidx(nch) : (WINDOW_VALLEY_END*OVERSAMPLING), nch),'min');
        if ~isempty(m)
            valley(nch) = m(1);
            valleyidx(nch) = im(1)+peakidx(nch)-1;
        else
            [m,im] = min(wU( peakidx(nch) : (WINDOW_VALLEY_END*OVERSAMPLING), nch));
            valley(nch) = m(1);
            valleyidx(nch) = im(1)+peakidx(nch)-1;
        end
        
        % Valley Pre-Threshold (Parameter a)
        [m,im] = sj_peakpick(wU( 1:(THRESHOLD_IND*OVERSAMPLING), nch),'min');
        if ~isempty(m)
            valleya(nch) = m(end);
            valleyaidx(nch) = im(end);
        else
            [m,im] = min(wU(1:(THRESHOLD_IND*OVERSAMPLING), nch));
            valleya(nch) = m(end);
            valleyaidx(nch) = im(end);
        end
        
        % Width - Full Width at Max(Half-Max, Trigger Thresh); No cannot use trigger thresh
        %halfmax = max(peak(nch)/2,triggers(nch));
        %idxs = lookup(wU( (WINDOW(1)*OVERSAMPLING) : peakidx(nch),nch), halfmax);
        
        halfmax = max(peak(nch)/2,0);
        idx1 = min( find( wU((WINDOW(1)*OVERSAMPLING):peakidx(nch)-1,nch)<=halfmax & wU((WINDOW(1)*OVERSAMPLING)+1:peakidx(nch),nch)>halfmax ) );
        if isempty(idx1)
            idx1 = WINDOW(1)*OVERSAMPLING;
        else
            idx1 = idx1 + (WINDOW(1)*OVERSAMPLING) - 1;
        end
        idx2 = min( find( wU(peakidx(nch):(WINDOW_VALLEY_END*OVERSAMPLING)-1,nch)>=halfmax & wU(peakidx(nch)+1:(WINDOW_VALLEY_END*OVERSAMPLING),nch)<halfmax ) )+1;
        if isempty(idx2)
            idx2 = WINDOW_VALLEY_END*OVERSAMPLING;
        else
            idx2 = idx2 + peakidx(nch) - 1;
        end
        width(nch) = max(0, (idx2-idx1)*1000/(Fs*4)); % in usec % Prevent negative - shouldnt be possible
        width(nch) = max(0, (idx2-idx1)); % in 0.1msec units
        
        % Peak-Valley Ratio - One channel at a time.
        if abs(valley(nch))<0.5,
            peak_val_ratio(nch) = peak(nch);  % To prevent blowup, set valley to 1 for ratio
        else
            peak_val_ratio(nch) = abs(peak(nch)./valley(nch));
        end
        
        if isinf(peak_val_ratio(nch)),
            peak_val_ratio(nch)=peak(nch);
        end
    end
    
    % Peak-Valley Latency (Parameter c) - All channels at once.
    peak_val_lat = (valleyidx - peakidx)*1000/(Fs*4); % in usec
    peak_val_lat = (valleyidx - peakidx); % in 0.1msec unit  
    
    % Put all additional params in
    addn_params(i,1:5*size(waves,2))=[peak_val_ratio, peak_val_lat, width, valley, valleya];
    
    % De-jitter/ Peak-align For PC Calculation
    
    % look for maximum in post-threshold window across all 4 channels
    sw = wU(((WINDOW(1))*OVERSAMPLING) : (WINDOW(2)*OVERSAMPLING), :);
    sw = sw - CM_THRESH;
    sw = sw .* (sw>0);
    
    % cm = [1:L]*sw ./ sum(sw); % by component center of mass
    cm = [1:L]*sum(sw,2) / sum(sum(sw)); % sum all components center of mass
    %center_of_mass(:,i) = cm./OVERSAMPLING;
    
    cm = round(cm);
    
    ind_start = max(cm - SNIPPET_PRE*OVERSAMPLING + WINDOW(1)*OVERSAMPLING,0);
    ind_start = min(ind_start, N*OVERSAMPLING - SNIPPET_LENGTH*OVERSAMPLING + 0);
    inds = ind_start + [1:OVERSAMPLING:SNIPPET_LENGTH*OVERSAMPLING];
    %ind_startret(i) = ind_start./OVERSAMPLING;
    new_waves(:,:,i) = wU(inds,:);
end

% Check Peak-Val Ratio Again
%  addn_params(find(isinf(addn_params(:,1))),1)=1;
%  addn_params(find(isinf(addn_params(:,2))),2)=1;
%  addn_params(find(isinf(addn_params(:,3))),3)=1;
%  addn_params(find(isinf(addn_params(:,4))),4)=1;

%toc;


if dopc==1
    
    disp(['      Doing PCs']);
    %tic
    addn_params(:,21:end) = sj_pcasvd_matspikes(new_waves,size(new_waves,2));
    %toc
end

% Calebs PCA
% K = size(waves,3); % nspikes
% ww = reshape(new_waves,[SNIPPET_LENGTH*4,K]);
% [v,d] = eig(ww*ww');
% out = zeros(K,9);
% for i = 1:9
%   out(:,i) = ww'*v(:,end - (i-1));
%   names{i} = sprintf('p%d',i);
% end


%pp = new_waves;

