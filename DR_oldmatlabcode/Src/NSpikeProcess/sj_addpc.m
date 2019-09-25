function [addn_params, new_waves] = sj_addpc(waves, triggers, dopc)

% Shantanu - May 2012
% From sj_addparams. Getting rid of additional parameters, and only keeping PCs
% Usually called from sj_nmakeparams_pc.m

%[addn_params] = sj_addpc (waves, triggers,dopc)

% ------------- Description from sj_addparams.m------------------
% 20 May 2011. Shantanu
% Caleb's Peak ALign by C-M. Using Over-sampled waveform to get additional
% parameters. Then Peak Aligning and getting PC if asked for.
% function pp = peak_align(waves)

if nargin<2 || isempty(triggers),
    triggers = 40*ones(size(waves,2),1);
end
if nargin<3
    dopc = 1;
end

% This code for PCs only. If 0, returns empty matrix for addn_params 
 
    
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

if dopc==1
    addn_params = nan(size(waves,3),(3)*size(waves,2)); % 3 PCs for each channel
    %center_of_mass = nan(1,size(waves,3));
    %ind_startret = nan(1,size(waves,3));
end

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
%toc;

if dopc==1
    disp(['      Doing PCs']);
    %tic
    addn_params = sj_pcasvd_matspikes(new_waves,size(new_waves,2));
    %toc
    
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
    
else
    
    addn_params = [];
    
    
    
end




