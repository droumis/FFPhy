function [heights, widths, peak_locs, mask] = thresholded_peaks_fortet(spikes, swaveforms)

%% SHANTANU - OCT 11, 2007.  SAME AS THRESHOLDED_PEAKS EXCEPT CALLED WITH 
%% SWAVEFORMS, NO ERROR IN ISTHREHSLO

% THRESHOLDED_PEAKS  Finds height/width of threshold crossing peaks.
%     [HEIGHTS, WIDTHS] = THRESHOLDED_PEAKS(SPIKES) takes a spike-sorting
%     SS object SPIKES with M spikes and returns two (M x 1) vectors
%     HEIGHTS and WIDTHS, containing measurements for each waveform in
%     the SPIKES object.
%
%     HEIGHTS and WIDTHS are determined from the central peak and 
%     require both spikes.threshT and spikes.threshV(1:2) to be defined
%     (waveforms from low and high thresholds are handled automatically).
%     The central peak is taken as the contiguous region starting at the
%     spikes.threshT sample and ending when threshold (either low or
%     high) is recrossed in the opposite direction.  Peak height
%     corresponds to the waveform value in this region that is farthest
%     from threshold and peak width is length of the region in samples.
%
%     [HEIGHTS, WIDTHS, PEAK_LOCS] = THRESHOLDED_PEAKS(SPIKES) returns
%     an (M x 3) matrix PEAK_LOCS in which the second column gives the
%     location (in samples) in the waveform that attained the HEIGHT, 
%     while the first and third columns give the threshold cross and
%     recross locations (in samples) used to compute the WIDTHS.
%
%     [HEIGHTS, WIDTHS, PEAK_LOCS, MASK] = THRESHOLDED_PEAKS(SPIKES) further
%     returns an (M x N) matrix MASK, where (M x N) are the dimensions of
%     the spikes.waveforms matrix.  MASK is a binary matrix with a 
%     contiguous block of 1's in each row marking the samples contained in
%     the thresholded central peak for that row.

%   Last Modified By: sbm on Thu Oct  6 20:30:16 2005
%                     Shantanu Oct 11, 2007: Change for my tetrode data    

%%%%%%%%%% ARGUMENT CHECKING
if (~isfield(spikes, 'threshT'))
	error('SS:threshT_Undefined', 'The SS object must define the sample index of threshold crossing.');
elseif (~isfield(spikes, 'threshV') | (size(spikes.threshV,2) ~= 2))
	error('SS:threshV_Error', 'The SS object must define both low and high thresholds.  Use +/- Inf if only one of the two was used.'); 
elseif ((spikes.threshV(1,2) - spikes.threshV(1,1)) < 0)
    error('SS:threshV_Illegal', 'The SS object high threshold must be greater than its low threshold.');
end

%%%%%%%%%% FIND PEAK REGIONS
% This could be done with some straightforward code in a big, ugly 'for' 
% loop, but, hey, what's Matlab for if not replacing FOR loops with big, 
% ugly, vectorized code?  At least its faster this way.

% Mark each time a waveform crosses threshold by taking the derivative
% of the binary thresholded spike data.  Zero-padding ensures that the
% first and last sample points can be marked, with the result that 
% the # of times the signal crosses threshold (marked with +1) equals
% the number of times it returns below it (marked with -1).

numspikes = size(swaveforms, 1);
numsamples = size(swaveforms, 2);
%Checking validity - for +ve thresholds for now - can incorporate -ve also
valid_ch = find(spikes.threshV(:,2)<1000); 
% Change thresholds appropriately
lothrs = sum(spikes.threshV(valid_ch,1));
hithrs = sum(spikes.threshV(valid_ch,2));


%figure; hold on; plot(swaveforms(1:1000,:)','b.'); plot(spikes.threshT+1*ones(size([-800:800])),[-800:800],'k-'); plot(0:40, hithrs*ones(size([0:40])),'r-');
%figure; hold on; plot(spikes.waveforms_ch2(1:1000,:)','b.'); plot(spikes.threshT+1*ones(size([-800:800])),[-800:800],'k-'); plot(0:40, 40*ones(size([0:40])),'r-');



lopeaks = diff(padmatrix([swaveforms <= lothrs], [1 1 0 0]), 1, 2); 
hipeaks = diff(padmatrix([swaveforms >= hithrs], [1 1 0 0]), 1, 2);

% The lo/hi markers can be combined; we'll remember which is which.
if ~any(spikes.threshV(:,1) ~= -Inf)
    isThreshLo = false(numspikes,1);
else
    isThreshLo = (swaveforms(:, spikes.threshT+2) <= lothrs);
end

allpeaks = zeros(size(lopeaks));
allpeaks( isThreshLo,:) = lopeaks( isThreshLo,:);
allpeaks(~isThreshLo,:) = hipeaks(~isThreshLo,:);

clear lopeaks hipeaks;     % manual garbage collection . . . 

% Read the marks into a list of peak start/stop indices.  Note that the
% +1's line up properly but the -1's correspond to the first sample _after_
% threshold recrossing; so we subtract 1 to get start/stop indices for the
% contiguous block of supra-threshold samples.
[start_r start_c] = find(allpeaks == 1);     % pts that first exceed threshold
[finis_r finis_c] = find(allpeaks == -1);    % pts that just sank below
start = sortrows([start_r start_c]);         % sort both starts and ...
finis = sortrows([finis_r finis_c]);         %  ... ends, so we can  ...
peak_inds = [start (finis(:,2)-1)];         %  ... combine them (since they're 1:1).
clear allpeaks;

% Next, select out the peaks that contain spikes.threshT.
%peaks_center = find((peak_inds(:,2) <= spikes.threshT+1) & (peak_inds(:,3) >= spikes.threshT));
peaks_center = find((peak_inds(:,2) <= spikes.threshT+1) & (peak_inds(:,3) >= spikes.threshT));
%peaks_center = find(peak_inds(:,3) >= spikes.threshT);

if (length(peaks_center) ~= numspikes)
    disp(['Warning! Some waveforms (' ...
          num2str(size(swaveforms, 1)-length(peaks_center)) ') do not have a central peak with these parameters.' ...
          ' Can happen with tets, but check carefully']);
                              
end

% Make the central peak mask; we need this to figure out the heights.

% %%% Shantanu - Oct11, 2007
% peaks_center=unique(peaks_center);
% idxwaves=peaks_center;  nwaves=length(peaks_center);

mask = zeros(size(swaveforms));
for k = 1:length(peaks_center)
    peakrow = peak_inds(peaks_center(k),:);
    mask(peakrow(1), peakrow(2):peakrow(3)) = 1; 
end

% OK, now heights are fairly easy.  Just need to keep lo/high thresholds straight.
peaks = mask .* swaveforms;
heights = zeros(size(swaveforms, 1), 1);
[loHeights, loWhere]= min(peaks, [], 2);
[hiHeights, hiWhere] = max(peaks, [], 2);
heights( isThreshLo) = loHeights( isThreshLo);
heights(~isThreshLo) = hiHeights(~isThreshLo);

% And widths are even easier.
widths = -1*ones(numspikes,1);
widths(peaks_center) = peak_inds(peaks_center, 3) - peak_inds(peaks_center, 2);

% Finally, we build the peak_locs by reorganizing information we already have.
if (nargout > 2)
    peak_locs = zeros(size(swaveforms, 1), 3);
    peak_locs(peaks_center, [1,3]) = peak_inds(peaks_center, [2,3]); % width info
    peak_locs( isThreshLo,2) = loWhere( isThreshLo);      % height info
    peak_locs(~isThreshLo,2) = hiWhere(~isThreshLo);
end

% Thats all, folks!
