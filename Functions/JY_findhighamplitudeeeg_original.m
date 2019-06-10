function [ event ] = JY_findhighamplitudeeeg( eeg,samprate,eegstart, threshold, baseline, duration )
%Find times when eeg signal exceeds a certain threshold
%   Use to filter out periods of high amplitude noise from eeg data
%   Based on code from extract ripples, uses mex code extractevents.mex
%   eeg - eeg trace
%   samprate -  eeg sampling rate
%   threshold - no. of std above mean to count as noise
%   baseline - time detection boundries for defining high amplitude periods
%   duration - (seconds) minimum duration amplitude must be above threshold


nstd=threshold;
starttime=eegstart;


% smooth the envelope:
% define the standard deviation for the Gaussian smoother which we
% apply before thresholding (this reduces sensitivity to spurious
% flucutations in the ripple envelope)

% smoothing_width = 0.002; % 4 ms
% kernel = gaussian(smoothing_width*samprate, ceil(8*smoothing_width*samprate));
% renv = smoothvect(renv, kernel);

% no smoothing used
renv=eeg;

% calculate the threshold in uV units
eegmean = mean(renv);
stdev = std(renv);
thresh = nstd;

% process data on z scored eeg, works better for some reason
% z score eeg

eegz=(eeg-eegmean)./stdev;

% find the events
% calculate the duration in terms of samples
mindur = round(duration * samprate);



renv = eeg;

% plot to check threshold and eeg
% plot(eeg,'r');
% hold on;
% plot(renv);
% 
% hold on;
% line([0 length(eeg)], [thresh thresh]);
% line([0 length(eeg)], [-thresh -thresh]);


% actual mex code to extract events
% -eegz works better than eegz
% 150 seems good for stringing together short events
tmpevent = extractevents(abs(eegz), nstd+1, baseline, 150, mindur, 0)';
% Assign the fields
% start and end indeces
event.startind = tmpevent(:,1);
event.endind = tmpevent(:,2);
% middle of energy index
event.midind = tmpevent(:,8);

% plot figure to check results
% figure;
% 
% for ii=1:size(event.startind,1)
%     line([event.startind(ii) event.endind(ii)],[0 0],'LineWidth',1000,'Color',[1 0 0]);
%     hold on;
% end
% 
% hold on;
% plot(eegz,'.b');
% hold on;
%  line([0 length(eegz)], [thresh thresh]);
%  line([0 length(eegz)], [-thresh -thresh]);

%convert the samples to times for the first three fields
event.starttime = starttime + event.startind / samprate;
event.endtime =starttime + event.endind / samprate;
event.midtime = starttime + event.midind / samprate;
event.peak = tmpevent(:,3);
event.energy = tmpevent(:,7);
event.maxthresh = (tmpevent(:,9) - baseline) / stdev;
event.excludetimes=[event.starttime event.endtime];





end

