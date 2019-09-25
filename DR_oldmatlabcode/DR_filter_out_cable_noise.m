%% this isn't really working at this point



% clear all
cd /data19/droumis/bob/bobstim_proc/EEG/
a = load('bobeeg25-4-08');
a = a.eeg{25}{4}{8}.data;

b = load('bobripple25-4-08');
b = b.ripple{25}{4}{08}.data;
c = double(b(:,3));


% events = JY_findhighamplitudeeeg(a, 1500, 0, .1, 1, 0.1);

% events = JY_findhighamplitudeeeg(c, 1500, 0, 1, 0, .5);
hold off
clear events
close all
%**************************
% event = JY_findhighamplitudeeeg(c, 1500, 0, .1, .1); % work on this
% ***********************

%function [ event ] = JY_findhighamplitudeeeg( eeg,samprate,eegstart, threshold, duration )
%[ event ] = JY_findhighamplitudeeeg( eeg,samprate,eegstart, threshold,
%baseline, duration ) old version

%   Use to filter out periods of high amplitude noise from eeg data
%   Based on code from extract ripples, uses mex code extractevents.mex
%   eeg - eeg trace
%   samprate -  eeg sampling rate
%   threshold - no. of std above mean to count as noise
%   baseline - time detection boundries for defining high amplitude periods
%   duration - (seconds) minimum duration amplitude must be above threshold

%yuck
% clear excludestartind
% excludestartind(:,1) = event.excludetimes(:,1);
% excludeendind(:,1) = event.excludetimes(:,2);

%use for a
% xinseconds = (0:1/1500:((max(size(a)))/1500)-1/1500);
% subplot(2,1,1); plot(xinseconds,a, 'Color', 'b'); hold on;

%use for c
xinseconds = (0:1/1500:((max(size(c)))/1500)-1/1500);
subplot(2,1,1); plot(xinseconds,c, 'Color', 'g'); hold on;

for i = 1:max(size(event.excludetimes));
line([event.excludetimes(i,1) event.excludetimes(i,1)], [-70 70], 'Color', 'r')
end
for i = 1:max(size(event.excludetimes));
line([event.excludetimes(i,2) event.excludetimes(i,2)], [-50 50], 'Color', 'k')
end
title('rippleband')
hold on;

subplot(2,1,2); plot(xinseconds,a, 'Color', 'b'); hold on;
for i = 1:max(size(event.excludetimes));
line([event.excludetimes(i,1) event.excludetimes(i,1)], [-70 70], 'Color', 'r')
end
for i = 1:max(size(event.excludetimes));
line([event.excludetimes(i,2) event.excludetimes(i,2)], [-50 50], 'Color', 'k')
end
title('unfiltered')


% %% make new eeg without extracttimes
% b = NaN(max(size(a)),1);
%     for rewrite = 1:length(b);
%         if b(rewrite,1) > 
%            
%         b(rewrite,1) = a(rewrite,1);
%     end
%     refaddedchan = newchan + ref;
%     figure
%     plot(newchan,1); hold on;
% else %equal length
%     refaddedchan = curchan + ref;
%     figure
% 
% end
% plot(refaddedchan, 'Color', 'r');
