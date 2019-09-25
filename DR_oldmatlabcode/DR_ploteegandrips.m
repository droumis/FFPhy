%% this isn't really working at this point
clear all; close all; clc;

% set data params
ANIMAL_DIR = '/data19/droumis/bob/bobstim_proc/';
PREFIX = 'bob';
DAY = 25;

% eeg08 = zeros(500000, 20);

% EPOCHS = [1:20];
% for EPOCHS = 1:20; % use for stairs plot
EPOCHS = 1;
% TET = [1:4 7:12];
% TET = [1 3 6 7 8 9];
TET = [9];
use_tet = TET;

eeg = loadeegstruct(ANIMAL_DIR, PREFIX, 'eeg', DAY, EPOCHS, TET);
ripple = loadeegstruct(ANIMAL_DIR, PREFIX, 'ripple', DAY, EPOCHS, TET);
theta = loadeegstruct(ANIMAL_DIR, PREFIX, 'theta', DAY, EPOCHS, TET);
ripples = loaddatastruct(ANIMAL_DIR, PREFIX, 'ripples', DAY);
multi = loaddatastruct(ANIMAL_DIR, PREFIX, 'multi', DAY);

fs = eeg{DAY(1)}{EPOCHS(1)}{TET(1)}.samprate;

% set analysis params
binsize = 3e-3;
window = 1;
backstep = 500e-3;
t = (1:round(window*fs))/fs - backstep;
t_theta = (1:round(window*150))/150 - backstep;
edges = binsize:binsize:window;

for ktet=TET
    count(ktet) = 0;
end
for kday = DAY
    for kepo = EPOCHS
        nrips = length(ripples{kday}{kepo}{use_tet}.startind);
        eeg_times = geteegtimes(eeg{kday}{kepo}{use_tet});
        theta_times = theta{kday}{kepo}{use_tet}.starttime + (1:length(theta{kday}{kepo}{use_tet}.data))/150;
        for krip = 1:nrips
            startind_rip = ripples{kday}{kepo}{use_tet}.midind(krip) - round(backstep*fs);
            endind_rip = startind_rip + round(window*fs)-1;
            starttime = ripples{kday}{kepo}{use_tet}.midtime(krip)-backstep;
            startind_eeg = find(eeg_times >starttime,1,'first');
            endind_eeg = startind_eeg + round(window*fs)-1;
            startind_theta = find(theta_times >starttime,1,'first');
            endind_theta = startind_theta + round(window*150)-1;
            if (~isempty(startind_eeg)...
                && startind_eeg > 1 ...
                && endind_eeg < length(eeg{kday}{kepo}{use_tet}.data) ...
                && endind_rip < length(ripple{kday}{kepo}{use_tet}.data)) % ...
                %&& ripples{kday}{kepo}{use_tet}.energy(krip) > 1e6)
                
                for ktet = TET    
                    count(ktet) = count(ktet) + 1;
                 
                    alligned{ktet}.eeg(count(ktet),:) = eeg{kday}{kepo}{ktet}.data(startind_eeg:endind_eeg);
                    alligned{ktet}.rip(count(ktet),:) = ripple{kday}{kepo}{ktet}.data(startind_rip:endind_rip,3);
                    alligned{ktet}.mu_rate(count(ktet),:) = histc(multi{kday}{kepo}{ktet}/1e4, edges + starttime)/binsize;
                    alligned{ktet}.theta(count(ktet),:) = theta{kday}{kepo}{ktet}.data(startind_theta:endind_theta,3);
                end
            end
        end
    end
end

%% inset lines on plot where rips detected

eegepoch = eeg{1,25}{1,EPOCHS}{1,use_tet}.data;

% dd = eegepoch(419300:422300, 1);

lines = ripples{1, 25}{1,EPOCHS}{1, use_tet}.startind;
lines = lines./1500;
figure
xinseconds = (0:1/1500:((max(size(eegepoch)))/1500)-1/1500);
xinseconds = xinseconds';
subplot(2,1,1); plot(xinseconds, eegepoch, 'Color', 'b');
yax = ylim;
for i = 1:length(lines)
patch([lines(i) lines(i) lines(i) lines(i)],[ yax(1,2) yax(1,2) yax(1,1) yax(1,1)], 'r', 'FaceAlpha', 0.8)
end
% for i = 1:length(linesEnd)
% patch([linesEnd(i) linesEnd(i) linesEnd(i)+1 linesEnd(i)],[ yax(1,2) yax(1,2) yax(1,1) yax(1,1)], 'r')
% end

% title( ['epoch' , num2str(EPOCHS) , 'tetrode' , num2str(use_tet)])
title( ['epoch' , num2str(EPOCHS) , 'Tetrode' , num2str(use_tet)])
hold on;
% % 
% 
%  set(gca,'xtick',1:20); hold on; 
% ymax = ylim;

%%

% clear all
% cd /data19/droumis/bob/bobstim_proc/EEG/
% a = load('bobeeg25-4-08');
% a = a.eeg{25}{4}{8}.data;
%% change the values of the loaded data
b = double(ripple{kday}{kepo}{use_tet}.data(:,3));

% xinseconds = (0:1/1500:((max(size(b)))/1500)-1/1500);

subplot(2,1,2); plot(xinseconds, b, 'Color', 'g'); hold on;

yax = ylim;
for i = 1:length(lines)
patch([lines(i) lines(i) lines(i) lines(i)],[ yax(1,2) yax(1,2) yax(1,1) yax(1,1)], 'r','FaceAlpha', 0.2)
end
title('rippleband')


% subplot(2,1,2); plot(xinseconds,a, 'Color', 'b'); hold on;
% for i = 1:max(size(event.excludetimes));
% line([event.excludetimes(i,1) event.excludetimes(i,1)], [-70 70], 'Color', 'r')
% end
% for i = 1:max(size(event.excludetimes));
% line([event.excludetimes(i,2) event.excludetimes(i,2)], [-50 50], 'Color', 'k')
% end
% title('unfiltered')
