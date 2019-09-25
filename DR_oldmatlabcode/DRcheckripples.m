%% clear all

% set data params
ANIMAL_DIR = '/data19/droumis/bob/bobstim_proc/';
PREFIX = 'bob';
DAY = 25;

% eeg08 = zeros(500000, 20);

% EPOCHS = [1:20];
% for EPOCHS = 1:20; % use for stairs plot
EPOCHS = 3;
% TET = [1:4 7:12];
% TET = [1 3 6 7 8 9];
TET = [9];
use_tet = TET;


% load data... right now i just need the eeg and ripples
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

% %% calculate cross correlations
% hippo = [8 9];
% ms = [1 3 6 7];
% 
% for khippo = hippo
%     cnt = 0;
%     for kms = ms
%         cnt = cnt + 1;
%         [site_corr{khippo}.corr_eeg(kms,:), eeglaggs] = xcov(mean(alligned{khippo}.eeg,1),mean(alligned{kms}.eeg,1));
%         [mx mxidx] = max(site_corr{khippo}.corr_eeg(kms,:));
%         site_corr{khippo}.max_eeg(cnt) = eeglaggs(mxidx)/1500;
%         
%         [site_corr{khippo}.corr_rip(kms,:), riplaggs] = xcov(mean(alligned{khippo}.rip,1),mean(alligned{kms}.rip,1));
%         [mx mxidx] = max(site_corr{khippo}.corr_rip(kms,:));
%         site_corr{khippo}.max_rip(cnt) = riplaggs(mxidx)/1500;
%         
%         
%         [site_corr{khippo}.corr_mu(kms,:), mulaggs] = xcov(mean(alligned{khippo}.mu_rate,1),mean(alligned{kms}.mu_rate,1));
%         [mx mxidx] = max(site_corr{khippo}.corr_mu(kms,:));
%         site_corr{khippo}.max_mu(cnt) = mulaggs(mxidx)*binsize;
%     end
% end
% 
% 
%%
% close all
% for ktet = TET
% fig = figure('color','w');
% subplot(3,1,1)
% hold on
% good_idxs = max(abs(alligned{ktet}.eeg),[],2)<3000;
% %plot(t,alligned{ktet}.eeg(1:300,:)','g','linewidth',1)
% plot(t,mean(alligned{ktet}.eeg(good_idxs,:),1),'k','linewidth',2)
% plot(t,mean(alligned{ktet}.eeg(good_idxs,:),1)+std(alligned{ktet}.eeg(good_idxs,:),[],1)/sqrt(sum(good_idxs)),'r','linewidth',2)
% plot(t,mean(alligned{ktet}.eeg(good_idxs,:),1)-std(alligned{ktet}.eeg(good_idxs,:),[],1)/sqrt(sum(good_idxs)),'r','linewidth',2)
% %ylim([-200 200])
% xlim([t(1) t(end)])
% ylabel('EEG')
% title(['Tetrode ' num2str(ktet)])
% 
% subplot(3,1,2)
% hold on
% %plot(t,alligned{ktet}.rip','g','linewidth',1)
% plot(t,mean(alligned{ktet}.rip(good_idxs,:),1),'k','linewidth',2)
% plot(t,mean(alligned{ktet}.rip(good_idxs,:),1)+std(double(alligned{ktet}.rip(good_idxs,:)),[],1)/sqrt(sum(good_idxs)),'r','linewidth',2)
% plot(t,mean(alligned{ktet}.rip(good_idxs,:),1)-std(double(alligned{ktet}.rip(good_idxs,:)),[],1)/sqrt(sum(good_idxs)),'r','linewidth',2)
% xlim([t(1) t(end)])
% ylabel('Ripple Band Envelope')
% 
% subplot(3,1,3)
% hold on
% plot(edges-backstep, mean(alligned{ktet}.mu_rate,1),'k','linewidth',2)
% plot(edges-backstep, mean(alligned{ktet}.mu_rate(good_idxs,:),1)+std(double(alligned{ktet}.mu_rate(good_idxs,:)),[],1)/sqrt(sum(good_idxs)),'r','linewidth',2)
% plot(edges-backstep, mean(alligned{ktet}.mu_rate(good_idxs,:),1)-std(double(alligned{ktet}.mu_rate(good_idxs,:)),[],1)/sqrt(sum(good_idxs)),'r','linewidth',2)
% ylabel('MU Rate (Hz)')
% 
% 
% 
% % 
% % subplot(4,1,4)
% % plot(t_theta,mean(alligned{ktet}.theta(good_idxs,:),1),'k','linewidth',2)
% % plot(t_theta,mean(alligned{ktet}.theta(good_idxs,:),1)+std(double(alligned{ktet}.theta(good_idxs,:)),[],1)/sqrt(sum(good_idxs)),'r','linewidth',2)
% % plot(t_theta,mean(alligned{ktet}.theta(good_idxs,:),1)-std(double(alligned{ktet}.theta(good_idxs,:)),[],1)/sqrt(sum(good_idxs)),'r','linewidth',2)
% % xlim([t(1) t(end)])
% % 
% 
% xlabel('time (s)')
% saveas(gcf, ['tetrode_' num2str(ktet) 'alligned_on' num2str(use_tet) '.fig'])
%  end
% %%
% for khippo=hippo
% fig = figure('color','w');
% subplot(3,1,1)
% plot(eeglaggs/1500,site_corr{khippo}.corr_eeg')
% xlim([-.3 .3])
% subplot(3,1,2)
% plot(riplaggs/1500,site_corr{khippo}.corr_rip')
% xlim([-.3 .3])
% subplot(3,1,3)
% plot(mulaggs/150,site_corr{khippo}.corr_mu')
% xlim([-.3 .3])
% saveas(gcf, ['hippo_' num2str(khippo) '_xcorr.fig'])
% end
% 
% fig = figure('color','w');
% for khippo = hippo
% subplot(1,2,1)
% hold on
% plot(site_corr{khippo}.max_rip*1000, khippo*ones(size(site_corr{khippo}.max_rip)), '.k')
% xlim([-10 10]);
% subplot(1,2,2)
% hold on
% plot(site_corr{khippo}.max_mu*1000, khippo*ones(size(site_corr{khippo}.max_mu)), '.k')
% xlim([-10 10]); 
% end

%% use this for ploting the stairs across epochs for a single tetrode
numripstoplot([EPOCHS]) = nrips;

%% generate variable to hold EEG for entire tetrode across epochs. yuck
% num2str(TET)+num2str(EPOCHS)
% 
% num2str(TET)+num2str(EPOCHS) = eeg{1,25}{1,EPOCHS}{1,TET}.data;
%% Use these for ploting the stairs

% 
% end
% 
%  figure 
%  stairs(numripstoplot,'LineWidth',4,'Color', 'm')
%  title 'tetrode 9 ripples over 20 epochs, even numbered are stim'
%  set(gca,'xtick',1:20); hold on; 
% ymax = ylim;
%  for i = 1:2:20
% patch([i i+1 i+1 i], [ymax(1,2) ymax(1,2)  ymax(1,1)  ymax(1,1)], 'c','FaceAlpha', 0.2);
% end
%  set(gca,'XTickLabel',{'2_P', 'rest', '8_P', 'rest', '20_P', 'rest', '60_P', 'rest', '100_P', 'rest', '2_C', 'rest', '8_C', 'rest', '20_C', 'rest', '60_C', 'rest', '100_C', 'rest'})
% 


%% inset lines on plot where rips detected

eegepoch = eeg{1,25}{1,EPOCHS}{1,use_tet}.data;

% dd = eegepoch(419300:422300, 1);

lines = ripples{1, 25}{1,EPOCHS}{1, use_tet}.startind;
figure

plot(eegepoch, 'Color', 'b');
for i = 1:length(lines)
line([lines(i) lines(i)],[-500 500], 'Color', 'r')
end
% title( ['epoch' , num2str(EPOCHS) , 'tetrode' , num2str(use_tet)])
title( ['epoch' , num2str(EPOCHS) , 'Tetrode' , num2str(use_tet)])
hold on
% 
 %% load other channel's rip to overlay
%vars replaced use_tet8
% set data params
ANIMAL_DIR = '/data19/droumis/bob/bobstim_proc/';
PREFIX = 'bob';
% DAY = 25;
%EPOCHS = [2 4 6];
%TET = [1:4 7:12];
TET = [1 3 6 7 8 9];
use_tet8 = 8;


% load data
eeg = loadeegstruct(ANIMAL_DIR, PREFIX, 'eeg', DAY, EPOCHS, TET);
ripple = loadeegstruct(ANIMAL_DIR, PREFIX, 'ripple', DAY, EPOCHS, TET);
theta = loadeegstruct(ANIMAL_DIR, PREFIX, 'theta', DAY, EPOCHS, TET);
ripples = loaddatastruct(ANIMAL_DIR, PREFIX, 'ripples', DAY);
multi = loaddatastruct(ANIMAL_DIR, PREFIX, 'multi', DAY);



for ktet=TET
    count(ktet) = 0;
end
for kday = DAY
    for kepo = EPOCHS
        nrips = length(ripples{kday}{kepo}{use_tet8}.startind);
        eeg_times = geteegtimes(eeg{kday}{kepo}{use_tet8});
        theta_times = theta{kday}{kepo}{use_tet8}.starttime + (1:length(theta{kday}{kepo}{use_tet8}.data))/150;
        for krip = 1:nrips
            startind_rip = ripples{kday}{kepo}{use_tet8}.midind(krip) - round(backstep*fs);
            endind_rip = startind_rip + round(window*fs)-1;
            starttime = ripples{kday}{kepo}{use_tet8}.midtime(krip)-backstep;
            startind_eeg = find(eeg_times >starttime,1,'first');
            endind_eeg = startind_eeg + round(window*fs)-1;
            startind_theta = find(theta_times >starttime,1,'first');
            endind_theta = startind_theta + round(window*150)-1;
            if (~isempty(startind_eeg)...
                && startind_eeg > 1 ...
                && endind_eeg < length(eeg{kday}{kepo}{use_tet8}.data) ...
                && endind_rip < length(ripple{kday}{kepo}{use_tet8}.data)) % ...
                %&& ripples{kday}{kepo}{use_tet8}.energy(krip) > 1e6)
                
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



%% check rips overlay of tet8

% eegepoch = eeg{1,25}{1,EPOCHS}{1,use_tet8}.data;
lines = ripples{1, 25}{1,EPOCHS}{1, use_tet8}.startind;
% plot(eegepoch, 'Color', 'b')
for i = 1:length(lines)
line([lines(i) lines(i)],[-450 450], 'Color', 'y')
end
% end


