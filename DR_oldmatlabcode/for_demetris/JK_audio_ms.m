%% clear all
clear; 

% set data params
ANIMAL_DIR = '/data19/droumis/bob/bobstim_proc/';
PREFIX = 'bob';
DAY = 25;
% EPOCHS = [2 4 6];
EPOCHS = [1];
%TET = [1:4 7:12];
TET = [1 3 6 7 8 9];
dio_channels = 1;


% load data
eeg = loadeegstruct(ANIMAL_DIR, PREFIX, 'eeg', DAY, EPOCHS, TET);
theta = loadeegstruct(ANIMAL_DIR, PREFIX, 'theta', DAY, EPOCHS, TET);
multi = loaddatastruct(ANIMAL_DIR, PREFIX, 'multi', DAY);
dio = loaddatastruct(ANIMAL_DIR, PREFIX, 'DIO', DAY);

fs = eeg{DAY(1)}{EPOCHS(1)}{TET(1)}.samprate;

% set analysis params
binsize = 20e-3;
window = 4;
backstep = 2;
t = (1:round(window*fs))/fs - backstep;
t_theta = (1:round(window*150))/150 - backstep;
edges = binsize:binsize:window;

for ktet=TET
    count(ktet) = 0;
end
for kday = DAY
    for kepo = EPOCHS
        for kstim = dio_channels
            ntrials = length(dio{kday}{kepo}{kstim}.pulsetimes);
            eeg_times = geteegtimes(eeg{kday}{kepo}{1});
            for ktrial = 1:500
                starttime =dio{kday}{kepo}{kstim}.pulsetimes(ktrial,1)/10000-backstep;
                starttime = starttime; % + rand - .5;
                startind_eeg = find(eeg_times <starttime,1,'last');
                endind_eeg = startind_eeg + round(window*fs)-1;

                if (~isempty(startind_eeg)...
                        && startind_eeg > 1 ...
                        && endind_eeg < length(eeg{kday}{kepo}{1}.data))
                    % ...
                    
                    
                    for ktet = TET
                        count(ktet) = count(ktet) + 1;
                        alligned{ktet}.stim{kstim}.eeg(count(ktet),:) = eeg{kday}{kepo}{ktet}.data(startind_eeg:endind_eeg);
                        alligned{ktet}.stim{kstim}.mu_rate(count(ktet),:) = histc(multi{kday}{kepo}{ktet}/1e4, edges + starttime)/binsize;
                        
                    end
                end
            end
        end
    end
end



%%
close all
for kstim = dio_channels
    for ktet = TET
        fig = figure('color','w');
        subplot(2,1,1)
        hold on
        good_idxs = max(abs(alligned{ktet}.stim{kstim}.eeg),[],2)<3000;
        %plot(t,alligned{ktet}.eeg(1:300,:)','g','linewidth',1)
        plot(t,mean(alligned{ktet}.stim{kstim}.eeg(good_idxs,:),1),'k','linewidth',2)
        plot(t,mean(alligned{ktet}.stim{kstim}.eeg(good_idxs,:),1)+std(alligned{ktet}.stim{kstim}.eeg(good_idxs,:),[],1)/sqrt(sum(good_idxs)),'r','linewidth',2)
        plot(t,mean(alligned{ktet}.stim{kstim}.eeg(good_idxs,:),1)-std(alligned{ktet}.stim{kstim}.eeg(good_idxs,:),[],1)/sqrt(sum(good_idxs)),'r','linewidth',2)
        %ylim([-200 200])
        xlim([t(1) t(end)])
        ylabel('EEG')
        title(['Tetrode ' num2str(ktet)])
        
        
        subplot(2,1,2)
        hold on
        plot(edges-backstep, mean(alligned{ktet}.stim{kstim}.mu_rate,1),'k','linewidth',2)
        plot(edges-backstep, mean(alligned{ktet}.stim{kstim}.mu_rate(good_idxs,:),1)+std(double(alligned{ktet}.stim{kstim}.mu_rate(good_idxs,:)),[],1)/sqrt(sum(good_idxs)),'r','linewidth',2)
        plot(edges-backstep, mean(alligned{ktet}.stim{kstim}.mu_rate(good_idxs,:),1)-std(double(alligned{ktet}.stim{kstim}.mu_rate(good_idxs,:)),[],1)/sqrt(sum(good_idxs)),'r','linewidth',2)
        ylabel('MU Rate (Hz)')
        
        
        
        %
        % subplot(4,1,4)
        % plot(t_theta,mean(alligned{ktet}.theta(good_idxs,:),1),'k','linewidth',2)
        % plot(t_theta,mean(alligned{ktet}.theta(good_idxs,:),1)+std(double(alligned{ktet}.theta(good_idxs,:)),[],1)/sqrt(sum(good_idxs)),'r','linewidth',2)
        % plot(t_theta,mean(alligned{ktet}.theta(good_idxs,:),1)-std(double(alligned{ktet}.theta(good_idxs,:)),[],1)/sqrt(sum(good_idxs)),'r','linewidth',2)
        % xlim([t(1) t(end)])
        %
        
        xlabel('time (s)')
        saveas(gcf, ['audiostim_' num2str(kstim) '_tetrode_' num2str(ktet) '.fig'])
    end
end

%%
