
function [DIO] = sj_stimresp1 (animdirect,prefix,day,epoch,tet,saveg1)

%%%% PLOT EXAMPLES AND MEANS OF STIMULUS ARTIFACT REMOVAL IN EEG AND RIPPLE
%%%% TRACES. Stimulation in vhC - Stim times are in DIO file.

%%
if nargin<6,
    saveg1 = 0;
end

%% Getfiles for day


% [daypath,dayname] = fileparts(daydirect);
% currdir = pwd;
if (animdirect(end) == '/')
   animdirect = animdirect(1:end-1);
end
cd(animdirect)
% eegdir = dir('EEG');

%% DIO FIle
DIOfile = sprintf('%s/%sDIO%02d.mat', animdirect,prefix,day);

%% EEG file
EEGfile = sprintf('%s/EEG/%seeg%02d-%01d-%02d.mat', animdirect,prefix,day,epoch,tet);

%% EEG Nostim file
EEGnostimfile = sprintf('%s/EEG/%seegnostim%02d-%01d-%02d.mat', animdirect,prefix,day,epoch,tet);

%% Ripple file
RIPfile = sprintf('%s/EEG/%sripple%02d-%01d-%02d.mat', animdirect,prefix,day,epoch,tet);

%% RIPnostim file
RIPnostimfile = sprintf('%s/EEG/%sripplenostim%02d-%01d-%02d.mat', animdirect,prefix,day,epoch,tet);


%% PROCESS DATA NOW

%%DIO
load(DIOfile);
stim = DIO{day}{epoch}{15};
stim_starttime = stim.pulsetimes(:,1);
stim_endtime = stim.pulsetimes(:,2);
stim_length = stim.pulselength;
stim_isi = stim.timesincelast(2:end);


%% EEG
load(EEGfile);
e = eeg{day}{epoch}{tet}; %%% NOTE: e.samprate = 1500;
t = geteegtimes(e);
pt = stim.pulsetimes ./ 10000;
eind = lookup(pt(:,1), t);

load(EEGnostimfile);
etmp = eeg{day}{epoch}{tet}.data;

%% Ripple band
load(RIPfile);
ripamp = ripple{day}{epoch}{tet}.data(:,1);
ripenv = ripple{day}{epoch}{tet}.data(:,3);

load(RIPnostimfile);
ripampfilt = ripple{day}{epoch}{tet}.data(:,1);
ripenvfilt = ripple{day}{epoch}{tet}.data(:,3);

% ripphase6 = ripple{4}{3}{6}.data(:,2);
%etmp6 = sj_removeartifacts1(e6,stim,50);
%ripamptmp6 = sj_removeartifacts2(ripple{4}{3}{6},1,stim, 30);
%ripenvtmp6 = sj_removeartifacts2(ripple{4}{3}{6},3,stim, 30);

%%
%%%% Loop over all stimulations
%% Skip some initial and final stimulations

nstim=length(eind);
cnt=0;
for i=5:nstim-5;
    cnt=cnt+1;
    %%% NOTE: e.samprate = 1500;
    e_stim(cnt,:)=e.data((eind(i+1)-round(0.5*e.samprate)):(eind(i+1)+round(e.samprate)));              
    etmp_stim(cnt,:)=etmp(eind(i+1)-round(0.5*e.samprate):eind(i+1)+round(e.samprate));
    ripamp_stim(cnt,:)=ripamp(eind(i+1)-round(0.5*e.samprate):eind(i+1)+round(e.samprate));
    ripampfilt_stim(cnt,:)=ripampfilt(eind(i+1)-round(0.5*e.samprate):eind(i+1)+round(e.samprate));
    ripenv_stim(cnt,:)=ripenv(eind(i+1)-round(0.5*e.samprate):eind(i+1)+round(e.samprate));
    ripenvfilt_stim(cnt,:)=ripenvfilt(eind(i+1)-round(0.5*e.samprate):eind(i+1)+round(e.samprate));
    
    %multi2_stim{i} = multi2(find(multi2>=pulsest(i,1)-7500):find(multi2>=pulsest(i,1)+15000));
    %multi5_stim{i} = multi5(find(multi5>=pulsest(i,1)-7500):find(multi5>=pulsest(i,1)+15000));
    %multi6_stim{i} = multi6(find(multi6>=pulsest(i,1)-7500):find(multi6>=pulsest(i,1)+15000));
end


figure; hold on;
redimscreen_land;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%subplot(2,1,1); hold on;
plot(e_stim(15,:),'k.-','Linewidth',2,'Markersize',6);
plot(etmp_stim(15,:),'r.','Linewidth',2,'Markersize',16);
title(['Example raw and filtered EEG (tet' num2str(tet) ') around stim: Day' num2str(day) ', Epoch' num2str(epoch)]);
xlim([0 e.samprate]);
if saveg1==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['ExampleEEGStimRemove_Day' num2str(day) 'Tet' num2str(tet) '_Epoch' num2str(epoch)],'jpg');
end


figure; hold on;
redimscreen_land;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%subplot(2,1,2); hold on;
plot(mean(e_stim),'k.-','Linewidth',2,'Markersize',6);
plot(mean(etmp_stim),'r.','Linewidth',2,'Markersize',16);
title(['Mean raw and filtered EEG (tet' num2str(tet) ') around stim: Day' num2str(day) ', Epoch' num2str(epoch)]);
xlim([0 e.samprate]);
if saveg1==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['MeanEEGStimRemove_Day' num2str(day) 'Tet' num2str(tet) '_Epoch' num2str(epoch)],'jpg');
end

figure; hold on;
redimscreen_land;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%subplot(2,1,1);hold on;
plot(ripamp_stim(15,:),'k.-','Linewidth',2,'Markersize',6);
plot(ripampfilt_stim(15,:),'r.-','Linewidth',2,'Markersize',6);
title(['Example raw and filtered ripple band (tet' num2str(tet) ') around stim: Day' num2str(day) ', Epoch' num2str(epoch)]);
xlim([0 e.samprate]);
if saveg1==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['ExampleRipAmpStimRemove_Day' num2str(day) 'Tet' num2str(tet) '_Epoch' num2str(epoch)],'jpg');
end

figure; hold on;
redimscreen_land;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%subplot(2,1,2); hold on;
plot(mean(ripamp_stim),'k.-','Linewidth',2,'Markersize',6);
plot(mean(ripampfilt_stim),'r.-','Linewidth',2,'Markersize',6);
title(['Mean raw and filtered ripple band (tet' num2str(tet) ') around stim: Day' num2str(day) ', Epoch' num2str(epoch)]);
xlim([0 e.samprate]);
if saveg1==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['MeanRipEnvStimRemove_Day' num2str(day) 'Tet' num2str(tet) '_Epoch' num2str(epoch)],'jpg');
end

figure; hold on;
redimscreen_land;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%subplot(2,1,1);hold on;
plot(ripenv_stim(15,:),'k.-','Linewidth',2,'Markersize',6);
plot(ripenvfilt_stim(15,:),'r.-','Linewidth',2,'Markersize',6);
title(['Example raw and filtered ripple envelope (tet' num2str(tet) ') around stim: Day' num2str(day) ', Epoch' num2str(epoch)]);
xlim([0 e.samprate]);
if saveg1==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['ExampleRipEnvStimRemove_Day' num2str(day) 'Tet' num2str(tet) '_Epoch' num2str(epoch)],'jpg');
end

figure; hold on;
redimscreen_land;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
%subplot(2,1,2); hold on;
plot(mean(ripenv_stim),'k.-','Linewidth',2,'Markersize',6);
plot(mean(ripenvfilt_stim),'r.-','Linewidth',2,'Markersize',6);
title(['Mean raw and filtered ripple envelope (tet' num2str(tet) ') around stim: Day' num2str(day) ', Epoch' num2str(epoch)]);
xlim([0 e.samprate]);
if saveg1==1,
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['MeanRipEnvStimRemove_Day' num2str(day) 'Tet' num2str(tet) '_Epoch' num2str(epoch)],'jpg');
end


e_stim = []; etmp_stim = [];
ripamp_stim = []; ripampfilt_stim = [];
ripenv_stim = []; ripenvfilt_stim = [];










%                 %% EEg on tet 2
%                 tet = 2;
%                 load /data25/sjadhav/SJStimC_direct/EEG/RE1eeg06-3-02.mat
%                 e6 = eeg{day}{epoch}{tet};
%                 t = geteegtimes(e6);
%                 pt6 = stim.pulsetimes ./ 10000;
%                 eind6 = lookup(pt6(:,1), t);
%
%                 load /data25/sjadhav/SJStimC_direct/EEG/RE1eegnostim06-3-02.mat
%                 etmp6 = eeg{day}{epoch}{tet}.data;
%
%                 %% Ripple processed EEG on 6
%                  load /data25/sjadhav/SJStimC_direct/EEG/sjcripple06-3-02.mat
%                  ripamp6 = ripple{day}{epoch}{tet}.data(:,1);
%                  ripenv6 = ripple{day}{epoch}{tet}.data(:,3);
%                  load /data25/sjadhav/SJStimC_direct/EEG/sjcripplenostim06-3-02.mat
%                  ripampfilt6 = ripple{day}{epoch}{tet}.data(:,1);
%                  ripenvfilt6 = ripple{day}{epoch}{tet}.data(:,3);
%
%                 % ripphase6 = ripple{4}{3}{6}.data(:,2);
%
%                  %etmp6 = sj_removeartifacts1(e6,stim,50);
%                  %ripamptmp6 = sj_removeartifacts2(ripple{4}{3}{6},1,stim, 30);
%                  %ripenvtmp6 = sj_removeartifacts2(ripple{4}{3}{6},3,stim, 30);
%
%                 for i =1:length(eind6)-1;
%                      e6_stim(i,:)=e6.data(eind6(i+1)-0.5*e6.samprate:eind6(i+1)+1*e6.samprate);
%                      etmp6_stim(i,:)=etmp6(eind6(i+1)-0.5*e6.samprate:eind6(i+1)+1*e6.samprate);
%                      ripamp6_stim(i,:)=ripamp6(eind6(i+1)-0.5*e6.samprate:eind6(i+1)+1*e6.samprate);
%                      ripampfilt6_stim(i,:)=ripampfilt6(eind6(i+1)-0.5*e6.samprate:eind6(i+1)+1*e6.samprate);
%                      ripenv6_stim(i,:)=ripenv6(eind6(i+1)-0.5*e6.samprate:eind6(i+1)+1*e6.samprate);
%                      ripenvfilt6_stim(i,:)=ripenvfilt6(eind6(i+1)-0.5*e6.samprate:eind6(i+1)+1*e6.samprate);
%
%                      %multi2_stim{i} = multi2(find(multi2>=pulsest(i,1)-7500):find(multi2>=pulsest(i,1)+15000));
%                      %multi5_stim{i} = multi5(find(multi5>=pulsest(i,1)-7500):find(multi5>=pulsest(i,1)+15000));
%                      %multi6_stim{i} = multi6(find(multi6>=pulsest(i,1)-7500):find(multi6>=pulsest(i,1)+15000));
%                 end
%
%
%                 figure; hold on; subplot(2,1,1); hold on;
%                 plot(e6_stim(10,:),'k.-'); plot(etmp6_stim(10,:),'r.');
%                 title(['Example raw and filtered EEG(tet 2) around stim: Sleep2']);
%                 subplot(2,1,2); hold on;
%                 plot(mean(e6_stim),'k.-'); plot(mean(etmp6_stim),'r.');
%                 title(['Mean raw and filtered EEG (tet 2)around stim: Sleep2']);
%
%                 figure; hold on; subplot(2,1,1);hold on;
%                 plot(ripamp6_stim(10,:),'k.-'); plot(ripampfilt6_stim(10,:),'r.-');
%                 title(['Example raw and filtered ripple band (tet 2) around stim: Sleep1']);
%                 subplot(2,1,2); hold on;
%                 plot(mean(ripamp6_stim),'k.-'); plot(mean(ripampfilt6_stim),'r.-');
%                 title(['Mean raw and filtered ripple band (tet 2)around stim: Sleep1']);
%
%                 figure; hold on;
%                 subplot(2,1,1);hold on;
%                 plot(ripenv6_stim(10,:),'k.-'); plot(ripenvfilt6_stim(10,:),'r.-');
%                 title(['Example raw and filtered ripple envelope (tet 2) around stim: Sleep2']);
%                 subplot(2,1,2); hold on;
%                 plot(mean(ripenv6_stim),'k.-'); plot(mean(ripenvfilt6_stim),'r.-');
%                 title(['Mean raw and filtered ripple envelope (tet 2) around stim: Sleep2']);
%


%
%                 e6_stim = []; etmp6_stim = [];
%                 ripamp6_stim = []; ripampfilt6_stim = [];
%                 ripenv6_stim = []; ripenvfilt6_stim = [];
%



%         %% Epoch 1 - Sleep 1
%         epoch = 1;
%         stim = DIO{day}{epoch}{15};
%         stim_starttime = stim.pulsetimes(:,1);
%         stim_endtime = stim.pulsetimes(:,2);
%         stim_length = stim.pulselength;
%         stim_isi = stim.timesincelast(2:end);
%
%         %% Multiunit spikes on Tet 2, 5, 6
%         load /data25/sjadhav/SJStimC_direct/sjcmulti04.mat
%         multi2 = multi{day}{epoch}{2} ;
%         multi5 = multi{day}{epoch}{5} ;
%         multi6 = multi{day}{epoch}{6} ;
%
%
%
%
%                 %% EEg on tet 6
%                 tet = 6;
%                 load /data25/sjadhav/SJStimC_direct/EEG/sjceeg04-1-06.mat
%                 e6 = eeg{day}{epoch}{tet};
%                 t = geteegtimes(e6);
%                 pt6 = stim.pulsetimes ./ 10000;
%                 eind6 = lookup(pt6(:,1), t);
%
%                 load /data25/sjadhav/SJStimC_direct/EEG/sjceegnostim04-1-06.mat
%                 etmp6 = eeg{day}{epoch}{tet}.data;

%                 %% Ripple processed EEG on 6
%                  load /data25/sjadhav/SJStimC_direct/EEG/sjcripple04-1-06.mat
%                  ripamp6 = ripple{day}{epoch}{tet}.data(:,1);
%                  ripenv6 = ripple{day}{epoch}{tet}.data(:,3);
%                  load /data25/sjadhav/SJStimC_direct/EEG/sjcripplenostim04-1-06.mat
%                  ripampfilt6 = ripple{day}{epoch}{tet}.data(:,1);
%                  ripenvfilt6 = ripple{day}{epoch}{tet}.data(:,3);
%
%                 % ripphase6 = ripple{4}{3}{6}.data(:,2);
%
%                  %etmp6 = sj_removeartifacts1(e6,stim,50);
%                  %ripamptmp6 = sj_removeartifacts2(ripple{4}{3}{6},1,stim, 30);
%                  %ripenvtmp6 = sj_removeartifacts2(ripple{4}{3}{6},3,stim, 30);
%                 cnt = 0;
%                 for i =5:length(eind6)-10;
%                      cnt = cnt + 1;
%                      e6_stim(cnt,:)=e6.data(eind6(i+5)-0.2*e6.samprate:eind6(i+5)+0.4*e6.samprate);
%                      etmp6_stim(cnt,:)=etmp6(eind6(i+5)-0.2*e6.samprate:eind6(i+5)+0.4*e6.samprate);
%                      ripamp6_stim(cnt,:)=ripamp6(eind6(i+5)-0.2*e6.samprate:eind6(i+5)+0.4*e6.samprate);
%                      ripampfilt6_stim(cnt,:)=ripampfilt6(eind6(i+5)-0.2*e6.samprate:eind6(i+5)+0.4*e6.samprate);
%                      ripenv6_stim(cnt,:)=ripenv6(eind6(i+5)-0.2*e6.samprate:eind6(i+5)+0.4*e6.samprate);
%                      ripenvfilt6_stim(cnt,:)=ripenvfilt6(eind6(i+5)-0.2*e6.samprate:eind6(i+5)+0.4*e6.samprate);
%
%                      %multi2_stim{i} = multi2(find(multi2>=pulsest(i,1)-7500):find(multi2>=pulsest(i,1)+15000));
%                      %multi5_stim{i} = multi5(find(multi5>=pulsest(i,1)-7500):find(multi5>=pulsest(i,1)+15000));
%                      %multi6_stim{i} = multi6(find(multi6>=pulsest(i,1)-7500):find(multi6>=pulsest(i,1)+15000));
%                 end


%                 figure; hold on; subplot(2,1,1); hold on;
%                 plot(e6_stim(10,:),'k.-'); plot(etmp6_stim(10,:),'r.');
%                 title(['Example raw and filtered EEG(tet 6) around stim: Run1']);
%                 subplot(2,1,2); hold on;
%                 plot(mean(e6_stim),'k.-'); plot(mean(etmp6_stim),'r.');
%                 title(['Mean raw and filtered EEG (tet 6)around stim: Run1']);
%
%                 figure; hold on; subplot(2,1,1);hold on;
%                 plot(ripamp6_stim(10,:),'k.-'); plot(ripampfilt6_stim(10,:),'r.-');
%                 title(['Example raw and filtered ripple band (tet 6) around stim: Run1']);
%                 subplot(2,1,2); hold on;
%                 plot(mean(ripamp6_stim),'k.-'); plot(mean(ripampfilt6_stim),'r.-');
%                 title(['Mean raw and filtered ripple band (tet 6)around stim: Run1']);
%
%                 figure; hold on;
%                 subplot(2,1,1);hold on;
%                 plot(ripenv6_stim(10,:),'k.-'); plot(ripenvfilt6_stim(10,:),'r.-');
%                 title(['Example raw and filtered ripple envelope (tet 6) around stim: Run1']);
%                 subplot(2,1,2); hold on;
%                 plot(mean(ripenv6_stim),'k.-'); plot(mean(ripenvfilt6_stim),'r.-');
%                 title(['Mean raw and filtered ripple envelope (tet 6) around stim: Run1']);
%
%
%
%                 e6_stim = []; etmp6_stim = [];
%                 ripamp6_stim = []; ripampfilt6_stim = [];
%                 ripenv6_stim = []; ripenvfilt6_stim = [];





%         %% EEg on tet 2
%                 tet = 2;
%                 load /data25/sjadhav/SJStimC_direct/EEG/sjceeg04-1-02.mat
%                 e6 = eeg{day}{epoch}{tet};
%                 t = geteegtimes(e6);
%                 pt6 = stim.pulsetimes ./ 10000;
%                 eind6 = lookup(pt6(:,1), t);
%
%                 load /data25/sjadhav/SJStimC_direct/EEG/sjceegnostim04-1-02.mat
%                 etmp6 = eeg{day}{epoch}{tet}.data;
%
%                 %% Ripple processed EEG on 6
%                  load /data25/sjadhav/SJStimC_direct/EEG/sjcripple04-1-02.mat
%                  ripamp6 = ripple{day}{epoch}{tet}.data(:,1);
%                  ripenv6 = ripple{day}{epoch}{tet}.data(:,3);
%                  load /data25/sjadhav/SJStimC_direct/EEG/sjcripplenostim04-1-02.mat
%                  ripampfilt6 = ripple{day}{epoch}{tet}.data(:,1);
%                  ripenvfilt6 = ripple{day}{epoch}{tet}.data(:,3);
%
%                 % ripphase6 = ripple{4}{3}{6}.data(:,2);
%
%                  %etmp6 = sj_removeartifacts1(e6,stim,50);
%                  %ripamptmp6 = sj_removeartifacts2(ripple{4}{3}{6},1,stim, 30);
%                  %ripenvtmp6 = sj_removeartifacts2(ripple{4}{3}{6},3,stim, 30);
%                 cnt = 0;
%                 for i =5:length(eind6)-10;
%                      cnt = cnt + 1;
%                      e6_stim(cnt,:)=e6.data(eind6(i+5)-0.17*e6.samprate:eind6(i+5)+0.33*e6.samprate);
%                      etmp6_stim(cnt,:)=etmp6(eind6(i+5)-0.17*e6.samprate:eind6(i+5)+0.33*e6.samprate);
%                      ripamp6_stim(cnt,:)=ripamp6(eind6(i+5)-0.17*e6.samprate:eind6(i+5)+0.33*e6.samprate);
%                      ripampfilt6_stim(cnt,:)=ripampfilt6(eind6(i+5)-0.17*e6.samprate:eind6(i+5)+0.33*e6.samprate);
%                      ripenv6_stim(cnt,:)=ripenv6(eind6(i+5)-0.17*e6.samprate:eind6(i+5)+0.33*e6.samprate);
%                      ripenvfilt6_stim(cnt,:)=ripenvfilt6(eind6(i+5)-0.17*e6.samprate:eind6(i+5)+0.33*e6.samprate);
%
%                      %multi2_stim{i} = multi2(find(multi2>=pulsest(i,1)-7500):find(multi2>=pulsest(i,1)+15000));
%                      %multi5_stim{i} = multi5(find(multi5>=pulsest(i,1)-7500):find(multi5>=pulsest(i,1)+15000));
%                      %multi6_stim{i} = multi6(find(multi6>=pulsest(i,1)-7500):find(multi6>=pulsest(i,1)+15000));
%                 end
%
%
%                 figure; hold on; subplot(2,1,1); hold on;
%                 plot(e6_stim(10,:),'k.-'); plot(etmp6_stim(10,:),'r.');
%                 title(['Example raw and filtered EEG(tet 2) around stim: Run1']);
%                 subplot(2,1,2); hold on;
%                 plot(mean(e6_stim),'k.-'); plot(mean(etmp6_stim),'r.');
%                 title(['Mean raw and filtered EEG (tet 2) around stim: Run1']);
%
%                 figure; hold on; subplot(2,1,1);hold on;
%                 plot(ripamp6_stim(10,:),'k.-'); plot(ripampfilt6_stim(10,:),'r.-');
%                 title(['Example raw and filtered ripple band (tet 2) around stim: Run1']);
%                 subplot(2,1,2); hold on;
%                 plot(mean(ripamp6_stim),'k.-'); plot(mean(ripampfilt6_stim),'r.-');
%                 title(['Mean raw and filtered ripple band (tet 2)around stim: Run1']);
%
%                 figure; hold on;
%                 subplot(2,1,1);hold on;
%                 plot(ripenv6_stim(10,:),'k.-'); plot(ripenvfilt6_stim(10,:),'r.-');
%                 title(['Example raw and filtered ripple envelope (tet 2) around stim: Run1']);
%                 subplot(2,1,2); hold on;
%                 plot(mean(ripenv6_stim),'k.-'); plot(mean(ripenvfilt6_stim),'r.-');
%                 title(['Mean raw and filtered ripple envelope (tet 2) around stim: Run1']);
%
%
