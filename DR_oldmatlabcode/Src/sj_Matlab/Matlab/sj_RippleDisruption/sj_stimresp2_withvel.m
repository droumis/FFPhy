
function [ep] = sj_stimresp2_withvel (animdirect,prefix,days,allepochs,tet,ref, figopt1,saveg1,figopt2,saveg2,savedata)

%%%% PLOT and compare Response to Probe Stimulation in different Sleep
%%%% Epochs across days, controlling for velocity of animal
% eg
% sj_stimresp2_withvel('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',1:9,[1 3 5],5,3,0,0,1,0,0);
% sj_stimresp2_withvel_longer('/data25/sjadhav/RippleInterruption/REd_direct','REd',5,[1 3 5],3,1,1,0,0,0,0);
% sj_stimresp2_withvel_longer('/data25/sjadhav/RippleInterruption/REd_direct','REd',5,[1 3 5],3,[],1,0,0,0,0);
% sj_stimresp2_withvel('/data25/sjadhav/RippleInterruption/REf_direct','REe',4,[1 3 5],11,9,1,0,0,0,0);
% sj_stimresp2_withvel('/data25/sjadhav/RippleInterruption/SJStimC_direct','sjc',3:7,[1 3],6,3,0,0,1,0,0);

% figopt1: Make plots for individual days
% saveg1: Save individual day graphs
% figopt2: Make final plot averaging across days
% saveg2: Save final graphs
% savedata: Save mat file with data

%%
if nargin<6 || isempty(ref),
    display('This is without Reference');
    
end
if nargin<7,
    figopt1 = 0;
end
if nargin<8,
    saveg1 = 0;
end
if nargin<9,
    figopt2 = 0;
end
if nargin<10,
    saveg2 = 0;
end
if nargin<11,
    savedata = 0;
end

savedir = '/data25/sjadhav/RippleInterruption/ProcessedData/EPSP';


% %directoryname = '/data25/sjadhav/SJStimC_direct';
% %prefix = 'sjc';
% directoryname = '/data25/sjadhav/RE1_direct';
% prefix = 'RE1';
% days = [1:9];
% tet=6;

%% Fixed parameters
Fs=30; %video sampling rate
tsamp=1/30; % in sec
velfiltlth = Fs/4; % Filter for smoothing velocity: 1 sec=Fs points, Std Dev = Lth/4

% Variable parameters
thrsvel=50;  %<x cm per sec is still
pastsecs = 1;

pre_frac=0.2; %eg. 0.2*e.samprate
post_frac=0.4; %eg. 0.4*e.samprate

indivbar=1;
indepsps=1;


e.samprate=1500;
startidx = round(pre_frac*e.samprate)+7;    %For 100 Hz sampling rate, 5ms = 7.5 points. Stim artifact over by 3 ms
endidx = round(pre_frac*e.samprate)+75;      % Peak generally at 10-15ms, EPSP over by 40-50 ms, 50 ms = 75 points
% WIndow for EPSP Ampl
% OR, better way - go by time
stimtime=pre_frac*e.samprate*(1000/e.samprate);
resp_win = stimtime + [7, 50]; % ms
resp_win_idx = floor(resp_win*e.samprate/1000);

base_win = stimtime - [105, 5];
base_win_idx = floor(base_win*e.samprate/1000);

startidx = resp_win_idx(1);
endidx = resp_win_idx(2);



%% Loop over days and load files
directoryname = animdirect;
if (animdirect(end) == '/')
    animdirect = animdirect(1:end-1);
end
cd(animdirect);
clr = {'b','r','g','c','m','y','k','b','r','g'};

e_pre = []; e_post = []; e_post2=[];
pre=[]; post1=[]; post2=[];

for d = 1:length(days)
    
    day = days(d);
    %% Day n
    DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
    load(DIOfile);
    
    %% Pos file
    posfile = sprintf('%s/%spos%02d.mat', directoryname, prefix, day);
    load(posfile);
    
    clr = {'k','r','g','c','b','m'};
    
    if figopt1==1
        figure(1); hold on;
        %redimscreen;
        redimscreen_land;
        orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','bold');
        set(0,'defaultaxeslinewidth',2);
    end
    
    %% Epoch 1,3,5... - Sleep 1,2,3...
    %%% Loop over epochs
    for ep = 1:length(allepochs)
        
        e_stim = []; etmp_stim = [];
        epoch = allepochs(ep);
        
        stim = DIO{day}{epoch}{15};
        if isempty(stim)
            stim = DIO{day}{epoch}{16};
        end
        
        stim_starttime = stim.pulsetimes(:,1);
        stim_endtime = stim.pulsetimes(:,2);
        stim_length = stim.pulselength;
        stim_isi = stim.timesincelast(2:end);
        
        
        
        %% EEg file for given tet
        
        EEGfile = sprintf('%s/EEG/%seeg%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tet);
        load(EEGfile);
        e = eeg{day}{epoch}{tet};
        t = geteegtimes(e);
        pt = stim.pulsetimes ./ 10000;
        eind = lookup(pt(:,1), t);
        
        % EEg File for Reference
        if ~isempty(ref)
            EEGfile_ref = sprintf('%s/EEG/%seeg%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,ref);
            load(EEGfile_ref);
            e_ref = eeg{day}{epoch}{ref};
        end
        
        %         EEGnostimfile = sprintf('%s/EEG/%seegnostim%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tet);
        %         load(EEGnostimfile);
        %         etmp = eeg{day}{epoch}{tet}.data;
        
        
        
        %%% Save response to stimulation
        %         for i =1:length(eind)-2;
        %             e_tempstim(i,:)=e.data(eind(i+1)-round(0.2*e.samprate):eind(i+1)+round(0.4*e.samprate));
        %             if ~isempty(ref)
        %                 e_stim_ref(i,:)=e_ref.data(eind(i+1)-round(0.2*e.samprate):eind(i+1)+round(0.4*e.samprate));
        %                 e_stim(i,:)=e_tempstim(i,:)+e_stim_ref(i,:);
        %             else
        %                 e_stim(i,:)=e_tempstim(i,:);
        %             end
        % %             etmp_stim(i,:)=etmp(eind(i+1)-round(0.2*e.samprate):eind(i+1)+round(0.4*e.samprate));
        %         end
        for i =1:length(eind)-2;
            e_tempstim(i,:)=e.data(eind(i+1)-round(pre_frac*e.samprate):eind(i+1)+round(post_frac*e.samprate));
            if ~isempty(ref)
                e_stim_ref(i,:)=e_ref.data(eind(i+1)-round(pre_frac*e.samprate):eind(i+1)+round(post_frac*e.samprate));
                e_stim(i,:)=e_tempstim(i,:)+e_stim_ref(i,:);
            else
                e_stim(i,:)=e_tempstim(i,:);
            end
            %             etmp_stim(i,:)=etmp(eind(i+1)-round(0.2*e.samprate):eind(i+1)+round(0.4*e.samprate));
        end
        taxis = [1:size(e_stim,2)]*1000/e.samprate;
        taxis = taxis - pre_frac*1000;
        
        
        %%% Get Velocity - DONT DO
        %         currvel = pos{day}{epoch}.data(:,5);
        %         %remove_ind_sec = find(currvel < 0);
        %         %tempvel = currvel(1:Fs*floor(length(currvel)/Fs));
        %         %x = mean(reshape(tempvel,Fs,floor(length(tempvel)/Fs)),1);  %% Get velocity estimate per second
        %         veltime =  pos{day}{epoch}.data(:,1);
        %         velind = lookup(pt(:,1), veltime);
        %         velindk = velind(2:end-1);
        %         %stim_veltime = veltime(velindk);
        %         %eindk = eind(2:end-1);
        %         %stim_eegtime = t(eindk)';
        %
        %         %% Look at instantaneous smoothed velocity
        %         % Smooth velocity over a second
        %         filt = gaussian(velfiltlth/4,velfiltlth);
        %         smooth_currvel = smoothvect(currvel,filt);
        %
        %         % look at instantaneous velocity at stimulation times
        %         stim_vel = smooth_currvel(velindk);
        %         % Divide into moving and still
        %         moving = find(stim_vel>thrsvel);
        %         still = find(stim_vel<=thrsvel);
        %         %Remove moving stimulations
        %         ekeep_stim = e_stim;
        %         ekeep_stim(moving,:)=[];
        
        ekeep_stim = e_stim;
        
        %%%% CAN ALSO LOOK AT VELOCITY FOR LARGE TIME IN PAST
        %         for i=1:length(velindk)
        %             vtemp = currvel(velindk(i)-Fs*pastsecs+1:velindk(i));
        %             stim_meanvel(i)=mean(vtemp);
        %         end
        %
        %         ekeep_stim = etmp_stim;
        %         remove = find(stim_meanvel>thrsvel);
        %         ekeep_stim(remove,:)=[];
        
        %%      %% Figure for Current Day
        
        %clr={'r','k','g'};
        if figopt1==1
            figure(1); hold on;
            plot(taxis, mean(ekeep_stim),[clr{ep} '.-'],'Linewidth',2,'Markersize',6);
            uperr = mean(ekeep_stim) + std(ekeep_stim);
            lowerr = mean(ekeep_stim) - std(ekeep_stim);
            %line([0.2*1000:0.2*1000],[min(lowerr):max(uperr)],'Color','b','LineWidth',2,'LineStyle','--');
            jbfill(taxis,lowerr, uperr,[clr{ep}],[clr{ep}],1,0.2);
            
            if ep==1
                yplot = min(lowerr):1:max(uperr);
                xplot = 0*1000*ones(size(yplot));
                plot(xplot,yplot,'b--','Linewidth',2);
                %                 title(['Day ' num2str(day) ': Evoked Response (Tet' num2str(tet) '): Pre-Sleep(black) and Post-sleep (red&green)'],...
                %                     'FontSize',24,'Fontweight','bold');
                xlabel('Time (ms)');
                ylabel('LFP (uV)');
            end
            
            %axis([-100 150 -5000 1000])
            %set(gca,'XTickLabel',[]);
            %set(gca,'YTickLabel',[]);
            %axis off;
            
            text(200, -1000-500*(ep-1), ['Nstim:Sleep' num2str(ep) ' = ' num2str(size(ekeep_stim,1))],'FontSize',14,'Fontweight','bold');
            
            
            % Indiv EPSPs
            if indepsps==1
                figure(2); hold on;
                %redimscreen;
                redimscreen_widevert;
                orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
                set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','bold');
                set(0,'defaultaxeslinewidth',2);
                
                subplot(length(allepochs),1,ep); hold on;
                plot(taxis', ekeep_stim',[clr{ep} '-'],'Linewidth',2,'Markersize',6);
            end
            
            
        end  % end epoch
        
        %% Save epoch responses for comparison across days
        
        all_estim{day}{ep} = ekeep_stim;
        
        if ep==1
            e_pre = [e_pre; mean(ekeep_stim)];
            n_pre(day) = size(ekeep_stim,1);
            %e_pre = e_pre(:,1:900);
        end
        if ep==2
            e_post = [e_post; mean(ekeep_stim)];
            n_post(day) = size(ekeep_stim,1);
            %e_post = e_post(:,1:900);
        end
        
        if ep==3
            e_post2 = [e_post2; mean(ekeep_stim)];
            n_post2(day) = size(ekeep_stim,1);
            %e_post = e_post(:,1:900);
        end
        
        
        e_stim = []; etmp_stim = []; ekeep_stim=[]; stim_meanvel=[];
        ripamp_stim = []; ripampfilt_stim = [];
        ripenv_stim = []; ripenvfilt_stim = [];
        
    end  % end epoch
    
    if saveg1==1,
        orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        if isempty(ref)
            saveas(gcf,[prefix '_Day' num2str(day) '_Tet' num2str(tet) '_SleepEPSPResponse'],'fig');
            saveas(gcf,[prefix '_Day' num2str(day) '_Tet' num2str(tet) '_SleepEPSPResponse'],'jpg');
        else
            saveas(gcf,[prefix '_Day' num2str(day) '_Tet' num2str(tet) '_SleepEPSPResponse_WithRef'],'fig');
            saveas(gcf,[prefix '_Day' num2str(day) '_Tet' num2str(tet) '_SleepEPSPResponse_WithRef'],'jpg');
        end
        
    end
    
    % Bar Graph for day if asked for
    if indivbar==1
        
        for i=1:size(all_estim{day}{1},1)
            base = mean(all_estim{day}{1}(i,base_win_idx(1):base_win_idx(2)));
            extr = [max(all_estim{day}{1}(i,startidx:endidx)) min(all_estim{day}{1}(i,startidx:endidx)) ];
            extr = extr(find(abs(extr)==max(abs(extr)))); extr=extr(1);
            pre(i) = abs(extr-base);
        end
        for i=1:size(all_estim{day}{2},1)
            base = mean(all_estim{day}{2}(i,base_win_idx(1):base_win_idx(2)));
            extr = [max(all_estim{day}{2}(i,startidx:endidx)) min(all_estim{day}{2}(i,startidx:endidx)) ];
            extr = extr(find(abs(extr)==max(abs(extr)))); extr=extr(1);
            post1(i) = abs(extr-base);
        end
        if length(allepochs)>2
            for i=1:size(all_estim{day}{3},1)
                base = mean(all_estim{day}{3}(i,base_win_idx(1):base_win_idx(2)));
                extr = [max(all_estim{day}{3}(i,startidx:endidx)) min(all_estim{day}{3}(i,startidx:endidx)) ];
                extr = extr(find(abs(extr)==max(abs(extr)))); extr=extr(1);
                post2(i) = abs(extr-base);
            end
        end
        
        figure; hold on; redimscreen_figforppt1;
        %redimscreen_land;
        orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        
        bar(1,abs(mean(pre)),'k'); errorbar(1,abs(mean(pre)),abs(sem(pre)),'k');
        bar(2,abs(mean(post1)),'r'); errorbar(2,abs(mean(post1)),abs(sem(post1)),'r');
        if length(allepochs)>2
            bar(3,abs(mean(post2)),'g'); errorbar(3,abs(mean(post2)),abs(sem(post2)),'g');
        end
        if length(allepochs)>2
            set(gca,'XTick',[1 2 3],'XTickLabel',{'Pre','Post1','Post2'});
        else
            set(gca,'XTick',[1 2],'XTickLabel',{'Pre','Post1'});
        end
        
        [h1,p1,ci1] = ttest2(pre, post1, 0.05);
        [hk1,pk1,cik1] = kstest2(pre, post1, 0.05);
        
        if length(allepochs)>2
            [h2,p2,ci2] = ttest2(pre, post2, 0.05);
            [hk2,pk2,cik2] = kstest2(pre, post2, 0.05);
        end
        
        if h1==1,
            mul = sign(mean(post1));
            plot(2, mean(post1)+1.4*mul*sem(post1), 'r*','MarkerSize',8);
        end
        
        if h2==1,
            mul = sign(mean(post2));
            plot(3, mean(post2)+1.4*mul*sem(post2), 'r*','MarkerSize',8);
        end

        pre=[]; post1=[]; post2=[];
        
    end
    
    
end % end day


%%%%%%%%%%%%

if figopt2==1
    
    figure; hold on;
    redimscreen_land;
    orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    
    e_pre = e_pre(:,1:length(taxis));
    plot(taxis, mean(e_pre),['k.-'],'Linewidth',2,'Markersize',6);
    uperr = mean(e_pre) + std(e_pre);
    lowerr = mean(e_pre) - std(e_pre);
    jbfill(taxis,lowerr, uperr,'k','k',1,0.2);
    
    plot(taxis, mean(e_post),['r.-'],'Linewidth',2,'Markersize',6);
    uperr = mean(e_post) + std(e_post);
    lowerr = mean(e_post) - std(e_post);
    jbfill(taxis,lowerr, uperr,'r','r',1,0.2);
    
    text(200, -1000, ['Mean Nstim:Sleep1 = ' num2str(round(mean(n_pre)))],'FontSize',16,'Fontweight','bold');
    text(200, -1500, ['Mean Nstim:Sleep2 = ' num2str(round(mean(n_post)))],'FontSize',16,'Fontweight','bold');
    
    if length(allepochs)>2
        plot(taxis, mean(e_post2),['g.-'],'Linewidth',2,'Markersize',6);
        uperr = mean(e_post2) + std(e_post2);
        lowerr = mean(e_post2) - std(e_post2);
        jbfill(taxis,lowerr, uperr,'g','g',1,0.2);
        text(200, -2000, ['Mean Nstim:Sleep3 = ' num2str(round(mean(n_post2)))],'FontSize',16,'Fontweight','bold');
    end
    
    yplot = min(lowerr):1:max(uperr);
    xplot = 0*1000*ones(size(yplot));
    plot(xplot,yplot,'b--','Linewidth',2);
    
    title(['Avg over Days ' num2str(min(days)) ':' num2str(max(days)) ': Evoked Response (Tet' num2str(tet) '): Pre-Sleep(black) and Post-sleep (red&green)'],...
        'FontSize',24,'Fontweight','bold');
    xlabel('Time (ms)');
    ylabel('LFP (uV)');
    
    
    
    if saveg2==1,
        orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        if isempty(ref)
            saveas(gcf,[prefix 'Days' num2str(min(days)) 'to' num2str(max(days)) '_Tet' num2str(tet) '_SleepEPSPResponse'],'jpg');
            saveas(gcf,[prefix 'Days' num2str(min(days)) 'to' num2str(max(days)) '_Tet' num2str(tet) '_SleepEPSPResponse'],'fig');
        else
            saveas(gcf,[prefix 'Days' num2str(min(days)) 'to' num2str(max(days)) '_Tet' num2str(tet) '_SleepEPSPResponse_WithRef'],'jpg');
            saveas(gcf,[prefix 'Days' num2str(min(days)) 'to' num2str(max(days)) '_Tet' num2str(tet) '_SleepEPSPResponse_WithRef'],'fig');
        end
    end
    
end

%%% Statistics on Peak Response

preall=[]; post1all=[]; post2all=[];

for d=1:length(days),
    
    % If range is -200 to 400, old method
    %     pre(i) = min(e_pre(i,301:345));
    %     post1(i) = min(e_post(i,301:345));
    %     if length(allepochs)>2
    %         post2(i) = min(e_post2(i,301:345));
    %     end
    day=days(d);
   
    
    for i=1:size(all_estim{day}{1},1)
        base = mean(all_estim{day}{1}(i,base_win_idx(1):base_win_idx(2)));
        extr = [max(all_estim{day}{1}(i,startidx:endidx)) min(all_estim{day}{1}(i,startidx:endidx)) ];       
        extr = extr(find(abs(extr)==max(abs(extr))));  extr=extr(1);    
        preall(d,i) = abs(extr-base);
        
%         if sign(extr)==sign(base)
%             preall(d,i)=abs(extr)-abs(base);
%         else
%             preall(d,i)=extr+base;
%         end
%preall(d,i) = max(abs(all_estim{day}{1}(i,startidx:endidx)));
    end
    for i=1:size(all_estim{day}{2},1)
        base = mean(all_estim{day}{2}(i,base_win_idx(1):base_win_idx(2)));
        extr = [max(all_estim{day}{2}(i,startidx:endidx)) min(all_estim{day}{2}(i,startidx:endidx)) ];
        extr=extr(find(abs(extr)==max(abs(extr))));  extr=extr(1);
        post1all(d,i) = abs(extr-base);
    end
    if length(allepochs)>2
        for i=1:size(all_estim{day}{3},1)
            base = mean(all_estim{day}{3}(i,base_win_idx(1):base_win_idx(2)));
            extr = [max(all_estim{day}{3}(i,startidx:endidx)) min(all_estim{day}{3}(i,startidx:endidx)) ];
            extr = extr(find(abs(extr)==max(abs(extr))));  extr=extr(1);
            post2all(d,i) = abs(extr-base);
        end
    end
    
end


if figopt2==1
    figure(22); hold on;
    %redimscreen_land;
    orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    
    bar(1,abs(mean(preall(:))),'k'); errorbar(1,abs(mean(preall(:))),abs(sem(preall(:))),'k');
    bar(2,abs(mean(post1all(:))),'r'); errorbar(2,abs(mean(post1all(:))),abs(sem(post1all(:))),'r');
    if length(allepochs)>2
        bar(3,abs(mean(post2all(:))),'g'); errorbar(3,abs(mean(post2all(:))),abs(sem(post2all(:))),'g');
    end
    if length(allepochs)>2
        set(gca,'XTick',[1 2 3],'XTickLabel',{'Pre','Post1','Post2'});
    else
        set(gca,'XTick',[1 2],'XTickLabel',{'Pre','Post1'});
    end
    
    [h1,p1,ci1] = ttest2(preall(:), post1all(:), 0.05);
    [hk1,pk1,cik1] = kstest2(preall(:), post1all(:), 0.05);
    
    if length(allepochs)>2
        [h2,p2,ci2] = ttest2(preall(:), post2all(:), 0.05, 'right');
        [hk2,pk2,cik2] = kstest2(preall(:), post2all(:), 0.05, 'smaller');
    end
    
    
    if h1==1,
        mul = sign(mean(post1all(:)));
        plot(2, mean(post1all(:))+1.4*mul*sem(post1all(:)), 'r*','MarkerSize',8);
    end
    
    if h2==1,
        mul = sign(mean(post2all(:)));
        plot(3, mean(post2all(:))+1.4*mul*sem(post2all(:)), 'r*','MarkerSize',8);
    end
    
end
%%%%%%%%%%%%%%% %%%%%%%%

if savedata==1,
    
    if length(days)>1
%         if isempty(ref)
%             savefile = sprintf('%s/%s_EPSP_Days%01dto%01d_Tet%02d.mat', savedir, prefix, min(days), max(days), tet);
%             %savefile = [prefix '_LTP_Days' num2str(min(days)) ':' num2str(max(days)) 'Tet' num2str(tet)];
%             save(savefile, 'all_estim');
%         else
%             savefile = sprintf('%s/%s_EPSP_Days%01dto%01d_Tet%02d_WithRef.mat', savedir, prefix, min(days), max(days), tet);
%             %savefile = [prefix '_LTP_Days' num2str(min(days)) ':' num2str(max(days)) 'Tet' num2str(tet)];
%             save(savefile, 'all_estim');
%         end
        
        % Change Names
        if isempty(ref)
            savefile = sprintf('%s/%s_EPSP_Tet%02d.mat', savedir, prefix, tet);
            %savefile = [prefix '_LTP_Days' num2str(min(days)) ':' num2str(max(days)) 'Tet' num2str(tet)];
            save(savefile, 'all_estim','days');
        else
            savefile = sprintf('%s/%s_EPSP_Tet%02d_WithRef.mat', savedir, prefix, tet);
            %savefile = [prefix '_LTP_Days' num2str(min(days)) ':' num2str(max(days)) 'Tet' num2str(tet)];
            save(savefile, 'all_estim','days');
        end
        
    end
    
    
end



