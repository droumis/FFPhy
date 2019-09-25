
function [ep] = sj_epsp_slope_day (animdirect,prefix,day,allepochs,tets,ref, figopt1,saveg1,figopt2,figopt3,savedata)

% From sj_epsp_slope_day
% From sj_epsp_amp_tet /sj_stimresp2_withvel2 on 18 May 2011
% EPSPs for given daye across tetrodes. Plots Summary Bar Graph of EPSPSlope.

% Ampl: Includes omission of non-existent stimulation, and accurate calc on EPSP ampl based
% on narrow window and which peak comes first and within window

% PLOT and compare Response to Probe Stimulation in different Sleep Epochs across tets
% Egs
% sj_epsp_slope_day('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',3,[1 3 5],5,[],1,0,1,1,0);
% sj_epsp_amp_day('/data25/sjadhav/RippleInterruption/REd_direct','REd',1,[1 3 5],[2 3 4 5],[],1,0,2,0,0);
% sj_stimresp2_withvel('/data25/sjadhav/RippleInterruption/RE1_direct','RE1',1:9,[1 3 5],5,3,0,0,1,0,0);
% sj_stimresp2_withvel_longer('/data25/sjadhav/RippleInterruption/REd_direct','REd',5,[1 3 5],3,1,1,0,0,0,0);
% sj_stimresp2_withvel_longer('/data25/sjadhav/RippleInterruption/REd_direct','REd',5,[1 3 5],3,[],1,0,0,0,0);
% sj_stimresp2_withvel('/data25/sjadhav/RippleInterruption/REf_direct','REe',4,[1 3 5],11,9,1,0,0,0,0);
% sj_stimresp2_withvel('/data25/sjadhav/RippleInterruption/SJStimC_direct','sjc',3:7,[1 3],6,3,0,0,1,0,0);

% figopt1: Make EPSP plot for individual tet
% saveg1: Save individual tet graphs
% figopt2: Make bar plot for individual tetr
% figopt3: Average for day across tetrodes
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
    figopt3 = 0;
end
if nargin<11,
    savedata = 0;
end


% Indiv plot options
indepsps=0;

% For slope, always need figure
if figopt1==0
    figopt1=1;
end


savedir = '/data25/sjadhav/RippleInterruption/ProcessedData/EPSP/Slope';

% Fixed parameters
Fs=30; %video sampling rate
tsamp=1/30; % in sec
velfiltlth = Fs/4; % Filter for smoothing velocity: 1 sec=Fs points, Std Dev = Lth/4

% Variable parameters
thrsvel=50;  %<x cm per sec is still
pastsecs = 1;

pre_frac=0.02; %eg. 0.2*e.samprate: corresponds to secs eg. 0.2 = 200ms
post_frac=0.03; %eg. 0.4*e.samprate



e.samprate=1500;
startidx = round(pre_frac*e.samprate)+7; %For 100 Hz sampling rate, 5ms = 7.5 points. Stim artifact over by 3 ms
endidx = round(pre_frac*e.samprate)+75;  % Peak generally at 10-15ms, EPSP over by 40-50 ms, 50 ms = 75 points
% WIndow for EPSP Ampl
% OR, better way - go by time
stimtime=pre_frac*e.samprate*(1000/e.samprate);
stimtime_idx=pre_frac*e.samprate;
resp_win = stimtime + [7, 14]; % ms
resp_win_idx = floor(resp_win*e.samprate/1000);

base_win = stimtime - [19, 1];
base_win_idx = floor(base_win*e.samprate/1000);

startidx = resp_win_idx(1);
endidx = resp_win_idx(2);


% Loop over tets and load files
directoryname = animdirect;
if (animdirect(end) == '/')
    animdirect = animdirect(1:end-1);
end
cd(animdirect);

% ------------------------------
% Figure and Font Sizes

forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/EPSP/Slope/';
summdir = figdir;
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

if forppr==1
    set(0,'defaultaxesfontsize',16);
    tfont = 18; % title font
    xfont = 16;
    yfont = 16;
else
    set(0,'defaultaxesfontsize',24);
    tfont = 28;
    xfont = 20;
    yfont = 20;
end

clr = {'b','r','g','c','m','y','k','b','r','g'};

% Saving figures
savefig1 = 0; % for single tet
savefig2 = 0; % for multiple tets


% ---------------------------------------


% Loop Across Tetrodes

for n = 1:length(tets)
    
    tet = tets(n);
    % Day n
    DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
    load(DIOfile);
    
    % Pos file
    posfile = sprintf('%s/%spos%02d.mat', directoryname, prefix, day);
    load(posfile);
    
    clr = {'k','r','g','c','b','m'};
    
    if figopt1==1
        figure(n); hold on;
        %redimscreen;
        %redimscreen_land;
        redimscreen_2versubplots;
    end
    
    % Epoch 1,3,5... - Sleep 1,2,3...
    % Loop over epochs
    for ep = 1:length(allepochs)
        
        e_stim = [];
        epoch = allepochs(ep);
        
        stim = DIO{day}{epoch}{15};
        if isempty(stim)
            stim = DIO{day}{epoch}{16};
        end
        
        stim_starttime = stim.pulsetimes(:,1);
        stim_endtime = stim.pulsetimes(:,2);
        stim_length = stim.pulselength;
        stim_isi = stim.timesincelast(2:end);
        
        
        
        %% EEg file for current tet
        
        EEGfile = sprintf('%s/EEG/%seeg%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tet);
        load(EEGfile);
        e = eeg{day}{epoch}{tet};
        ti = geteegtimes(e);
        pt = stim.pulsetimes ./ 10000; % in sec
        eind = lookup(pt(:,1), ti);
        
        % EEg File for Reference
        if ~isempty(ref)
            EEGfile_ref = sprintf('%s/EEG/%seeg%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,ref);
            load(EEGfile_ref);
            e_ref = eeg{day}{epoch}{ref};
        end
        
        for i =1:length(eind)-2;
            e_tempstim=e.data(eind(i+1)-round(pre_frac*e.samprate):eind(i+1)+round(post_frac*e.samprate));
            if ~isempty(ref)
                e_stim_ref=e_ref.data(eind(i+1)-round(pre_frac*e.samprate):eind(i+1)+round(post_frac*e.samprate));
                e_stim(i,:)=e_tempstim+e_stim_ref;
            else
                e_stim(i,:)=e_tempstim;
            end
        end
        taxis = [1:size(e_stim,2)]*1000/e.samprate;
        taxis = taxis - pre_frac*1000;
        ekeep_stim = e_stim;
        
        
        % Figure for Current Tet
        
        if figopt1==1
            figure(n); hold on;
            subplot(2,1,1); hold on;
            % Plot with time axis already on
            plot(taxis, mean(ekeep_stim),[clr{ep} '.-'],'Linewidth',2,'Markersize',6);
            uperr = mean(ekeep_stim) + std(ekeep_stim);
            lowerr = mean(ekeep_stim) - std(ekeep_stim);
            jbfill(taxis,lowerr, uperr,[clr{ep}],[clr{ep}],1,0.2);
            
            if ep==1
                % Plots with time axis
                yplot = min(lowerr):1:max(uperr);
                
                xplot = 0*ones(size(yplot));
                plot(xplot,yplot,'b--','Linewidth',2);
                xplot = (resp_win(1)-pre_frac*1000)*ones(size(yplot));
                
                plot(xplot,yplot,'c--','Linewidth',2);
                xplot = (resp_win(2)-pre_frac*1000)*ones(size(yplot));
                plot(xplot,yplot,'c--','Linewidth',2);
                
                xplot = (base_win(1)-pre_frac*1000)*ones(size(yplot));
                plot(xplot,yplot,'g--','Linewidth',2);
                xplot = (base_win(2)-pre_frac*1000)*ones(size(yplot));
                plot(xplot,yplot,'g--','Linewidth',2);
                
                xlabel('Time (ms)');
                ylabel('LFP (uV)');
                
                % Plot Peak
                % Get Direction of Mean EPSP and Plot Peak Locn/Ampl for Current Tetrode
                % Find out whether it is a maxima or minima for given day from mean response in EACH epoch
                meanresp = mean(ekeep_stim(:,startidx:endidx));
                
                % Method of Peaks - better
                peakmaxidx=[peakpick(meanresp, 'max')];
                peakminidx=[peakpick(meanresp, 'min')];
                if isempty(peakmaxidx), dirn='neg'; peakidx=peakminidx; meanmax = abs(min(meanresp)); mul=-1; end;
                if isempty(peakminidx), dirn='pos'; peakidx=peakmaxidx; meanmax = abs(max(meanresp)); mul=1; end;
                
                % Find bigger peak if both exist - not earliest - that is problematic. In the small response window, there should be only 1 big peak
                if (~isempty(peakmaxidx)) && (~isempty(peakminidx)),
                    [~,idx]= max([abs(meanresp(peakmaxidx)) abs(meanresp(peakminidx))]);
                    if idx==1,
                        dirn='pos'; meanmax = abs(max(meanresp)); peakidx=peakmaxidx; mul=1;
                    else
                        dirn='neg'; meanmax = abs(min(meanresp)); peakidx=peakminidx; mul=-1;
                    end
                end
                
                if length(peakidx)>1, peakidx=peakidx(find(abs(meanresp(peakidx))==meanmax)); end
                plot((peakidx*1000/e.samprate)+(resp_win(1)-stimtime-1),meanmax*mul,'bo','MarkerSize',12,'Linewidth',4);
            end
            
            text(200, -1000-500*(ep-1), ['Nstim:Sleep' num2str(ep) ' = ' num2str(size(ekeep_stim,1))],'FontSize',14,'Fontweight','bold');
            title (['Tet ' num2str(tet)],'FontSize',14,'Fontweight','bold')
            
            % Get Window for slope with 2nd figure on subplot2 with idx axis
            subplot(2,1,2); hold on;
            plot(mean(ekeep_stim),[clr{ep} '.-'],'Linewidth',1,'Markersize',18);
            
            if ep == length(allepochs)
                
                xlabel('Idx');
                ylabel('LFP (uV)');
                grid on;
                set(gca,'XLim',[stimtime_idx-2 stimtime_idx+30]);
                yplot = round(min(mean(ekeep_stim))):1:round(max(mean(ekeep_stim)));
                
                xplot = stimtime_idx*ones(size(yplot));
                plot(xplot,yplot,'b--','Linewidth',2);
                xplot = resp_win_idx(1)*ones(size(yplot));
                plot(xplot,yplot,'c--','Linewidth',1);
                xplot = resp_win_idx(2)*ones(size(yplot));
                plot(xplot,yplot,'c--','Linewidth',2);
                
                disp ('Choose Window for Slope');
                [rng,lims]=ginput;
                if exist('rng')==0
                    disp ('You did not pick window correctly');
                    return
                else
                    rng = sort(round(rng));
                    lims=sort(lims); lims=lims';
                end
            end
 
        end  % end Figure(t)
        
        % Plot Indiv EPSPs to check non-existent/delayed stimulation
        if indepsps==1
            figure(200+n); hold on;
            %redimscreen;
            redimscreen_widevert;
            orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
            set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','bold');
            set(0,'defaultaxeslinewidth',2);
            
            subplot(length(allepochs),1,ep); hold on;
            plot(taxis', ekeep_stim',[clr{ep} '-'],'Linewidth',2,'Markersize',6);
        end
        
        
        %Save epoch responses for comparison across tets
        all_estim{day}{ep}{tet} = ekeep_stim;
        
        %Reset
        e_stim = []; ekeep_stim=[];
        
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
    
    
    % SLOPE CALCULATION FOR CURRENT TETR and Statistics on Slope Response
    
    % KEEP AMP CALCULN TO SKIP FAILED STIMULATIONS
    % Find out whether it is a maxima or minima for current tet on given day from mean response in 1ST epoch
    meanresp = mean(all_estim{day}{1}{tet}(:,startidx:endidx));
    
    % Old method of extrema - some issues
    extr = [max(meanresp) min(meanresp)];
    extridx = abs(extr)==max(abs(extr));
    meanmax=abs(extr(extridx));
    
    % Method of Peaks - better
    peakmaxidx=[peakpick(meanresp, 'max')];
    peakminidx=[peakpick(meanresp, 'min')];
    
    if isempty(peakmaxidx), dirn='neg'; peakidx=peakminidx; meanmax = abs(min(meanresp)); end;
    if isempty(peakminidx), dirn='pos'; peakidx=peakmaxidx; meanmax = abs(max(meanresp)); end;
    % Find bigger peak if both exist - not earliest - that is problematic. In the small response window, there should be only 1 big peak
    if (~isempty(peakmaxidx)) && (~isempty(peakminidx)),
        [~,idx]= max([abs(meanresp(peakmaxidx)) abs(meanresp(peakminidx))]);
        if idx==1,
            dirn='pos'; meanmax = abs(max(meanresp)); peakidx=peakmaxidx;
        else
            dirn='neg'; meanmax = abs(min(meanresp)); peakidx=peakminidx;
        end
    end

    % Get x-range for slope for current tetrode
    x1 = rng(1):rng(2);
    X = [ones(size(x1))' x1'];
    
    %% Actual Calculations
    pre=[]; post1=[]; post2=[];
    cnt=0;
    for i=1:size(all_estim{day}{1}{tet},1)
        % Amp
        base = mean(all_estim{day}{1}{tet}(i,base_win_idx(1):base_win_idx(2)));
        if dirn=='pos'
            extr = [max(all_estim{day}{1}{tet}(i,startidx:endidx))];
        else
            extr = [min(all_estim{day}{1}{tet}(i,startidx:endidx)) ];
        end
        val=abs(extr-base);
        
        % Slope
        % Direct
        del_y = all_estim{day}{1}{tet}(i,rng(2)) - all_estim{day}{1}{tet}(i,rng(1));
        del_x = (rng(2) - rng(1))*1000./e.samprate; % idx to ms
        slope1 = abs(del_y./del_x);
        
        % By Regression     
        b = regress(all_estim{day}{1}{tet}(i,rng(1):rng(2))',X);
        slope = abs(b(2))*e.samprate/1000; % mV/ms
        
        % Check for failed stimulation
        if (val>=0.1*meanmax) && (val>100)
            cnt=cnt+1;
            pre(cnt) = slope;
            preall{tet}(cnt) = slope;
            all_eslope{day}{1}{tet}(cnt) = slope;
        end
        
    end
    cnt=0;
    for i=1:size(all_estim{day}{2}{tet},1)
        base = mean(all_estim{day}{2}{tet}(i,base_win_idx(1):base_win_idx(2)));
        extr = [max(all_estim{day}{2}{tet}(i,startidx:endidx)) min(all_estim{day}{2}{tet}(i,startidx:endidx)) ];
        %             extr = extr(find(abs(extr)==max(abs(extr)))); extr=extr(1);
        if dirn=='pos', extr=extr(1); else extr=extr(2); end
        val=abs(extr-base);
        
         % Slope
         % Direct
        del_y = all_estim{day}{2}{tet}(i,rng(2)) - all_estim{day}{2}{tet}(i,rng(1));
        del_x = (rng(2) - rng(1))*1000./e.samprate; % idx to ms
        slope1 = abs(del_y./del_x);
        
        % By Regression
        b = regress(all_estim{day}{2}{tet}(i,rng(1):rng(2))',X);
        slope = abs(b(2))*e.samprate/1000;
        
        if (val>=0.1*meanmax) && (val>100)
            cnt=cnt+1;
            post1(cnt) = slope;
            post1all{tet}(cnt) = slope;
            all_eslope{day}{2}{tet}(cnt) = slope;
        end
        
    end
    cnt=0;
    if length(allepochs)>2
        for i=1:size(all_estim{day}{3}{tet},1)
            base = mean(all_estim{day}{3}{tet}(i,base_win_idx(1):base_win_idx(2)));
            extr = [max(all_estim{day}{3}{tet}(i,startidx:endidx)) min(all_estim{day}{3}{tet}(i,startidx:endidx)) ];
            %                 extr = extr(find(abs(extr)==max(abs(extr)))); extr=extr(1);
            if dirn=='pos', extr=extr(1); else extr=extr(2); end
            val=abs(extr-base);
            
             % Slope
             % Direct
            del_y = all_estim{day}{3}{tet}(i,rng(2)) - all_estim{day}{3}{tet}(i,rng(1));
            del_x = (rng(2) - rng(1))*1000./e.samprate; % idx to ms
            slope1 = abs(del_y./del_x);
        
            % By Regression
            b = regress(all_estim{day}{3}{tet}(i,rng(1):rng(2))',X);
            slope = abs(b(2))*e.samprate/1000;
        
            if (val>=0.1*meanmax) && (val>100)
                cnt=cnt+1;
                post2(cnt) = slope;
                post2all{tet}(cnt) = slope;
                all_eslope{day}{3}{tet}(cnt) = slope;
            end
        end
    end
    
    % Bar Graph for tet if asked for - Can Plot Raw EPSP Slope (1) or Z-score (2)
    if figopt2==1 || figopt2==2
        figure(100+n); hold on; 
        redimscreen_2horsubplots
        %redimscreen_figforppt1;
        %redimscreen_land;
        
        % Raw EPSP Amp
        if figopt2==1
            
            subplot(1,2,1); hold on;
            bar(1,abs(mean(pre)),'r'); errorbar(1,abs(mean(pre)),abs(sem(pre)),'k');
            bar(2,abs(mean(post1)),'r'); errorbar(2,abs(mean(post1)),abs(sem(post1)),'k');
            if length(allepochs)>2
                bar(3,abs(mean(post2)),'r'); errorbar(3,abs(mean(post2)),abs(sem(post2)),'k');
            end
            if length(allepochs)>2
                set(gca,'XTick',[1 2 3],'XTickLabel',{'Pre','Post1','Post2'});
            else
                set(gca,'XTick',[1 2],'XTickLabel',{'Pre','Post1'});
            end
            ylabel('Slope (mV/ms)');
            [h1,p1,ci1] = ttest2(pre, post1, 0.05);
            [hk1,pk1,cik1] = kstest2(pre, post1, 0.05);
            if length(allepochs)>2
                [h2,p2,ci2] = ttest2(pre, post2, 0.05);
                [hk2,pk2,cik2] = kstest2(pre, post2, 0.05);
            end
            title (['Tet ',num2str(tet),': p1=',num2str(roundn(p1,-3)),'; p2=',num2str(roundn(p2,-3))],'FontSize',tfont,'Fontweight','normal')
            if h1==1,
                mul = sign(mean(post1));
                plot(2, abs(mean(post1))+1.4*sem(post1), 'r*','MarkerSize',8);
            end
            if h2==1,
                mul = sign(mean(post2));
                plot(3, abs(mean(post2))+1.4*sem(post2), 'r*','MarkerSize',8);
            end
            
            subplot(1,2,2); hold on;
            post=[post1,post2];
            bar(1,abs(mean(pre)),'r'); errorbar(1,abs(mean(pre)),abs(sem(pre)),'k');
            bar(2,abs(mean(post)),'r');
            errorbar(2,abs(mean(post)),abs(sem(post)),'k');
            set(gca,'XTick',[1 2],'XTickLabel',{'Pre','Post'});
            ylabel('Slope (uV/ms)');
            [h1,p1,ci1] = ttest2(pre, post, 0.05);
            [hk1,pk1,cik1] = kstest2(pre, post, 0.05);
            title (['Tet ',num2str(tet),': p = ',num2str(roundn(p1,-3))],'FontSize',tfont,'Fontweight','normal')
            if h1==1,
                mul = sign(mean(post));
                plot(2, abs(mean(post))+1.4*abs(sem(post)), 'r*','MarkerSize',8);
            end
            if savefig1==1,
                figfile = [figdir,'EpspSlope_',prefix,'d',num2str(day),'tet',num2str(tets)];
                print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
            end
        else % Z-score
            base_mean=mean(pre); base_std = std(pre);
            preZ = (pre - repmat(base_mean,size(pre)))./repmat(base_std,size(pre));
            post1Z = (post1 - repmat(base_mean,size(post1)))./repmat(base_std,size(post1));
            post2Z = (post2 - repmat(base_mean,size(post2)))./repmat(base_std,size(post2));
            postZ = [post1Z,post2Z];
            
            subplot(1,2,1); hold on;
            bar(1,mean(preZ),'r'); errorbar(1,mean(preZ),sem(preZ),'k');
            bar(2,mean(post1Z),'r'); errorbar(2,mean(post1Z),sem(post1Z),'k');
            if length(allepochs)>2
                bar(3,mean(post2Z),'r'); errorbar(3,mean(post2Z),sem(post2Z),'k');
            end
            if length(allepochs)>2
                set(gca,'XTick',[1 2 3],'XTickLabel',{'Pre','Post1','Post2'});
            else
                set(gca,'XTick',[1 2],'XTickLabel',{'Pre','Post1'});
            end
            [h1,p1,ci1] = ttest2(preZ, post1Z, 0.05);
            [hk1,pk1,cik1] = kstest2(preZ, post1Z, 0.05);
            if length(allepochs)>2
                [h2,p2,ci2] = ttest2(preZ, post2Z, 0.05);
                [hk2,pk2,cik2] = kstest2(preZ, post2Z, 0.05);
            end
            title (['Tet ',num2str(tet),': p1=',num2str(roundn(p1,-3)),'; p2=',num2str(roundn(p2,-3))],'FontSize',tfont,'Fontweight','normal');
            ylabel('Z-score Slope');
            if h1==1,
                mul = sign(mean(post1Z));
                plot(2, mean(post1Z)+1.4*mul*sem(post1Z), 'r*','MarkerSize',8);
            end
            if h2==1,
                mul = sign(mean(post2Z));
                plot(3, mean(post2Z)+1.4*mul*sem(post2Z), 'r*','MarkerSize',8);
            end
            
            subplot(1,2,2); hold on;
            bar(1,abs(mean(preZ)),'r'); errorbar(1,abs(mean(preZ)),abs(sem(preZ)),'k');
            bar(2,abs(mean(postZ)),'r');
            errorbar(2,abs(mean(postZ)),abs(sem(postZ)),'k');
            set(gca,'XTick',[1 2],'XTickLabel',{'Pre','Post'});
            ylabel('Slope (uV/ms)');
            [h1,p1,ci1] = ttest2(preZ, postZ, 0.05);
            [hk1,pk1,cik1] = kstest2(preZ, postZ, 0.05);
            title (['Tet ',num2str(tet),': p = ',num2str(roundn(p1,-3))],'FontSize',tfont,'Fontweight','normal')
            if h1==1,
                mul = sign(mean(postZ));
                plot(2, abs(mean(postZ))+1.4*abs(sem(postZ)), 'r*','MarkerSize',8);
            end
            if savefig1==1,
                figfile = [figdir,'EpspSlopeZscore_',prefix,'d',num2str(day),'tet',num2str(tets)];
                print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
            end
        end % end Z-score
    end % figopt2
end % end tet


%% Plot Summary Bar Graph of EPSP Slope Across Tetrodes for Given Day

% Vectorize raw slope across-tet data and get Z-scores also
preallvec=[]; post1allvec=[]; post2allvec=[]; postallvec=[];
respZep1=[]; respZep2=[]; respZep3=[]; respZep23=[];

for n=1:length(tets),
    
    tet=tets(n); respZ=[];
    
    %Raw Ampl
    preallvec=[preallvec;preall{tet}(:)];
    post1allvec=[post1allvec;post1all{tet}(:)];
    post2allvec=[post2allvec;post2all{tet}(:)];
    postallvec=[postallvec;post1all{tet}(:);post2all{tet}(:)];
    
    % Z-scores
    base_mean = mean(all_eslope{day}{1}{tet});
    base_std = std(all_eslope{day}{1}{tet});
    
    % For current tet, get z-scores
    for ep = 1:length(allepochs)
        currslope = all_eslope{day}{ep}{tet};
        respZ{ep} = (currslope - repmat(base_mean,size(currslope)))./repmat(base_std,size(currslope));
        respZsave{day}{ep}{tet} = (currslope - repmat(base_mean,size(currslope)))./repmat(base_std,size(currslope));
    end % end epoch
    
    respZep1 = [respZep1, respZ{1}];
    respZep2 = [respZep2, respZ{2}];
    respZep3 = [respZep3, respZ{3}];
    respZep23 = [respZep23, respZ{2}, respZ{3}];
end


if figopt3==1
    
    % Plot Summary Bar Graph of Raw Slope For Day Across Tetrodes
    figure; hold on; %redimscreen_figforppt1;
    redimscreen_2horsubplots
    subplot(1,2,1); hold on;
    bar(1,abs(mean(preallvec)),'r'); errorbar(1,abs(mean(preallvec)),abs(sem(preallvec)),'k');
    bar(2,abs(mean(post1allvec)),'r'); errorbar(2,abs(mean(post1allvec)),abs(sem(post1allvec)),'k');
    if length(allepochs)>2
        bar(3,abs(mean(post2allvec)),'r'); errorbar(3,abs(mean(post2allvec)),abs(sem(post2allvec)),'k');
    end
    if length(allepochs)>2
        set(gca,'XTick',[1 2 3],'XTickLabel',{'Pre','Post1','Post2'});
    else
        set(gca,'XTick',[1 2],'XTickLabel',{'Pre','Post'});
    end
    ylabel('EPSP Slope (uV/ms)','FontSize',yfont,'Fontweight','normal');
    [h1,p1,ci1] = ttest2(preallvec, post1allvec, 0.05);
    [hk1,pk1,cik1] = kstest2(preallvec, post1allvec, 0.05);
    if length(allepochs)>2
        [h2,p2,ci2] = ttest2(preallvec, post2allvec, 0.05, 'right');
        [hk2,pk2,cik2] = kstest2(preallvec, post2allvec, 0.05, 'smaller');
    end
    title (['All Tets ',num2str(tet),': p1=',num2str(roundn(p1,-3)),'; p2=',num2str(roundn(p2,-3))],'FontSize',tfont,'Fontweight','normal');
    if h1==1,
        mul = sign(mean(post1allvec));
        plot(2, mean(post1allvec)+1.4*mul*sem(post1allvec), 'r*','MarkerSize',8);
    end
    if h2==1,
        mul = sign(mean(post2allvec));
        plot(3, mean(post2allvec)+1.4*mul*sem(post2allvec), 'r*','MarkerSize',8);
    end
    
    subplot(1,2,2); hold on;
    bar(1,abs(mean(preallvec)),'r'); errorbar(1,abs(mean(preallvec)),abs(sem(preallvec)),'k');
    bar(2,abs(mean(postallvec)),'r'); errorbar(2,abs(mean(postallvec)),abs(sem(postallvec)),'k');
    set(gca,'XTick',[1 2],'XTickLabel',{'Pre','Post'});
    ylabel('EPSP Slope (uV/ms)','FontSize',yfont,'Fontweight','normal');
    [h1,p1,ci1] = ttest2(preallvec, post1allvec, 0.05);
    [hk1,pk1,cik1] = kstest2(preallvec, post1allvec, 0.05);
    title (['All Tets ',num2str(tet),': p1=',num2str(roundn(p1,-3))],'FontSize',tfont,'Fontweight','normal');
    if h1==1,
        mul = sign(mean(post1allvec));
        plot(2, mean(post1allvec)+1.4*mul*sem(post1allvec), 'r*','MarkerSize',8);
    end
    if h2==1,
        mul = sign(mean(post2allvec));
        plot(3, mean(post2allvec)+1.4*mul*sem(post2allvec), 'r*','MarkerSize',8);
    end
    
    if savefig2==1,
        figfile = [figdir,'EpspSlope_',prefix,'d',num2str(day),'alltets'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
    
    %% Plot Summary Bar Graph of Z-scores For Day Across Tetrodes
    
    figure; hold on;  %redimscreen_figforppt1;
    redimscreen_2horsubplots
    subplot(1,2,1); hold on;
    bar(1,mean(respZep1),'r'); errorbar(1,mean(respZep1),sem(respZep1),'r');
    bar(2,mean(respZep2),'r'); errorbar(2,mean(respZep2),sem(respZep2),'r');
    bar(3,mean(respZep3),'r'); errorbar(3,mean(respZep3),sem(respZep3),'r');
    set(gca,'XTick',[1 2 3],'XTickLabel',{'Pre','Post1','Post2'},'FontSize',xfont,'Fontweight','normal');
    ylabel('EPSP Slope Z-score','FontSize',yfont,'Fontweight','normal');
    axis([0 4 -0.5 0.5])
    [ht12m,pt12m,cit12m] = ttest2(respZep1, respZep2, 0.05);
    [ht13m,pt13m,cit13m] = ttest2(respZep1, respZep3, 0.05);
    [hk12m, pk12m, cik12m] = kstest2(respZep1, respZep2, 0.05, 'unequal');
    [hk13m, pk13m, cik13m] = kstest2(respZep1, respZep3, 0.05, 'unequal');
    title (['All Tets ',num2str(tet),': pt12m=',num2str(roundn(p1,-3)),'; p2=',num2str(roundn(pt13m,-3))],'FontSize',tfont,'Fontweight','normal');
    if ht12m==1,
        mul = sign(mean(respZep2));
        plot(2, mean(respZep2)+1.4*mul*sem(respZep2), 'r*','MarkerSize',8);
    end
    if ht13m==1,
        mul = sign(mean(respZep3));
        plot(3, mean(respZep3)+1.4*mul*sem(respZep3), 'r*','MarkerSize',8);
    end
    
    subplot(1,2,2); hold on;
    bar(1,mean(respZep1),'r'); errorbar(1,mean(respZep1),sem(respZep1),'r');
    bar(2,mean(respZep23),'r'); errorbar(2,mean(respZep23),sem(respZep23),'r');
    set(gca,'XTick',[1 2],'XTickLabel',{'Pre','Post'},'FontSize',xfont,'Fontweight','normal');
    ylabel('EPSP Slope Z-score','FontSize',yfont,'Fontweight','normal');
    axis([0 4 -0.5 0.5]);
    [ht12m,pt12m,cit12m] = ttest2(respZep1, respZep23, 0.05);
    [hk12m, pk12m, cik12m] = kstest2(respZep1, respZep23, 0.05, 'unequal');
    title (['All Tets ',num2str(tet),': p1=',num2str(roundn(pt12m,-3))],'FontSize',tfont,'Fontweight','normal');
    if ht12m==1,
        mul = sign(mean(respZep2));
        plot(2, mean(respZep2)+1.4*mul*sem(respZep2), 'r*','MarkerSize',8);
    end
    if ht13m==1,
        mul = sign(mean(respZep3));
        plot(3, mean(respZep3)+1.4*mul*sem(respZep3), 'r*','MarkerSize',8);
    end
    
    if savefig2==1,
        figfile = [figdir,'EpspSlopeZscore_',prefix,'d',num2str(day),'alltets'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
end
%% Save Data

if savedata==1,
    
    if isempty(ref)
        savefile = sprintf('%s/%s_EPSP_Day%02d.mat', savedir, prefix, day);
        save(savefile, 'all_estim','day','tets','all_eslope','respZsave');
    else
        savefile = sprintf('%s/%s_EPSP_Day%02d_WithRef.mat', savedir, prefix, day);
        save(savefile, 'all_estim','day','tets','all_eslope','respZsave');
    end
    
end


%keyboard;

