
figopt1=1; saveg1=0; saveg2=0;
%%% Compare Sleep Epochs - Response to stimulation ###

%directoryname = '/data25/sjadhav/RippleInterruption/SJStimC_direct';
%prefix = 'sjc';
directoryname = '/data25/sjadhav/RippleInterruption/REb_direct';
prefix = 'REb';
days = [15];
tet=3;

e_pre = []; e_post = []; e_post2=[];
for d = 1:length(days)
    
    day = days(d)
    allepochs = [1 3 5];
    allepochs = 1; 
    
    %% Day n
    DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
    load(DIOfile);
    clr = {'k','r','g','c','b','m'};
    
    if figopt1==1
        figure; hold on;
        %redimscreen;
        redimscreen_land;
        orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','bold');
        set(0,'defaultaxeslinewidth',2);
    end
    
    %% Epoch 2 and 3 - Sleep 1 and 2
    
    
    for ep = 1:length(allepochs)
        
        e_stim = []; etmp_stim = [];
        epoch = allepochs(ep);
        
        stim = DIO{day}{epoch}{15};
        stim_starttime = stim.pulsetimes(:,1);
        stim_endtime = stim.pulsetimes(:,2);
        stim_length = stim.pulselength;
        stim_isi = stim.timesincelast(2:end);
        
      
        
        
        %% EEg filon tet 
        
        EEGfile = sprintf('%s/EEG/%seeg%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tet);
        load(EEGfile);
        e = eeg{day}{epoch}{tet};
        t = geteegtimes(e);
        pt = stim.pulsetimes ./ 10000;
        eind = lookup(pt(:,1), t);
        
%         EEGnostimfile = sprintf('%s/EEG/%seegnostim%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tet);
%         load(EEGnostimfile);
%         etmp = eeg{day}{epoch}{tet}.data;
        
        %% Ripple processed EEG 
        
        %     ripfile = sprintf('%s/EEG/%sripple%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tet);
        %     load(ripfile);
        %     ripamp = ripple{day}{epoch}{tet}.data(:,1);
        %     ripenv = ripple{day}{epoch}{tet}.data(:,3);
        %
        %     ripnostimfile = sprintf('%s/EEG/%sripplenostim%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tet);
        %     load(ripnostimfile);
        %     ripampfilt = ripple{day}{epoch}{tet}.data(:,1);
        %     ripenvfilt = ripple{day}{epoch}{tet}.data(:,3);
        
        
        for i =1:length(eind)-2;
            e_stim(i,:)=e.data(eind(i+1)-round(0.2*e.samprate):eind(i+1)+round(0.4*e.samprate));
            % etmp_stim(i,:)=etmp(eind(i+1)-round(0.2*e.samprate):eind(i+1)+round(0.4*e.samprate));
            %         ripamp_stim(i,:)=ripamp(eind(i+1)-0.5*e.samprate:eind(i+1)+1*e.samprate);
            %         ripampfilt_stim(i,:)=ripampfilt(eind(i+1)-0.5*e.samprate:eind(i+1)+1*e.samprate);
            %         ripenv_stim(i,:)=ripenv(eind(i+1)-0.5*e.samprate:eind6(i+1)+1*e.samprate);
            %         ripenvfilt_stim(i,:)=ripenvfilt(eind(i+1)-0.5*e.samprate:eind(i+1)+1*e.samprate);
        end
        
        taxis = [1:size(e_stim,2)]*1000/e.samprate;
        taxis = taxis - 0.2*1000;
        
        if figopt1==1
            plot(taxis, mean(e_stim),[clr{ep} '.-'],'Linewidth',2,'Markersize',6);
            uperr = mean(e_stim) + std(e_stim);
            lowerr = mean(e_stim) - std(e_stim);
            %line([0.2*1000:0.2*1000],[min(lowerr):max(uperr)],'Color','b','LineWidth',2,'LineStyle','--');
            jbfill(taxis,lowerr, uperr,[clr{ep}],[clr{ep}],1,0.2);
            
            if ep==1
                yplot = min(lowerr):1:max(uperr);
                xplot = 0*1000*ones(size(yplot));
                plot(xplot,yplot,'b--','Linewidth',2);
                title(['Day ' num2str(day) ': Evoked Response (Tet' num2str(tet) '): Pre-Sleep(black) and Post-sleep (red&green)'],...
                    'FontSize',24,'Fontweight','bold');
                xlabel('Time (ms)');
                ylabel('LFP (uV)');
            end
            
        end
        
        if ep==1
            e_pre = [e_pre; mean(e_stim)];
            %e_pre = e_pre(:,1:900);
        end
        if ep==2
            e_post = [e_post; mean(e_stim)];
            %e_post = e_post(:,1:900);
        end
        
        if ep==3
            e_post2 = [e_post2; mean(e_stim)];
            %e_post = e_post(:,1:900);
        end
        
        
        e_stim = []; etmp_stim = [];
        ripamp_stim = []; ripampfilt_stim = [];
        ripenv_stim = []; ripenvfilt_stim = [];
        
    end
    
    
    if saveg1==1,
        orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['Day' num2str(day) 'Tet' num2str(tet) '_Sleep1and2_StimResponse'],'jpg');
    end
    
    
end


%%%%%%%%%%%%
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

plot(taxis, mean(e_post2),['g.-'],'Linewidth',2,'Markersize',6);
uperr = mean(e_post2) + std(e_post2);
lowerr = mean(e_post2) - std(e_post2);
jbfill(taxis,lowerr, uperr,'g','g',1,0.2);

yplot = min(lowerr):1:max(uperr);
xplot = 0*1000*ones(size(yplot));
plot(xplot,yplot,'b--','Linewidth',2);

title(['Avg over Days ' num2str(min(days)) ':' num2str(max(days)) ': Evoked Response (Tet' num2str(tet) '): Pre-Sleep(black) and Post-sleep (red&green)'],...
    'FontSize',24,'Fontweight','bold');
xlabel('Time (ms)');
ylabel('LFP (uV)');



if saveg2==1,
    orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['Days' num2str(min(days)) ':' num2str(max(days)) '_Tet' num2str(tet) '_Sleep1and2_StimResponse'],'jpg');
    saveas(gcf,['Days' num2str(min(days)) ':' num2str(max(days)) '_Tet' num2str(tet) '_Sleep1and2_StimResponse'],'fig');
    
end


%%% Statistics on Peak Response

pre=[]; post1=[]; post2=[];

for i=1:length(days),
    
    pre(i) = min(e_pre(i,301:345));
    post1(i) = min(e_post(i,301:345));
    post2(i) = min(e_post2(i,301:345));
end
    
figure; hold on;
%redimscreen_land;
orient(gcf,'portrait'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    
bar(1,abs(mean(pre)),'k'); errorbar(1,abs(mean(pre)),abs(sem(pre)),'k');
bar(2,abs(mean(post1)),'r'); errorbar(2,abs(mean(post1)),abs(sem(post1)),'r');
bar(3,abs(mean(post2)),'g'); errorbar(3,abs(mean(post2)),abs(sem(post2)),'g');
set(gca,'XTick',[1 2 3],'XTickLabel',{'Pre','Post1','Post2'});

[h1,p1,ci1] = ttest2(pre, post1, 0.05, 'right'),
[h2,p2,ci2] = ttest2(pre, post2, 0.05, 'right'),

[hk1,pk1,cik1] = kstest2(pre, post1, 0.05, 'smaller'),
[hk2,pk2,cik2] = kstest2(pre, post2, 0.05, 'smaller'),

%%%%%%%%%%%%%%% %%%%%%%%


