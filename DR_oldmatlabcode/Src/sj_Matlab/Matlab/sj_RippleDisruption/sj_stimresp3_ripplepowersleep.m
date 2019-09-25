
saveg1=0; saveg2=0;
%%% Ripple ampl and power before and after removing stimulation artifact: Post-Sleep Epoch - Response to stimulation ###

directoryname = '/data25/sjadhav/RippleInterruption/SJStimC_direct';
prefix = 'sjc';
directoryname = '/data25/sjadhav/RippleInterruption/RE1_direct';
prefix = 'RE1';
days = [1 2];
tet=5;

e_amp = []; e_env = []; e_ampfilt = []; e_envfilt = [];
for d = 1:length(days)
        
    day = days(d)
    allepochs = 3;
    
    %% Day n
    DIOfile = sprintf('%s/%sDIO%02d.mat', directoryname, prefix, day);
    load(DIOfile);
    clr = {'k','r','g','c','b','m'};
    
    figure(10*d+1); hold on;
    %redimscreen;
    redimscreen_land;
    orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','bold');
    set(0,'defaultaxeslinewidth',2);
    
    figure(100*d+1); hold on;
    %redimscreen;
    redimscreen_land;
    orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','bold');
    set(0,'defaultaxeslinewidth',2);
    
    
    
    %% Epoch  3 - Sleep 2
    
    
    ep = 1;
    
    e_stim = []; etmp_stim = [];
    epoch = allepochs(ep);
    
    stim = DIO{day}{epoch}{15};
    stim_starttime = stim.pulsetimes(:,1);
    stim_endtime = stim.pulsetimes(:,2);
    stim_length = stim.pulselength;
    stim_isi = stim.timesincelast(2:end);
    
    %% Multiunit spikes on Tet 2, 5, 6
    multifile = sprintf('%s/%smulti%02d.mat', directoryname, prefix, day);
    load(multifile);
    %         multi2 = multi{day}{epoch}{2} ;
             multi5 = multi{day}{epoch}{5} ;
    %         multi6 = multi{day}{epoch}{6} ;
    
    
    %% EEg filon tet 6
    
    EEGfile = sprintf('%s/EEG/%seeg%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tet);
    load(EEGfile);
    e = eeg{day}{epoch}{tet};
    t = geteegtimes(e);
    pt = stim.pulsetimes ./ 10000;
    eind = lookup(pt(:,1), t);
    
    EEGnostimfile = sprintf('%s/EEG/%seegnostim%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tet);
    load(EEGnostimfile);
    etmp = eeg{day}{epoch}{tet}.data;
    
    %% Ripple processed EEG on 6
    
    ripfile = sprintf('%s/EEG/%sripple%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tet);
    load(ripfile);
    ripamp = ripple{day}{epoch}{tet}.data(:,1);
    ripenv = ripple{day}{epoch}{tet}.data(:,3);
    
    ripnostimfile = sprintf('%s/EEG/%sripplenostim%02d-%01d-%02d.mat', directoryname, prefix, day,epoch,tet);
    load(ripnostimfile);
    ripampfilt = ripple{day}{epoch}{tet}.data(:,1);
    ripenvfilt = ripple{day}{epoch}{tet}.data(:,3);
    
    
    for i =1:length(eind)-2;
        e_stim(i,:)=e.data(eind(i+1)-0.2*e.samprate:eind(i+1)+0.4*e.samprate);
        etmp_stim(i,:)=etmp(eind(i+1)-0.2*e.samprate:eind(i+1)+0.4*e.samprate);
        ripamp_stim(i,:)=double(ripamp(eind(i+1)-0.2*e.samprate:eind(i+1)+0.4*e.samprate));
        ripampfilt_stim(i,:)=double(ripampfilt(eind(i+1)-0.2*e.samprate:eind(i+1)+0.4*e.samprate));
        ripenv_stim(i,:)=double(ripenv(eind(i+1)-0.2*e.samprate:eind(i+1)+0.4*e.samprate));
        ripenvfilt_stim(i,:)=double(ripenvfilt(eind(i+1)-0.2*e.samprate:eind(i+1)+0.4*e.samprate));
    end
    
    taxis = [1:size(ripamp_stim,2)]*1000/e.samprate;
    taxis = taxis - 0.2*1000;
    
    figure(10*d+1); hold on;
    plot(taxis, mean(ripamp_stim),['k.-'],'Linewidth',2,'Markersize',6);
    uperr = mean(ripamp_stim) + std(ripamp_stim);
    lowerr = mean(ripamp_stim) - std(ripamp_stim);
    jbfill(taxis,lowerr, uperr,'k','k',1,0.2);
    yplot = min(lowerr):1:max(uperr);
    xplot = 0*1000*ones(size(yplot));
    
    plot(taxis, mean(ripampfilt_stim),['r.-'],'Linewidth',2,'Markersize',6);
    uperr = mean(ripampfilt_stim) + std(ripampfilt_stim);
    lowerr = mean(ripampfilt_stim) - std(ripampfilt_stim);
    jbfill(taxis,lowerr, uperr,'r','r',1,0.2);
    
    
    plot(xplot,yplot,'b--','Linewidth',2);
    title(['Day ' num2str(day) ': Ripple Ampl (Tet' num2str(tet) '): Pre-Filter(black) and Post-Filter (red)'],...
        'FontSize',24,'Fontweight','bold');
    xlabel('Time (ms)');
    ylabel('Ripple Amplitude');
    
    if saveg1==1,
        orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['Day' num2str(day) 'Tet' num2str(tet) '_Sleep2_RipplAmp'],'jpg');
    end

    
    
    figure(100*d+1); hold on;
       
    plot(taxis, mean(ripenv_stim),['k.-'],'Linewidth',2,'Markersize',6);
    uperr = mean(ripenv_stim) + std(ripenv_stim);
    lowerr = mean(ripenv_stim) - std(ripenv_stim);
    jbfill(taxis,lowerr, uperr,'k','k',1,0.2);
    yplot = min(lowerr):1:max(uperr);
    xplot = 0*1000*ones(size(yplot));
    
    plot(taxis, mean(ripenvfilt_stim),['r.-'],'Linewidth',2,'Markersize',6);
    uperr = mean(ripenvfilt_stim) + std(ripenvfilt_stim);
    lowerr = mean(ripenvfilt_stim) - std(ripenvfilt_stim);
    jbfill(taxis,lowerr, uperr,'r','r',1,0.2);
    
    
    plot(xplot,yplot,'b--','Linewidth',2);
    title(['Day ' num2str(day) ': Ripple Env(Tet' num2str(tet) '): Pre-Filter(black) and Post-Filter (red)'],...
        'FontSize',24,'Fontweight','bold');
    xlabel('Time (ms)');
    ylabel('Ripple Envelope');
    
     if saveg1==1,
        orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
        saveas(gcf,['Day' num2str(day) 'Tet' num2str(tet) '_Sleep2_RipplEnv'],'jpg');
    end

      
    e_amp = [e_amp; mean(ripamp_stim)];
    e_envfilt = [e_envfilt; mean(ripenvfilt_stim)];
    e_ampfilt = [e_ampfilt; mean(ripampfilt_stim)];
    e_env = [e_env; mean(ripenv_stim)];
    
    e_stim = []; etmp_stim = [];
    %ripamp_stim = []; ripampfilt_stim = [];
    %ripenv_stim = []; ripenvfilt_stim = [];
    
end







%%%%%%%%%%%%
figure; hold on;
redimscreen_land;
orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');

plot(taxis, mean(e_amp),['k.-'],'Linewidth',2,'Markersize',6);
uperr = mean(e_amp) + std(e_amp);
lowerr = mean(e_amp) - std(e_amp);
jbfill(taxis,lowerr, uperr,'k','k',1,0.2);
yplot = min(lowerr):1:max(uperr);
xplot = 0*1000*ones(size(yplot));

plot(taxis, mean(e_ampfilt),['r.-'],'Linewidth',2,'Markersize',6);
uperr = mean(e_ampfilt) + std(e_ampfilt);
lowerr = mean(e_ampfilt) - std(e_ampfilt);
jbfill(taxis,lowerr, uperr,'r','r',1,0.2);


plot(xplot,yplot,'b--','Linewidth',2);

title(['Avg over Days 4:7 : Ripple Ampl (Tet' num2str(tet) '): Pre-Filter(black) and Post-filter (red)'],...
    'FontSize',24,'Fontweight','bold');
xlabel('Time (ms)');
ylabel('Rip Ampl');


if saveg2==1,
    orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['Days4to7_Tet' num2str(tet) '_Sleep2_RipAmp'],'jpg');
    saveas(gcf,['Days4to7_Tet' num2str(tet) '_Sleep2_RipAmp'],'fig');
    
end


%%%%%%%%%%%%%%%5

figure; hold on;
redimscreen_land;
orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');

plot(taxis, mean(e_env),['k.-'],'Linewidth',2,'Markersize',6);
uperr = mean(e_env) + std(e_env);
lowerr = mean(e_env) - std(e_env);
jbfill(taxis,lowerr, uperr,'k','k',1,0.2);
yplot = min(lowerr):1:max(uperr);
xplot = 0*1000*ones(size(yplot));

plot(taxis, mean(e_envfilt),['r.-'],'Linewidth',2,'Markersize',6);
uperr = mean(e_envfilt) + std(e_envfilt);
lowerr = mean(e_envfilt) - std(e_envfilt);
jbfill(taxis,lowerr, uperr,'r','r',1,0.2);


plot(xplot,yplot,'b--','Linewidth',2);

title(['Avg over Days 4:7 : Ripple Env (Tet' num2str(tet) '): Pre-Filter(black) and Post-filter (red)'],...
    'FontSize',24,'Fontweight','bold');
xlabel('Time (ms)');
ylabel('Rip Power');



if saveg2==1,
    orient(gcf,'landscape'); hold on; set(gcf, 'PaperPositionMode', 'auto');
    saveas(gcf,['Days4to7_Tet' num2str(tet) '_Sleep2_RipEnv'],'jpg');
    saveas(gcf,['Days4to7_Tet' num2str(tet) '_Sleep2_RipEnv'],'fig');
    
end





%%%%%%%%%%%%%%% %%%%%%%%


