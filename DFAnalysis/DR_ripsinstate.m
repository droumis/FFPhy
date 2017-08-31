

%get ripples from areas and plot by behave/trial state

clear all
directory = '/mnt/data19/droumis/MEC_Expt/D05/D05_d04/testdata_ep04/';
% directory = '/mnt/data19/droumis/MEC_Expt/D05/D05_d04/testdata_ep02/';
% filename = 'D05_d04_e02_09-20-2015(12_27_15)'; %full rec name
filename = 'D05_d04_e04_09-20-2015(13_13_53)';
animalID = 'D05';
day = 4;
ep = 4;
centerportIN = 2;
centerportOUT = 5 +15;
portsIN = [1 2 3]; % need to use portins to assign performance
portsOUT = [4 5 6];
ports = [portsIN (portsOUT)+15]; %trodes dio assignment starts Douts at 16

linpos = loaddatastruct(directory, animalID, 'linpos', day);
load(sprintf('%s%sBehaveState%02d.mat',directory, animalID, day));
load(sprintf('%sdio/%sdio%02d.mat',directory, animalID, day));
%% get the ripples from each area
clear ripout
tetfilter =  '( isequal($area, ''ca1''))';
animalinfo = animaldef(animalID);
areas = [{'ca1'};{'mec'};{'por'}];
clear ripout
for i = 1:length(areas)
    ripout{i} = DR_kk_getconstimes(animalinfo{2},animalinfo{3},[day ep],sprintf('%sripplekons',areas{i}),1,'consensus_numtets', 3, 'maxvelocity',4);
    % kk_getconstimes(animaldir,animalprefix, epochs, eventconsname, tetfilter, varargin)
    % ripout = kk_getconstimes(animalinfo{2},animalinfo{3},[day ep], 'ripplesdons',1,'consensus_numtets',3);
    [periodtimes_rip2 periodripsinds] = dr_vec2list(ripout{i}{day}{ep}.cons,ripout{i}{day}{ep}.time);
    ripout{i}{day}{ep}.periodtimes_rip2 = periodtimes_rip2;
    ripout{i}{day}{ep}.periodripsinds = periodripsinds;
    plot(ripout{i}{day}{ep}.periodtimes_rip2', repmat([i], size(ripout{i}{day}{ep}.periodtimes_rip2')),'-', 'Color', rand(1,3), 'LineWidth',30); hold on;
end

%% plot lindist and overlay allthe rip times 
close all
time = linpos{day}{ep}.statematrix.time;
segmentnum = linpos{day}{ep}.statematrix.traj;
xvec = time;
lindist = linpos{day}{ep}.statematrix.linearDistanceToWells(:,1);
plot(xvec,lindist,'linewidth',1.5,'color',[.8 .8 .8]);
hold on
% clr = [.5 1 .8; .5 1 .8; .8 .9 .6; .8 .9 .6];
clr1 = [1 .9 .8; 1 .9 .8; .8 1 .8; .8 1 .8];
for segno = 1:4
    inds = [];
    inds = double(segmentnum == segno);
    inds(find(inds==0))=nan;
    % plot animal's travel lindist
    plot(xvec.*inds,lindist.*inds,'-','linewidth',2,'color',clr1(segno,:));
end

clr2 = ([.5 .5 1; .5 .7 .7; .3 .8 .3; .5 .5 1; .5 .7 .7; .3 .8 .3;]);
for i = 1:length(ports)
    if ports(i) == centerportIN || ports(i) == centerportOUT %is center
        yv = min(lindist) - 7; %plot triggers near bottom
    else
        yv = max(lindist) + 7 ; %plot triggers near top
    end
    if ports(i) > 15 %is a portOUT
        plot(dio{day}{ep}{ports(i)}.time(dio{day}{ep}{ports(i)}.state == 1), yv, '^', 'linewidth', 1, 'color', 'k'); hold on;
        plot(dio{day}{ep}{ports(i)}.time(dio{day}{ep}{ports(i)}.state == 0), yv, 'v', 'linewidth', 1, 'color', 'k'); hold on;
    else
        plot(dio{day}{ep}{ports(i)}.time(dio{day}{ep}{ports(i)}.state == 1), yv, '.', 'linewidth', 1, 'color', clr2(i,:)); hold on;
    end
    axis tight
end
fcclr = ([1 1 1 ; 0.2 .6 .9; 0 0 0]);
edgclr = ([1 0 0; 0.2 .6 .9; 0 0 0]);
lintyp = (['o' '*' 'x']);
linwid = ([2 3 1.5]);
marksz = ([7 4 7]);
for i = 1:length(areas)
    riptimesind = knnsearch(time, ripout{i}{day}{ep}.periodtimes_rip2(:,1));
    riptimesLIN = time(riptimesind);
    ripYLIN = lindist(riptimesind);
    plot(riptimesLIN, ripYLIN,lintyp(i), 'MarkerSize',marksz(i), 'LineWidth',linwid(i), 'MarkerEdgeColor',edgclr(i,:), 'MarkerFaceColor',fcclr(i,:)); hold on;%'MarkerFaceColor',[.49 1 .63] 'LineWidth',1,
end

    behavstorig = BehaveState{day}{ep}.statechangeseq;
    behavst = [behavstorig(1:end-1,2) behavstorig(2:end,:)]; %hack fix.. put the last last state transition time in front
    behavstOUTcorr = behavst(behavst(:,7) == 1 & behavst(:,9) == 1,[1,3]);
    behavstOUTinc = behavst(behavst(:,7) == 0 & behavst(:,9) == 1,[1,3]);
    behavstINcorr = behavst(behavst(:,7) == 1 & behavst(:,8) == 1,[1,3]);
    behavstINinc = behavst(behavst(:,7) == 0 & behavst(:,8) == 1,[1,3]);
    yl = ylim;
    xl = xlim;
for i = 1:length(behavstOUTcorr(:,1)); %patch correct outbound trials
    patch([behavstOUTcorr(i,1) behavstOUTcorr(i,2) behavstOUTcorr(i,2) behavstOUTcorr(i,1)],  [yl(2) yl(2) yl(1) yl(1)] , 'k', 'edgecolor','none', 'FaceAlpha',.1); 
end
for i = 1:length(behavstOUTinc(:,1)); %patch incorrect outbound trials
    patch([behavstOUTinc(i,1) behavstOUTinc(i,2) behavstOUTinc(i,2) behavstOUTinc(i,1)],  [yl(2) yl(2) yl(1) yl(1)] , 'b', 'edgecolor','none', 'FaceAlpha',.1); 
end
% set(gca,'Color',[0 0 0]);
%% collect the ripples that occur during either correct/incorrect inbound/outbound times using the behave state struct
clear ripOUTcorr ripOUTinc ripINcorr ripINinc ripbehav ripcountOUTcorr ripcountOUTinc ripbehav

ca1riptimestart = ripout{1}{day}{ep}.periodtimes_rip2(:,1); %make sure ca1 is areas(1)
ca1riptimeend = ripout{1}{day}{ep}.periodtimes_rip2(:,2); %make sure ca1 is areas(1)
ca1ripOUTcorr = logical(isExcluded(ca1riptimestart, behavstOUTcorr)); %all the ca1ripples that started during this behavior state
ca1ripOUTinc = logical(isExcluded(ca1riptimestart, behavstOUTinc));
ca1ripINcorr = logical(isExcluded(ca1riptimestart, behavstINcorr));
ca1ripINinc = logical(isExcluded(ca1riptimestart, behavstINinc));


for i = 1:length(areas)
    riptimestart = ripout{i}{day}{ep}.periodtimes_rip2(:,1);
    riptimeend = ripout{i}{day}{ep}.periodtimes_rip2(:,2);
    riptimes = ripout{i}{day}{ep}.periodtimes_rip2(:,:);
    ripOUTcorr{i} = riptimestart(logical(isExcluded(riptimestart, behavstOUTcorr)));
    ripOUTinc{i} = riptimestart(logical(isExcluded(riptimestart, behavstOUTinc)));
    ripINcorr{i} = riptimestart(logical(isExcluded(riptimestart, behavstINcorr)));
    ripINinc{i} = riptimestart(logical(isExcluded(riptimestart, behavstINinc)));
    
    for j = 1:length(behavstOUTcorr(:,1)); %cycle through every correct outbound trial period
        ripcountOUTcorr{i}(j) = sum(isExcluded(riptimestart, behavstOUTcorr(j,:))); %get the number of rips in each period
    end
    for j = 1:length(behavstOUTinc(:,1)); %cycle through every incorrect outbound trial period
        ripcountOUTinc{i}(j) = sum(isExcluded(riptimestart, behavstOUTinc(j,:))); %get the number of rips in each period
    end


    ripCA1toArea = logical(isExcluded(ca1riptimestart, riptimes)) | logical(isExcluded(ca1riptimeend, riptimes)); %if either the ca1ripstart or ca1ripend time overlaps with a rip period from this area
    

    ripbehav.ripCA1toArea(:,i) = ripCA1toArea;
%     ripbehav.OUTcorrincripCA1toArea(:,i) = OUTcorrincripCA1toArea;
    ripbehav.ripcountOUTcorr{i} = ripcountOUTcorr{i};
    ripbehav.ripcountOUTinc{i} = ripcountOUTinc{i};
    ripbehav.meannumOUTcorr(i) = length(ripOUTcorr{i})/length(behavstOUTcorr(:,1));
    ripbehav.meannumOUTinc(i) = length(ripOUTinc{i})/length(behavstOUTinc(:,1));
    ripbehav.meannumINcorr(i) = length(ripINcorr{i})/length(behavstINcorr(:,1));
    ripbehav.meannumINinc(i) = length(ripINinc{i})/length(behavstINinc(:,1));
    ripbehav.numOUTcorr(i) = length(behavstOUTcorr(:,1));
    ripbehav.numOUTinc(i) = length(behavstOUTinc(:,1));
    ripbehav.numINcorr(i) = length(behavstINcorr(:,1));
    ripbehav.numINinc(i) = length(behavstINinc(:,1));
    ripbehav.numripOUTcorr(i) = length(ripOUTcorr{i});
    ripbehav.numripOUTinc(i) = length(ripOUTinc{i});
    ripbehav.numripINcorr(i) = length(ripINcorr{i});
    ripbehav.numripINinc(i) = length(ripINinc{i});

    
    % get rippleamplitudes, smaller std ripples

end
    ripbehav.diffmeanOUTcorrinc = ripbehav.meannumOUTinc./ripbehav.meannumOUTcorr;
    % hmmm interesting.. there was a slightly larger increase to preceding incorrect than correct outbounds for POR than ca1/mec
    
    ripbehav.OUTcorrCA1MECratio = ripbehav.ripcountOUTcorr{1}./ripbehav.ripcountOUTcorr{2};
    ripbehav.OUTincCA1MECratio = ripbehav.ripcountOUTinc{1}./ripbehav.ripcountOUTinc{2};
    ripbehav.OUTcorrincCA1MECratiostats = [mean(ripbehav.OUTcorrCA1MECratio) std(ripbehav.OUTcorrCA1MECratio)/sqrt(length(ripbehav.OUTcorrCA1MECratio)) mean(ripbehav.OUTincCA1MECratio) std(ripbehav.OUTincCA1MECratio)/sqrt(length(ripbehav.OUTincCA1MECratio))];
    
    ripbehav.OUTcorrCA1PORratio = ripbehav.ripcountOUTcorr{1}./ripbehav.ripcountOUTcorr{3};
    ripbehav.OUTincCA1PORratio = ripbehav.ripcountOUTinc{1}./ripbehav.ripcountOUTinc{3};
    ripbehav.OUTcorrincCA1PORratiostats = [mean(ripbehav.OUTcorrCA1PORratio) std(ripbehav.OUTcorrCA1PORratio)/sqrt(length(ripbehav.OUTcorrCA1PORratio)) mean(ripbehav.OUTincCA1PORratio) std(ripbehav.OUTincCA1PORratio)/sqrt(length(ripbehav.OUTincCA1PORratio))];

    ripbehav.OUTcorrMECPORratio = ripbehav.ripcountOUTcorr{2}./ripbehav.ripcountOUTcorr{3};
    ripbehav.OUTincMECPORratio = ripbehav.ripcountOUTinc{2}./ripbehav.ripcountOUTinc{3};
    ripbehav.OUTcorrincMECPORratiostats = [mean(ripbehav.OUTcorrMECPORratio) std(ripbehav.OUTcorrMECPORratio)/sqrt(length(ripbehav.OUTcorrMECPORratio)) mean(ripbehav.OUTincMECPORratio) std(ripbehav.OUTincMECPORratio)/sqrt(length(ripbehav.OUTincMECPORratio))];
    %% bar plotting
    figure(1)
    bar(1:2,[ripbehav.OUTcorrincMECPORratiostats(1) ripbehav.OUTcorrincMECPORratiostats(3)]); hold on;
    errorbar(1:2,[ripbehav.OUTcorrincMECPORratiostats(1) ripbehav.OUTcorrincMECPORratiostats(3)],[ripbehav.OUTcorrincMECPORratiostats(2) ripbehav.OUTcorrincMECPORratiostats(4)],'.')
    
    figure(2)
    bar(1:2,[ripbehav.OUTcorrincCA1PORratiostats(1) ripbehav.OUTcorrincCA1PORratiostats(3)]); hold on;
    errorbar(1:2,[ripbehav.OUTcorrincCA1PORratiostats(1) ripbehav.OUTcorrincCA1PORratiostats(3)],[ripbehav.OUTcorrincCA1PORratiostats(2) ripbehav.OUTcorrincCA1PORratiostats(4)],'.')
    
    figure(3)
    bar(1:2,[ripbehav.OUTcorrincCA1MECratiostats(1) ripbehav.OUTcorrincCA1MECratiostats(3)]); hold on;
    errorbar(1:2,[ripbehav.OUTcorrincCA1MECratiostats(1) ripbehav.OUTcorrincCA1MECratiostats(3)],[ripbehav.OUTcorrincCA1MECratiostats(2) ripbehav.OUTcorrincCA1MECratiostats(4)],'.')
    
    %%
    
%     ripcountOUTcorr = [ca1ripOUTcorr ripbehav.ripCA1toArea];
    OUTcorrincripCA1toArea = [ca1ripOUTcorr ca1ripOUTinc ripbehav.ripCA1toArea];
    corrchange = diff(OUTcorrincripCA1toArea(:,1));
    try
        corrtrialind = [find(corrchange == 1)+1 find(corrchange == -1)];
    catch %if the epoch ends before a correct trial is completed
        booo = find(corrchange == 1);
        boo = booo(1:end-1); %truncate the trouble maker!
        corrtrialind = [boo+1 find(corrchange == -1)];
    end
    incchange = diff(OUTcorrincripCA1toArea(:,2));
    inctrialind = [find(incchange == 1)+1 find(incchange == -1)];
    
    for i = 1:length(corrtrialind(:,1))
        corrtrialripinfo(i,1) = (corrtrialind(i,2) - corrtrialind(i,1))+1; %num of ripples in this trial
        corrtrialripinfo(i,2) = sum((sum(OUTcorrincripCA1toArea(corrtrialind(i,1):corrtrialind(i,2),3:5),2) == 3),1); %num of ripples in ca1mecpor     
        corrtrialripinfo(i,3) = sum((sum(OUTcorrincripCA1toArea(corrtrialind(i,1):corrtrialind(i,2),[3 4]),2) == 2),1); %num of ripples in ca1mec     
        corrtrialripinfo(i,4) = sum((sum(OUTcorrincripCA1toArea(corrtrialind(i,1):corrtrialind(i,2),[3 5]),2) == 2),1); %num of ripples in ca1por     
        corrtrialripinfo(i,5) = sum((sum(OUTcorrincripCA1toArea(corrtrialind(i,1):corrtrialind(i,2),[4 5]),2) == 2),1); %num of ripples in ca1por 
        
        corrtrialripstats(i,1) = corrtrialripinfo(i,2)/corrtrialripinfo(i,1);
        corrtrialripstats(i,2) = corrtrialripinfo(i,3)/corrtrialripinfo(i,1);
        corrtrialripstats(i,3) = corrtrialripinfo(i,4)/corrtrialripinfo(i,1);
        corrtrialripstats(i,4) = corrtrialripinfo(i,5)/corrtrialripinfo(i,1);
    end
    
    for i = 1:length(inctrialind(:,1))
        inctrialripinfo(i,1) = (inctrialind(i,2) - inctrialind(i,1))+1; %num of ripples in this trial
        inctrialripinfo(i,2) = sum((sum(OUTcorrincripCA1toArea(inctrialind(i,1):inctrialind(i,2),3:5),2) == 3),1); %num of ripples in ca1mecpor     
        inctrialripinfo(i,3) = sum((sum(OUTcorrincripCA1toArea(inctrialind(i,1):inctrialind(i,2),[3 4]),2) == 2),1); %num of ripples in ca1mec     
        inctrialripinfo(i,4) = sum((sum(OUTcorrincripCA1toArea(inctrialind(i,1):inctrialind(i,2),[3 5]),2) == 2),1); %num of ripples in ca1por     
        inctrialripinfo(i,5) = sum((sum(OUTcorrincripCA1toArea(inctrialind(i,1):inctrialind(i,2),[4 5]),2) == 2),1); %num of ripples in ca1por
        
        inctrialripstats(i,1) = inctrialripinfo(i,2)/inctrialripinfo(i,1);
        inctrialripstats(i,2) = inctrialripinfo(i,3)/inctrialripinfo(i,1);
        inctrialripstats(i,3) = inctrialripinfo(i,4)/inctrialripinfo(i,1);
        inctrialripstats(i,4) = inctrialripinfo(i,5)/inctrialripinfo(i,1);
    end
    
    meancorrtrialripstats = mean(corrtrialripstats);
    stderrcorrtrialripstats = std(corrtrialripstats)/sqrt(length(corrtrialripstats(:,1)));
    stdcorrtrialripstats = std(corrtrialripstats);
    meaninctrialripstats = mean(inctrialripstats);
    stderrinctrialripstats = std(inctrialripstats)/sqrt(length(inctrialripstats(:,1)));
    stdinctrialripstats = std(inctrialripstats);
    
    figure(1)
    h = bar([meancorrtrialripstats(2)' meaninctrialripstats(2)'],'k'); hold on
%     errorbar([.8 1.2; 1.8 2.2; 2.8 3.2; 3.8 4.2], [meancorrtrialripstats' meaninctrialripstats'], [stderrcorrtrialripstats' stderrinctrialripstats'], '.');
        errorbar([1 2], [meancorrtrialripstats(2)' meaninctrialripstats(2)'], [stderrcorrtrialripstats(2)' stderrinctrialripstats(2)'], '.');
        
       [h,p,ci,stats] = ttest2(corrtrialripstats(:,2), inctrialripstats(:,2));
       
        B = {'Correct';'Incorrect'};
       set(gca,'XTickLabel',B, 'fontweight','bold','fontsize',9);
       title('norm joint CA1-MEC ripple rate by trial type', 'fontweight','bold','fontsize',12)
%% save data so that i can combine epochs by hand
       if 1 
           ripsinbehavstatedata.corrtrialripinfo = corrtrialripinfo;
           ripsinbehavstatedata.corrtrialripstats = corrtrialripstats;
           ripsinbehavstatedata.inctrialripinfo = inctrialripinfo;
           ripsinbehavstatedata.inctrialripstats = inctrialripstats;
           save(sprintf('%sdio/%sripsinbehavstatedata%02d.mat',directory, animalID, day), 'ripsinbehavstatedata');
       end
       %% combine epochs 2,4
       dir1 = '/mnt/data19/droumis/MEC_Expt/D05/D05_d04/testdata_ep04/';
       dir2 = '/mnt/data19/droumis/MEC_Expt/D05/D05_d04/testdata_ep02/';
       data4 = load(sprintf('%sdio/%sripsinbehavstatedata%02d.mat',dir1, animalID, day));
       data2 = load(sprintf('%sdio/%sripsinbehavstatedata%02d.mat',dir2, animalID, day));
       
       if 0 %combine
       corrtrialripstats = [data2.ripsinbehavstatedata.corrtrialripstats; data4.ripsinbehavstatedata.corrtrialripstats];
       corrtrialripinfo = [data2.ripsinbehavstatedata.corrtrialripinfo; data4.ripsinbehavstatedata.corrtrialripinfo];
       inctrialripstats = [data2.ripsinbehavstatedata.inctrialripstats; data4.ripsinbehavstatedata.inctrialripstats];
%        inctrialripstats = [data2.ripsinbehavstatedata.inctrialripstats; data4.ripsinbehavstatedata.inctrialripstats(1:end-1,:)];
       inctrialripinfo = [data2.ripsinbehavstatedata.inctrialripinfo; data4.ripsinbehavstatedata.inctrialripinfo];
       
       end 
       
       if 1 %don't combine
           corrtrialripstats = [data2.ripsinbehavstatedata.corrtrialripstats];
           inctrialripstats = [data2.ripsinbehavstatedata.inctrialripstats];
       end
       
       meancorrtrialripstats = mean(corrtrialripstats);
       stderrcorrtrialripstats = std(corrtrialripstats)/sqrt(length(corrtrialripstats(:,1)));
       stdcorrtrialripstats = std(corrtrialripstats);
       meaninctrialripstats = mean(inctrialripstats);
       stderrinctrialripstats = std(inctrialripstats)/sqrt(length(inctrialripstats(:,1)));
       stdinctrialripstats = std(inctrialripstats);
       
       figure(1)
       h = bar([meancorrtrialripstats(2)' meaninctrialripstats(2)'],'k'); hold on
%     errorbar([.8 1.2; 1.8 2.2; 2.8 3.2; 3.8 4.2], [meancorrtrialripstats' meaninctrialripstats'], [stderrcorrtrialripstats' stderrinctrialripstats'], '.');
        errorbar([1 2], [meancorrtrialripstats(2)' meaninctrialripstats(2)'], [stderrcorrtrialripstats(2)' stderrinctrialripstats(2)'], '.');
        
       [h,p,ci,stats] = ttest2(corrtrialripstats(:,2), inctrialripstats(:,2));
       
        B = {'Correct';'Incorrect'};
       set(gca,'XTickLabel',B, 'fontweight','bold','fontsize',9);
       title('norm joint CA1-MEC ripple rate by trial type', 'fontweight','bold','fontsize',12)


%% plot lindist and overlay the rip times coded by behav state





%% plot histograms of rips by behav state


    %try lowering/raising ripple std. 
    %try to get amplitude
    %1pm start processing a sleep session for figure 2    

    %then process an early day and a later day epoch from the single W track.
    
    %then process a sleep epoch from day 4 to get the example trace of up/down states with ripples
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
