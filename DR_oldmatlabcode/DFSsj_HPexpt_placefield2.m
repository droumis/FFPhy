
% Changing how I plot openfield rate from version 1. I call twodoccupanct instead of openfieldrate
% Open field is now separated by trajectory and a separate plot is made.

% Plotting of Open Field Rate and Linearized
% Place Fields if desired

clear;
runscript = 1;
savedata = 1; % save data option - only works if runscript is also on
figopt1 = 1; % Figure Options
savelinfields = 1; % To save trajdata in [prefix-linfields-day] file for each day
                   % Also mapdata in [prefix-mapfields-day] file for day
                   % Save mapdata as opnfield rate, as well as separate trajectories

savedir = '/data25/sjadhav/HPExpt/ProcessedDataDRtemp/';
savefile = [savedir 'HPc_PFCfields']; % HPa and HPb fields  %%DR
%savefile = [savedir 'HPa_PlaceFieldsTrajMaps']; %Day1-4
%savefile = [savedir 'HPb_PlaceFieldsTrajMaps']; %Day1-4
%savefile = [savedir 'PlaceFieldsTrajMaps_Day8'];

minabsvel = 3;  % cm/sec - Most conservative for runs and place fields
minlinvel = 5;

plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days
plotanimidx = 1; % To pick animals for plotting


% If runscript, run Datafilter and save data
if runscript == 1
    
    
    %Animal selection
    %-----------------------------------------------------
    animals = {'HPc'};  %%DR
    
    %Filter creation
    %--------------------------------------------------------
    
    % Epoch filter
    % -------------
    dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    %epochfilter = 'isequal($type, ''run'')';
    epochfilter = 'isequal($environment, ''wtr1'')';

    
    % Cell filter
    % -----------
    placecellfilter = '( strcmp($tag, ''PFC'') && ($numspikes > 100))';
    %placecellfilter = '( strcmp($tag, ''CA1Pyr'') || strcmp($tag, ''iCA1Pyr'') || strcmp($tag, ''PFC'') )';
    %placecellfilter = '( strcmp($tag2, ''CA1Pyr'') || strcmp($tag2, ''iCA1Pyr'') || strcmp($tag2, ''PFC'') ) && ($numspikes > 100)';

    
    % Time filter
    % -----------
    
    % Either use tetlist for riple detection for ripfilter,
    % or use ripple tet filter based on tag in tetinfo marking elec as ripple tet
 
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    % Abs linear velocity(thrs 5)/ velocity(thrs 3) time filter
    %timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getvelpos', '(($absvel >= 3))'}, {'DFTFsj_getriptimes','($nripples == 0)',tetlist,'minthresh',2} }
    %timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getvelpos', '(($absvel >= 3))'}, {'DFTFsj_getriptimes','($nripples == 0)',[],'tetfilter',riptetfilter,'minthresh',2} }
    
    % Linear velocity time filter implemented in getlinstate
    %timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getriptimes','($nripples == 0)',tetlist,'minthresh',2} };
    
    timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 3))', 6}, {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };
    %timefilter = { {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',3} };
    
    % Note - For includestates=1:6, state will never be -1
    
    
    % Iterator
    % --------
    iterator = 'singlecellanal';
    
    % Filter creation
    % ----------------
    psf = createfilter('animal',animals,'days',dayfilter,'epochs',epochfilter,'cells',placecellfilter,'excludetime', timefilter,'iterator', iterator);
    % psf = createfilter('animal',animals,'days',dayfilter,'epochs',epochfilter,'cells',placecellfilter,'iterator', iterator);
    
    % Set analysis function
    % ----------------------
    psf = setfilterfunction(psf, 'DFAsj_filtercalclinfields_tf', {'spikes', 'linpos'}, 'binsize', 2);
    %psf = setfilterfunction(psf, 'DFAsj_filtercalclinfields', {'spikes', 'linpos', 'pos'}, 'binsize', 2, 'minabsvel', 3);
    pmf = setfilterfunction(psf, 'DFAsj_openfieldrate_tf', {'spikes', 'linpos', 'pos'}, 'binsize', 1, 'std', 2);
    pof = setfilterfunction(psf, 'DFAsj_twodoccupancy', {'spikes', 'linpos','pos'}, 'binsize', 1, 'std', 2);
    
    
   disp('Finished filter creation');
    
    % Run analysis
    % ------------
    psf = runfilter(psf);  % Place Field Stability
    pmf = runfilter(pmf);  % Place Field Map
    %pof = runfilter(pof);  % Place Field Map - Separate trajectories
    
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear figopt1 runscript plotdays plotanimidx savedata savelinfields
        save(savefile);
    end
    
else
    
    load(savefile);
    
end  % end runscript

if ~exist('savedata')
    return
end

% ----------------------------------

gatherdata = 1; savegatherdata = 1;
gatherdatafile = [savedir 'HPa_PFCfields_gather']; 


% To Control Plotting, enter parameters here

if ~isempty(plotanimidx)
    useanim = plotanimidx;
else
    useanim = 1:length(psf); % Use all animals
end
if ~isempty(plotdays)
    usedays = plotdays;
else
    usedays = [];   % Get for each animal separately
end

% if savedata == 1
%     figopt1 = 0; plotdays = [];
% end
% if exist('runscript') & exist('figopt1')
%     if (runscript == 0) && (figopt1==1)
%         days = plotdays;
%     end
% end


% ---------------------------------


% Get trajdata and days and epochs
trajdata = []; index=[];

allanimindex=[]; allmaps=[]; alltrajs=[];
for anidx = 1:length(useanim)
    an = useanim(anidx);
    for i=1:length(psf(an).output{1}),
        index{an}(i,:)=psf(an).output{1}(i).index;
        alltrajdata{an}{i}=psf(an).output{1}(i).trajdata;
        allmapdata{an}{i}=pmf(an).output{1}(i);
        %allmapdata_sep{an}{i}=pof(an).output{1}(i);
        
        % Only indexes
        animindex=[an modf(an).output{1}(i).index]; % Put animal index in front
        allanimindex = [allanimindex; animindex]; % Collect all Anim Day Epoch Tet Cell Index
        allmaps{i} = pmf(an).output{1}(i);
        alltrajs{i} = psf(an).output{1}(i).trajdata;
    end
end

if savegatherdata == 1
        save(gatherdatafile);
end

% Initialize
clr = {'b','r','g','m','c','y','k','r'};
novel=1:3;
fam=4:8;


for anidx = 1:length(useanim)
    an = useanim(anidx);
    prefix = psf(an).animal{1};
    animdirect = psf(an).animal{2};
    if (animdirect(end) == '/')
        animdirect = animdirect(1:end-1);
    end
    % Get days and epochs for current animal
    if isempty(usedays) 
        days = unique(psf(an).epochs{1}(:,1));
    else
        days = usedays;
    end
    allepochs = unique(psf(an).epochs{1}(:,2));
    %days = unique(index(:,1));
    %allepochs = unique(index(:,2));
    
    %Initialize
    cnt_allneu=0;  % Count neurons across days
    %%%allcorr=[];
    %%%allcorr_day=[];
    cntneu=0; % These are also across days
    cntneu_ep2=0; % These are also across days
    
    for d = 1:length(days)
        
        % Reset linfields in case you are going to save it
        linfields = [];
        mapfields = [];
        
        day = days(d)
        cnt_dayneu = 0;
        
        if day==1, 
            allepochs = [4 6];
        else
            allepochs = [2 4];
        end
        
        
        
        dayidxs = find(index{an}(:,1)==day);
        % Get tetlist
        tet = unique(index{an}(dayidxs,3));
        
        for elec = 1:length(tet)
            currtet = tet(elec);
            % Get celllist
            daytetidxs = find(index{an}(:,1)==day & index{an}(:,3)==currtet);
            cells = unique(index{an}(daytetidxs,4));
            
            for neuron = 1:length(cells)
                currcell = cells(neuron);
                cnt_dayneu = cnt_dayneu + 1;
                cnt_allneu = cnt_allneu + 1;  % Neuron Count across days - Not reset
                
                neuidx{an}{day}(cnt_dayneu) = cnt_allneu; % Very Imp - To index into place cell data later
                
                for ep = 1:length(allepochs)
                    epoch = allepochs(ep);
                    curridx = find( index{an}(:,1)==day & index{an}(:,2)==epoch & index{an}(:,3)==currtet & index{an}(:,4)==currcell);
                    
                    if ~isempty(curridx), 
                        trajdata = alltrajdata{an}{curridx};
                        mapdata = allmapdata{an}{curridx}.smoothedspikerate;
                        %mapdata_sep = allmapdata_sep{an}{curridx}.smoothedspikerate; 
                
                     maxrate=[];
                    % Save in linfields in case you are going to save it
                    linfields{day}{epoch}{currtet}{currcell}=trajdata;
                    mapfields{day}{epoch}{currtet}{currcell}=allmapdata{an}{curridx}; % Save everything .smoothedspikerate is what you want
                    %mapfields_sep{day}{epoch}{currtet}{currcell}=allmapdata_sep{an}{curridx}; % .smoothedspikerate will have 4 cell arrays                     
                    % Get traj with maximum peak fir rate
                    for i=1:length(trajdata),
                        maxrate(i) = max(trajdata{i}(5:end-5,5));
                    end
                    
                    if ~exist('figopt1'), figopt1=0; end
                    if figopt1==1,
                        
                        % Linearized traj
                        figure(str2num([num2str(day) num2str(currtet) num2str(currcell) num2str(0)]));
                        %redimscreen_2horsubplots;
                        redimscreen;
                        hold on;
                        subplot(2,length(allepochs),ep); hold on;
                        for i=1:length(trajdata),
                            plot(trajdata{i}(:,5),[clr{i} '.-'],'Linewidth',2);
                            %maxrate(i) = max(trajdata{i}(5:end-5,5));
                        end
                        if ep==1
                            title(['Day ' num2str(day) '  Epoch ' num2str(epoch) '  Tet ' num2str(currtet) '  Cell ' num2str(currcell)],'FontSize',14,'Fontweight','bold');
                            %legend('InRealRight','OutRealLeft','OutRealRight','InRealLeft');
                            legend('OutRealLeft','InRealLeft','OutRealRight','InRealRight');
                            xlabel ('Position along linear trajectory (cm)','FontSize',14,'Fontweight','bold');
                            ylabel ('Firing Rate (Hz)','FontSize',14,'Fontweight','bold');
                            
                        else
                            title(['Epoch ' num2str(epoch)],'FontSize',14,'Fontweight','bold');
                            xlabel ('Position along linear trajectory (cm)','FontSize',14,'Fontweight','bold');
                        end
                        
                        % Maps -
                        % 0) Entire Map
                        subplot(2,length(allepochs),length(allepochs)+ep); hold on;
                        imagesc(flipud(mapdata)); colorbar
                        set(gca,'YLim',[0 110]);
                        set(gca,'XLim',[0 110]);
                        title(['Epoch ' num2str(epoch)], 'FontSize',24,'Fontweight','bold');
                        if ep==1
                            xlabel ('X-position (cm)','FontSize',24,'Fontweight','bold');
                            ylabel ('Y-position (cm)','FontSize',24,'Fontweight','bold');
                        else
                            xlabel ('X-position (cm)','FontSize',24,'Fontweight','bold');
                        end
                        
                        % 1) Each traj spearately - 1 fig for each epoch
%                                             figure(str2num([num2str(day) num2str(currtet) num2str(currcell) num2str(ep)]));
%                                             titlestr = ['Out RealLeft';'In  RealLeft';'OutRealRight';'In RealRight'];
%                                             for tr = 1:length(mapdata_sep)
%                                                 subplot(2,length(allepochs),tr); hold on;
%                                                 imagesc(flipud(mapdata_sep{tr})); colorbar
%                                                 title(['Epoch' num2str(epoch) ': ' titlestr(tr,:)], 'FontSize',14,'Fontweight','bold');
%                                                 if ep==1
%                                                     xlabel ('X-position (cm)','FontSize',14,'Fontweight','bold');
%                                                     ylabel ('Y-position (cm)','FontSize',14,'Fontweight','bold');
%                                                 else
%                                                     xlabel ('X-position (cm)','FontSize',14,'Fontweight','bold');
%                                                 end
%                                             end
                        
                        % 2) Combine Outs and Ins / Or alternatively, Rights and Lefts. Just 1 fig
                        
                        % Rightmap
%                         if isempty(mapdata_sep{3}),
%                             rightmap = mapdata_sep{4};
%                         elseif isempty(mapdata_sep{4}),
%                             rightmap = mapdata_sep{3};
%                         else
%                             a=min(size(mapdata_sep{3},1),size(mapdata_sep{4},1));
%                             b=min(size(mapdata_sep{3},2),size(mapdata_sep{4},2));
%                             rightmap = (mapdata_sep{3}(1:a,1:b) + mapdata_sep{4}(1:a,1:b))./2;
%                         end
%                         
%                         % Leftmap
%                         if isempty(mapdata_sep{1}),
%                             leftmap = mapdata_sep{2};
%                         elseif isempty(mapdata_sep{2}),
%                             leftmap = mapdata_sep{1};
%                         else
%                             a=min(size(mapdata_sep{1},1),size(mapdata_sep{2},1));
%                             b=min(size(mapdata_sep{1},2),size(mapdata_sep{2},2));
%                             leftmap = (mapdata_sep{1}(1:a,1:b) + mapdata_sep{2}(1:a,1:b))./2;
%                         end
%                         
%                         figure(str2num([num2str(day) num2str(currtet) num2str(currcell) num2str(1)]));
%                         redimscreen
%                         subplot(2,length(allepochs),2*ep-1); hold on;
%                         imagesc(flipud(rightmap)); colorbar
%                         title(['Epoch' num2str(epoch) ': Rightmap'], 'FontSize',14,'Fontweight','bold');
%                         xlabel ('X-position (cm)','FontSize',14,'Fontweight','bold');
%                         ylabel ('Y-position (cm)','FontSize',14,'Fontweight','bold');
%                         subplot(2,length(allepochs),2*ep); hold on;
%                         imagesc(flipud(leftmap)); colorbar
%                         title(['Epoch' num2str(epoch) ': Leftmap'], 'FontSize',14,'Fontweight','bold');
%                         xlabel ('X-position (cm)','FontSize',14,'Fontweight','bold');
%                         ylabel ('Y-position (cm)','FontSize',14,'Fontweight','bold');
                        
                        
                    end % end figopt
                    
                    if ep==1
                        cntneu=cntneu+1;
                        allplacemaps_ep1{an}{cntneu} = allmapdata{an}{curridx};
                        alllintrajs_ep1{an}{cntneu} = trajdata;
                        allpeakrates_ep1{an}{cntneu} = maxrate;
                    end
                    if ep==2
                        cntneu_ep2=cntneu_ep2+1;
                        allplacemaps_ep2{an}{cntneu_ep2} = allmapdata{an}{curridx};
                        alllintrajs_ep2{an}{cntneu_ep2} = trajdata;
                        allpeakrates_ep2{an}{cntneu_ep2} = maxrate;
                    end
                    
                    trajplace{an}{day}{ep}{currtet}{currcell}=trajdata;
                    mapplace{an}{day}{ep}{currtet}{currcell}=allmapdata{an}{curridx};
                    %mapplace_sep{an}{day}{ep}{currtet}{currcell}=allmapdata_sep{an}{curridx};
                    
                    
                    else % if there is no cell corresponding to these indices
                        day, epoch, currtet, currcell
                        keyboard;                 
                        trajdata=[];
                        mapdata=[];
                    
                    end % if ~ismpty curridx
                    
                end % end epoch
                
                if figopt1==1,
                    keyboard;  % Pause and return control after neurons plot - 2 or 3 figs in all
                end
                
            end % end neuron
            
        end % end tet
        
        
        % FOR DAY
        % Get correlation data for current day neurons
%        for nc=1:cnt_dayneu
%             n = neuidx{an}{day}(nc);
%             
%             % Find which trajectory had the max firing rate in all epochs
%             peakrate=[]; peaktraj=[];
%             [peakrate,peaktraj1] = max(allpeakrates_ep1{an}{n});
%             [peakrate,peaktraj2] = max(allpeakrates_ep2{an}{n});
%             
%             if peaktraj1(1)==peaktraj2(1),
%                 allmaxtraj{an}{n}=peaktraj1(1);
%             else
%                 allmaxtraj{an}{n}=[peaktraj1(1) peaktraj2(1)];
%             end
%             
%             % Get correlations between place fields
%             single_traj=0; double_traj=0;
%             
%             %         % 2d corr
%             sizex = min([length(allplacemaps_ep1{an}{n}.yticks),length(allplacemaps_ep2{an}{n}.yticks)]);
%             sizey = min([length(allplacemaps_ep1{an}{n}.xticks),length(allplacemaps_ep2{an}{n}.xticks)]);
% %            corr2d = corr2(allplacemaps_ep1{an}{n}.smoothedspikerate(1:sizex,1:sizey), allplacemaps_ep2{an}{n}.smoothedspikerate(1:sizex,1:sizey));
%             
%             % Method 0: Get correlation across all trajectories
%             place_ep1all=[]; place_ep2all=[]; % Reset for each neuron
%             for tr=1:4
%                 place_ep1 = alllintrajs_ep1{an}{n}{tr}(:,5);
%                 place_ep2 = alllintrajs_ep2{an}{n}{tr}(:,5);
%                 le = min(length(place_ep1),length(place_ep2));
%                 place_ep1 = place_ep1(1:le);
%                 place_ep2 = place_ep2(1:le);
%                 
%                 invalid1 = find(isnan(place_ep1)==1);
%                 invalid2 = find(isnan(place_ep2)==1);
%                 invalid=unique([invalid1; invalid2]);
%                 place_ep1(invalid)=[];
%                 place_ep2(invalid)=[];
%                 
%                 place_ep1all = [place_ep1all; place_ep1];
%                 place_ep2all = [place_ep2all; place_ep2];
%             end              
%             currcorr = corrcoef(place_ep1all,place_ep2all);
%             currcorr = currcorr(1,2);
%             
%             allcorr{an}(n)=currcorr;          
%             allcorr_day{an}{day}(nc)=currcorr;
%             
%         end % end cnt_dayneu
        
        
         % Save linfields and mapfields for day if asked for
        %savefile = [animdirect,prefix,'linfields',num2str(day)];
        if savelinfields==1
            savefile = sprintf('%s/%slinfields%02d.mat', animdirect, prefix, day);
            save(savefile,'linfields');
            savefile = sprintf('%s/%smapfields%02d.mat', animdirect, prefix, day);
            save(savefile,'mapfields');
            %savefile = sprintf('%s/%smapfields_sep%02d.mat', animdirect, prefix, day);
            %save(savefile,'mapfields_sep');
        end
     
    end % end day
    
    if savelinfields==1
        disp(['Saved linfields and mapfields for animal ',prefix]);
    end
    
end % end animal



allcorrvec=[]; 
for anidx = 1:length(useanim)
    an = useanim(anidx);
    allcorrvec = [allcorrvec; allcorr{an}(:)];
end

disp('Mean LinTrajCorr');
meancorr = nanmean(allcorrvec),
errcorr = nansem(allcorrvec),

% disp('Mean 2dCorr');
% nanmean(all2dcorr),
% nansem(all2dcorr),

%% Summarize correlation - 1) Bar plot for all, novel, familiar days, & 2) Corrln across days

set(0,'defaultaxesfontsize',24);set(0,'defaultaxesfontweight','normal');
set(0,'defaultaxeslinewidth',2);

% 1) Bar plot for all, novel, familiar days
figure; hold on; redimscreen_figforppt1;
bar(1,meancorr,'k'); errorbar(1,meancorr,errcorr,'k','LineWidth',3);
set(gca,'XTick',[1],'XTickLabel',{'All Days'},'FontSize',24,'Fontweight','normal');
axis([0 2 0 1.0])

novcorrvec=[]; 
for anidx = 1:length(useanim)
    an = useanim(anidx);
    if isempty(usedays) 
        days = unique(psf(an).epochs{1}(:,1));
    else
        days = usedays;
    end
    if ~isempty(intersect(days,novel))
        noveldays = intersect(days,novel);
        novcorr{an} = [];        
        for nd=noveldays
            novcorr{an} = [novcorr{an}, allcorr_day{an}{nd}];
        end  
        
        novcorrvec = [novcorrvec; novcorr{an}(:)]; 
    end   
end
if ~isempty(novcorrvec)
    bar(2,nanmean(novcorrvec),'r'); errorbar(2,nanmean(novcorrvec),nansem(novcorrvec),'r','LineWidth',3);
    set(gca,'XTick',[1 2],'XTickLabel',{'All Days','Novel Days'},'FontSize',24,'Fontweight','normal');
    axis([0 3 0 1.0])
end


famcorrvec=[];
for anidx = 1:length(useanim)
    an = useanim(anidx);
    if isempty(usedays) 
        days = unique(psf(an).epochs{1}(:,1));
    else
        days = usedays;
    end
    if ~isempty(intersect(days,fam))
        famdays = intersect(days,fam);
        famcorr{an} = [];        
        for fd=famdays
            famcorr{an} = [famcorr{an}, allcorr_day{an}{fd}];
        end
        
        famcorrvec = [famcorrvec; famcorr{an}(:)];
    end
end
    
    if isempty(novcorrvec)
        bar(2,nanmean(famcorrvec),'g'); errorbar(2,nanmean(famcorrvec),nansem(famcorrvec),'g','LineWidth',3);
        set(gca,'XTick',[1 2],'XTickLabel',{'All Days','Fam Days'},'FontSize',24,'Fontweight','normal');
        axis([0 3 0 1.0])
    else
        
        bar(3,nanmean(famcorrvec),'g'); errorbar(3,nanmean(famcorrvec),nansem(famcorrvec),'g','LineWidth',3);
        set(gca,'XTick',[1 2 3],'XTickLabel',{'All Days','Novel Days','Fam Days'},'FontSize',24,'Fontweight','normal');
        axis([0 4 0 1.0])
    end

title(['Place Field Stability' ],'FontSize',28,'Fontweight','normal');
ylabel('Place Field Correlation across epochs','FontSize',28,'Fontweight','normal');


% 2) Corrln vs days

% check and get days with first animal
if isempty(usedays)
    days = unique(psf(1).epochs{1}(:,1));
else
    days = usedays;
end
cumdaycorr = cell(length(days),1);

for anidx = 1:length(useanim)
    an = useanim(anidx);
    if isempty(usedays)
        days = unique(psf(an).epochs{1}(:,1));
    else
        days = usedays;
    end
    
    for nd = 1:length(days)
        d = days(nd);
        cumdaycorr{d} = [ cumdaycorr{d}, allcorr_day{an}{d}];
    end
    
end

% check and get days with first animal
if isempty(usedays)
    days = unique(psf(1).epochs{1}(:,1));
else
    days = usedays;
end
startday = [];
for i=1:length(days)   
    d = days(i);
    cumdaycorrmean(d) = nanmean(cumdaycorr{d});
    cumdaycorrerr(d) = nansem(cumdaycorr{d});
    startday = [startday, d];
end

figure; hold on; redimscreen_figforppt1;
errorbar(startday, cumdaycorrmean, cumdaycorrerr,'ko-','Linewidth',2,'MarkerSize',8);
title('Place Field Stability Across Days','FontSize',28,'Fontweight','normal');
ylabel('Place Field Correlation across epochs','FontSize',24,'Fontweight','normal');
xlabel('Day','FontSize',28,'Fontweight','normal');
axis([startday(1)-0.5 startday(end)+0.5 0 1.0])


% % 3) Corrln vs pairs of days

% check and get days with first animal
if isempty(usedays)
    days = unique(psf(1).epochs{1}(:,1));
else
    days = usedays;
end
cumdaycorr = cell(length(days),1);

if length(days)>=3
    
    for anidx = 1:length(useanim)
        an = useanim(anidx);
        if isempty(usedays)
            days = unique(psf(an).epochs{1}(:,1));
        else
            days = usedays;
        end
        
        for nd = 1:length(days)-1            
            d = days(nd);
            currdays = [d, d+1];
            cumdaycorr{d} = [cumdaycorr{d}, allcorr_day{an}{d}, allcorr_day{an}{d+1}];
        end
    end
    
    % check and get days with first animal
    if isempty(usedays)
        days = unique(psf(1).epochs{1}(:,1));
    else
        days = usedays;
    end
    startday = [];

    for i=1:length(days)
        d = days(i);
        cumdaycorrmean(d) = nanmean(cumdaycorr{d});
        cumdaycorrerr(d) = nansem(cumdaycorr{d});
        startday = [startday, d];
    end
    
    
    figure; hold on; redimscreen_figforppt1;
    errorbar(startday, cumdaycorrmean, cumdaycorrerr,'ko-','Linewidth',2,'MarkerSize',8);
    title('Place Field Stability Across Days','FontSize',28,'Fontweight','normal');
    ylabel('Place Field Correlation across epochs','FontSize',24,'Fontweight','normal');
    xlabel('Day: Doublets','FontSize',28,'Fontweight','normal');
    axis([startday(1)-0.5 startday(end)+0.5 0 1.0])
    
end


% % 4) Corrln vs triplets of days


% check and get days with first animal
if isempty(usedays)
    days = unique(psf(1).epochs{1}(:,1));
else
    days = usedays;
end
cumdaycorr = cell(length(days),1);

if length(days)>=4
    
    for anidx = 1:length(useanim)
        an = useanim(anidx);
        if isempty(usedays)
            days = unique(psf(an).epochs{1}(:,1));
        else
            days = usedays;
        end
        
        for nd = 1:length(days)-2            
            d = days(nd);
            currdays = [d, d+1, d+2];
            cumdaycorr{d} = [cumdaycorr{d}, allcorr_day{an}{d}, allcorr_day{an}{d+1}, allcorr_day{an}{d+2}];
        end
    end
    
    % check and get days with first animal
    if isempty(usedays)
        days = unique(psf(1).epochs{1}(:,1));
    else
        days = usedays;
    end
    startday = [];

    for i=1:length(days)
        d = days(i);
        cumdaycorrmean(d) = nanmean(cumdaycorr{d});
        cumdaycorrerr(d) = nansem(cumdaycorr{d});
        startday = [startday, d];
    end

    figure; hold on; redimscreen_figforppt1;
    errorbar(startday, cumdaycorrmean, cumdaycorrerr,'ko-','Linewidth',2,'MarkerSize',8);
    title('Place Field Stability Across Days','FontSize',24,'Fontweight','normal');
    ylabel('Place Field Corrln across epochs','FontSize',20,'Fontweight','normal');
    xlabel('Day: Triplets','FontSize',16,'Fontweight','normal');
    axis([startday(1)-0.5 startday(end)+0.5 0 1.0])
    
end











% ---------- Old Way To Plot without an field

% disp('Mean LinTrajCorr');
% meancorr = nanmean(allcorr),
% errcorr = nansem(allcorr),
% 
% 
% %% Summarize correlation - 1) Bar plot for all, novel, familiar days, & 2) Corrln across days
% 
% set(0,'defaultaxesfontsize',20);set(0,'defaultaxesfontweight','normal');
% set(0,'defaultaxeslinewidth',2);
% 
% % 1) Bar plot for all, novel, familiar days
% figure; hold on; redimscreen_figforppt1;
% bar(1,meancorr,'k'); errorbar(1,meancorr,errcorr,'k');
% set(gca,'XTick',[1],'XTickLabel',{'All Days'},'FontSize',16,'Fontweight','normal');
% axis([0 2 0 1.0])
% 
% if ~isempty(intersect(days,novel))
%     noveldays = intersect(days,novel);
%     novelcorr = [];
%     for nd=noveldays
%         novelcorr = [novelcorr, allcorr_day{nd}];
%     end
%     
%     bar(2,nanmean(novelcorr),'r'); errorbar(2,nanmean(novelcorr),nansem(novelcorr),'r');
%     set(gca,'XTick',[1 2],'XTickLabel',{'All Days','Novel Days'},'FontSize',16,'Fontweight','normal');
%     axis([0 3 0 1.0])
% end
% 
% if ~isempty(intersect(days,fam))
%     famdays = intersect(days,fam);
%     famcorr = [];
%     for fd=famdays
%         famcorr = [famcorr, allcorr_day{fd}];
%     end
%     
%     if isempty(intersect(days,novel))
%         bar(2,nanmean(famcorr),'g'); errorbar(2,nanmean(famcorr),nansem(famcorr),'g');
%         set(gca,'XTick',[1 2],'XTickLabel',{'All Days','Fam Days'},'FontSize',16,'Fontweight','normal');
%         axis([0 3 0 1.0])
%     else
%         
%         bar(3,nanmean(famcorr),'g'); errorbar(3,nanmean(famcorr),nansem(famcorr),'g');
%         set(gca,'XTick',[1 2 3],'XTickLabel',{'All Days','Novel Days','Fam Days'},'FontSize',16,'Fontweight','normal');
%         axis([0 4 0 1.0])
%     end
% end
% 
% title(['Place Field Stability - Days ' num2str(days') ],'FontSize',20,'Fontweight','normal');
% ylabel('Place Field Corrln across epochs','FontSize',20,'Fontweight','normal');
% 
% % 2) Corrln vs triplets of days
% 
% if length(days)>=4
%     
%     cnt=0;
%     for nd = 1:length(days)-2
%         
%         d = days(nd);
%         cnt=cnt+1;
%         dayscorr = [];
%         currdays = [d, d+1, d+2];
%         dayscorr = [dayscorr, allcorr_day{d}, allcorr_day{d+1}, allcorr_day{d+2}];
%         cumdaycorr{nd} = dayscorr;
%         cumdaycorrmean(nd) = nanmean(dayscorr);
%         cumdaycorrerr(nd) = nansem(dayscorr);
%         startday(nd) = d;
%     end
%     
%     figure; hold on; redimscreen_figforppt1;
%     errorbar(startday, cumdaycorrmean, cumdaycorrerr,'ko-','Linewidth',2,'MarkerSize',8);
%     title('Place Field Stability Across Days','FontSize',24,'Fontweight','normal');
%     ylabel('Place Field Corrln across epochs','FontSize',20,'Fontweight','normal');
%     xlabel('Day: Doublets','FontSize',16,'Fontweight','normal');
%     axis([startday(1)-0.5 startday(end)+0.5 0 1.0])
%     
% end
% 
% % 3) Corrln vs pairs of days
% 
% if length(days)>=4
%     
%     cnt=0;
%     for nd = 1:length(days)-2
%         
%         d = days(nd);
%         cnt=cnt+1;
%         dayscorr = [];
%         currdays = [d, d+1];
%         dayscorr = [dayscorr, allcorr_day{d}, allcorr_day{d+1}];
%         cumdaycorr{cnt} = dayscorr;
%         cumdaycorrmean(cnt) = nanmean(dayscorr);
%         cumdaycorrerr(cnt) = nansem(dayscorr);
%         startday(cnt) = d;
%     end
%     
%     figure; hold on; redimscreen_figforppt1;
%     errorbar(startday, cumdaycorrmean, cumdaycorrerr,'ko-','Linewidth',2,'MarkerSize',8);
%     title('Place Field Stability Across Days','FontSize',28,'Fontweight','normal');
%     ylabel('Place Field Correlation across epochs','FontSize',24,'Fontweight','normal');
%     xlabel('Day: Doublets','FontSize',28,'Fontweight','normal');
%     axis([startday(1)-0.5 startday(end)+0.5 0 1.0])
%     
% end
% 
% % 4) Corrln vs days
% 
% if length(days)>=3
%     
%     cnt=0;
%     for nd = 1:length(days)
%         
%         d = days(nd);
%         cnt=cnt+1;
%         dayscorr = [];
%         currdays = [d];
%         dayscorr = [dayscorr, allcorr_day{d}];
%         cumdaycorr{cnt} = dayscorr;
%         cumdaycorrmean(cnt) = nanmean(dayscorr);
%         cumdaycorrerr(cnt) = nansem(dayscorr);
%         startday(cnt) = d;
%     end
%     
%     figure; hold on; redimscreen_figforppt1;
%     errorbar(startday, cumdaycorrmean, cumdaycorrerr,'ko-','Linewidth',2,'MarkerSize',8);
%     title('Place Field Stability Across Days','FontSize',28,'Fontweight','normal');
%     ylabel('Place Field Correlation across epochs','FontSize',24,'Fontweight','normal');
%     xlabel('Day','FontSize',28,'Fontweight','normal');
%     axis([startday(1)-0.5 startday(end)+0.5 0 1.0])
%     
% end
% 
