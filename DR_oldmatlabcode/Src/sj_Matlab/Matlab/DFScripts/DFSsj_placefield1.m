
% Place Field Stability and
% Plotting of Open Field Rate and Linearized Place Fields if desired
% Dont Plot all at once for all days - too many plots!

clear; %close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 1; % Figure Options - Plot individual cell place plots
savelinfields = 0; % To save trajdata in [prefix-linfields-day] file for each day
                   % Also mapdata in [prefix-mapfields-day] file for day

savedir = '/data25/sjadhav/RippleInterruption/ProcessedData/';
savefile = [savedir 'REGrp_PlaceFieldsn'];  % n is after std is changed for mapfields
%savefile = [savedir 'RCGrp_PlaceFieldsn'];

% Parameters
minabsvel = 3;  % cm/sec - Most conservative for runs and place fields
minlinvel = 5;

% Plot options
plotanimidx = 3; % To pick animals for plotting
plotdays = 7; % If you only load data when runscript=0 and savedata=0, then this field will supplant days



% If runscript, run Datafilter and save data
if runscript == 1
    
    
    %Animal selection
    %-----------------------------------------------------
    %animals = {'REc,''REd','REe','REf'};
    %animals = {'RCa','RCb','RCc','RCd'}
    %animals = {'REd','REe','REf'};
    %animals = {'RCa','RCb','RCc','RE1'};
    %animals = {'REe'};
    
    %Filter creation
    %-----------------------------------------------------
    
    % Epoch filter
    % -------------
    dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    epochfilter = 'isequal($type, ''run'')';
    
    % Cell filter
    % -----------
    %placecellfilter = '( strcmp($tag, ''PyrSR'') || strcmp($tag, ''PyrSRe'') || strcmp($tag, ''PyrR'') )';
    placecellfilter = 'strcmp($tag, ''PyrSR'')';
    
    % Time filter
    % -----------
    
    % Either use tetlist for riple detection for ripfilter,
    % or use ripple tet filter based on tag in tetinfo marking elec as ripple tet
    
    %tetlist = [1,5,9,10,11,12]; % REf
    %tetlist = [3,4,6,11,12,13]; % REe
    %riptetfilter = '(isequal($area, ''riptet''))';
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    % Abs linear velocity(thrs 5)/ velocity(thrs 3) time filter
    %timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getvelpos', '(($absvel >= 3))'}, {'DFTFsj_getriptimes','($nripples == 0)',tetlist,'minthresh',2} }
    %timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getvelpos', '(($absvel >= 3))'}, {'DFTFsj_getriptimes','($nripples == 0)',[],'tetfilter',riptetfilter,'minthresh',2} }
    
    % Linear velocity time filter implemented in getlinstate
    %timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getriptimes','($nripples == 0)',tetlist,'minthresh',2} };
    timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',2} };
    
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
    pmf = setfilterfunction(psf, 'DFAsj_openfieldrate_tf', {'spikes', 'linpos','pos'}, 'binsize', 1, 'std', 2);
    %pmf = setfilterfunction(psf, 'DFAsj_openfieldrate', {'spikes', 'linpos', 'pos'}, 'binsize', 1, 'minabsvel', 3, 'std', 2);
    
    % Run analysis
    % ------------
    psf = runfilter(psf);  % Place Field Stability
    pmf = runfilter(pmf);  % Place Field Map
    
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

% ---------------------------------

% ------------------------------
% Figure and Font Sizes

forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/PlaceFields/';
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

clr = {'b','r','g','c','m','y','k','r'};

% Saving processed data and figures
savefig1 = 0;
saveprocdatafile = [savedir,'PlaceFieldStabData'];
saveprocdata = 0; saveproctag = 'exp';

% ---------------------------------------

% Get trajdata and days and epochs
trajdata = []; index=[];

for anidx = 1:length(useanim)
    an = useanim(anidx);
    for i=1:length(psf(an).output{1}),
        index{an}(i,:)=psf(an).output{1}(i).index;
        alltrajdata{an}{i}=psf(an).output{1}(i).trajdata;
        allmapdata{an}{i}=pmf(an).output{1}(i);
    end
end

% Initialize
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
    
    % Initialize
    cnt_allneu=0;  % Count neurons across days for current animal
    %%%allcorr=[];
    %%%allcorr_day=[];
    cntneu=0; % These are also across days for current animal
    cntneu_ep2=0; % These are also across days for current animal
    
    for d = 1:length(days)
        
        % Reset linfields in case you are going to save it
        linfields = [];
        mapfields = [];
        
        day = days(d);
        cnt_dayneu = 0;    % This is reset for each day
        
        dayidxs = find(index{an}(:,1)==day);
        % Get tetlist
        tet = unique(index{an}(dayidxs,3));
        
        for elec = 2:length(tet)
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
                    trajdata = alltrajdata{an}{curridx};
                    mapdata = allmapdata{an}{curridx}.smoothedspikerate;
                    maxrate=[];
                    
                    % Save in linfields in case you are going to save it
                    linfields{day}{epoch}{currtet}{currcell}=trajdata;
                    mapfields{day}{epoch}{currtet}{currcell}=allmapdata{an}{curridx}; % Save everything .smoothedspikerate is what you want
                    
                    % Get traj with maximum peak fir rate
                    for i=1:length(trajdata),
                        maxrate(i) = max(trajdata{i}(5:end-5,5));
                    end
                    
                    if figopt1==1,
                        figure(str2num([num2str(day) num2str(currtet) num2str(currcell)]));
                        redimscreen;
                        hold on;
                        subplot(2,length(allepochs),ep); hold on;
                        for i=1:length(trajdata),
                            xax = 2*([1:length(trajdata{i}(:,5))]); % trajdata is in 2cm bins
                            if rem(i,2)==0 % inbound trajectories
                                xax=xax(end:-1:1);
                            end
                            plot(xax,trajdata{i}(:,5),[clr{i} '.-'],'Linewidth',2);
                            %maxrate(i) = max(trajdata{i}(5:end-5,5));
                        end
                        if ep==1
                            title(['Day ' num2str(day) '  Epoch ' num2str(epoch) '  Tet ' num2str(currtet) '  Cell ' num2str(currcell)],'FontSize',24,'Fontweight','bold');
                            %legend('InRealRight','OutRealLeft','OutRealRight','InRealLeft');
                            legend('OutRealLeft','InRealLeft','OutRealRight','InRealRight');
                            xlabel ('Position along linear trajectory (cm)','FontSize',20,'Fontweight','bold');
                            ylabel ('Firing Rate (Hz)','FontSize',24,'Fontweight','bold');
                            
                        else
                            title(['Epoch ' num2str(epoch)],'FontSize',24,'Fontweight','bold');
                            xlabel ('Position along linear trajectory (cm)','FontSize',20,'Fontweight','bold');
                        end
                        subplot(2,length(allepochs),length(allepochs)+ep); hold on;
                        imagesc(flipud(mapdata)); colorbar
                        set(gca,'YLim',[0 100]);
                        set(gca,'XLim',[0 100]);
                        title(['Epoch ' num2str(epoch)], 'FontSize',24,'Fontweight','bold');
                        if ep==1
                            xlabel ('X-position (cm)','FontSize',24,'Fontweight','bold');
                            ylabel ('Y-position (cm)','FontSize',24,'Fontweight','bold');
                        else
                            xlabel ('X-position (cm)','FontSize',24,'Fontweight','bold');
                        end
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
                    
                    
                end % end epoch
            end % end neuron
            
        end % end tet
        
        % FOR DAY
        % Get correlation data for current day neurons
        % Get correlations between place fields
        for nc=1:cnt_dayneu
            n = neuidx{an}{day}(nc);
            
            % Find which trajectory had the max firing rate in all epochs
            peakrate=[]; peaktraj=[];
            [peakrate,peaktraj1] = max(allpeakrates_ep1{an}{n});
            [peakrate,peaktraj2] = max(allpeakrates_ep2{an}{n});
            
            if peaktraj1(1)==peaktraj2(1),
                allmaxtraj{an}{n}=peaktraj1(1);
            else
                allmaxtraj{an}{n}=[peaktraj1(1) peaktraj2(1)];
            end
            
            % 2d correlation between maps
           sizex = min([length(allplacemaps_ep1{an}{n}.yticks),length(allplacemaps_ep2{an}{n}.yticks)]);
           sizey = min([length(allplacemaps_ep1{an}{n}.xticks),length(allplacemaps_ep2{an}{n}.xticks)]);
           corr2d = corr2(allplacemaps_ep1{an}{n}.smoothedspikerate(1:sizex,1:sizey), allplacemaps_ep2{an}{n}.smoothedspikerate(1:sizex,1:sizey));      
            
            % Method 0: Get correlation across all trajectories
            
            place_ep1all=[]; place_ep2all=[]; % Reset fpr each neuron
            for tr=1:4
                place_ep1 = alllintrajs_ep1{an}{n}{tr}(:,5);
                place_ep2 = alllintrajs_ep2{an}{n}{tr}(:,5);
                le = min(length(place_ep1),length(place_ep2));
                place_ep1 = place_ep1(1:le);
                place_ep2 = place_ep2(1:le);
                
                invalid1 = find(isnan(place_ep1)==1);
                invalid2 = find(isnan(place_ep2)==1);
                invalid=unique([invalid1; invalid2]);
                place_ep1(invalid)=[];
                place_ep2(invalid)=[];
                
                place_ep1all = [place_ep1all; place_ep1];
                place_ep2all = [place_ep2all; place_ep2];
            end              
            currcorr = corrcoef(place_ep1all,place_ep2all);
            currcorr = currcorr(1,2);
            
            
            
            % Method1: Use 1st epoch maximal trajectory only
            % ----------------------------------------------          
%             
%             maxtraj = allmaxtraj{an}{n}(1);
%             place_ep1 = alllintrajs_ep1{an}{n}{maxtraj}(:,5);
%             place_ep2 = alllintrajs_ep2{an}{n}{maxtraj}(:,5);
%             le = min(length(place_ep1),length(place_ep2));
%             place_ep1 = place_ep1(1:le);
%             place_ep2 = place_ep2(1:le);
%             
%             invalid1 = find(isnan(place_ep1)==1);
%             invalid2 = find(isnan(place_ep2)==1);
%             invalid=unique([invalid1; invalid2]);
%             place_ep1(invalid)=[];
%             place_ep2(invalid)=[];
%             
%             currcorr = corrcoef(place_ep1,place_ep2);
%             currcorr = currcorr(1,2);
            
            
            % Meth 2: Original: Avg of maximal trajectories if double traj  
            % -------------------------------------------------------------
%             single_traj=0; double_traj=0;            
%             %         % 2d corr
%             
%             % Lin Traj Corr
%             maxtraj = allmaxtraj{an}{n};
%             if length(maxtraj)==1
%                 single_traj = single_traj+1;
%                 place_ep1 = alllintrajs_ep1{an}{n}{maxtraj}(:,5);
%                 place_ep2 = alllintrajs_ep2{an}{n}{maxtraj}(:,5);
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
%                 currcorr = corrcoef(place_ep1,place_ep2);
%                 currcorr = currcorr(1,2);
%             else
%                 double_traj = double_traj+1;
%                 place1_ep1 = alllintrajs_ep1{an}{n}{maxtraj(1)}(:,5);
%                 place1_ep2 = alllintrajs_ep2{an}{n}{maxtraj(1)}(:,5);
%                 le = min(length(place1_ep1),length(place1_ep2));
%                 place1_ep1 = place1_ep1(1:le);
%                 place1_ep2 = place1_ep2(1:le);
%                 
%                 invalid1 = find(isnan(place1_ep1)==1);
%                 invalid2 = find(isnan(place1_ep2)==1);
%                 invalid=unique([invalid1; invalid2]);
%                 place1_ep1(invalid)=[];
%                 place1_ep2(invalid)=[];
%                 
%                 currcorr1 = corrcoef(place1_ep1,place1_ep2);
%                 currcorr1 = currcorr1(1,2);
%                 
%                 
%                 place2_ep1 = alllintrajs_ep1{an}{n}{maxtraj(2)}(:,5);
%                 place2_ep2 = alllintrajs_ep2{an}{n}{maxtraj(2)}(:,5);
%                 le = min(length(place2_ep1),length(place2_ep2));
%                 place2_ep1 = place2_ep1(1:le);
%                 place2_ep2 = place2_ep2(1:le);
%                 
%                 
%                 invalid1 = find(isnan(place2_ep1)==1);
%                 invalid2 = find(isnan(place2_ep2)==1);
%                 invalid=unique([invalid1; invalid2]);
%                 place2_ep1(invalid)=[];
%                 place2_ep2(invalid)=[];
%                 
%                 currcorr2 = corrcoef(place2_ep1,place2_ep2);
%                 currcorr2 = currcorr2(1,2);
%                 currcorr = mean([currcorr1,currcorr2]);
%                 %         place_ep1=sum(place1_ep1,place2_ep1);
%                 %         place_ep2=sum(place1_ep2,place2_ep2);
%                 %         currcorr = corrcoef(place_ep1,place_ep2);
%                 %         currcorr = currcorr(1,2)
%                 
%             end
            
            allcorr{an}(n)=currcorr;
            all2dcorr{an}(n)=corr2d;
            
            allcorr_day{an}{day}(nc)=currcorr;
            all2dcorr_day{an}{day}(nc)=corr2d;
            
        end % end cnt_dayneu
        
        % Save linfields and mapfields for day if asked for
        %savefile = [animdirect,prefix,'linfields',num2str(day)];
        if savelinfields==1
            savefile = sprintf('%s/%slinfields%02d.mat', animdirect, prefix, day);
            save(savefile,'linfields');
            savefile = sprintf('%s/%smapfields%02d.mat', animdirect, prefix, day);
            save(savefile,'mapfields');
        end
    end % end day
    
    
    if savelinfields==1
        disp(['Saved linfields for animal ',prefix]);
    end
    
    if prefix=='RCa'
        for setd=3:8
            allcorr_day{an}{setd}=[]; % It is empty for this animal
        end
    end
        
    
    
    
end % end animal



%--------------------------------------

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

% 1) Bar plot for all, novel, familiar days
figure; hold on;
if forppr==1
    redimscreen_figforppr1;
else
    redimscreen_figforppt1;
end
bar(2,meancorr,'b'); errorbar(2,meancorr,errcorr,'k','LineWidth',3);
set(gca,'XTick',[1.5],'XTickLabel',{'All Days'},'FontSize',xfont,'Fontweight','normal');
axis([0 3 0 1.0])

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
            nd;
            if ((length(allcorr_day{an})>=nd) && (~isempty(allcorr_day{an})))
                novcorr{an} = [novcorr{an}, allcorr_day{an}{nd}];
            end
        end  
        
        novcorrvec = [novcorrvec; novcorr{an}(:)]; 
    end   
end
if ~isempty(novcorrvec)
    bar(2,nanmean(novcorrvec),'r'); errorbar(2,nanmean(novcorrvec),nansem(novcorrvec),'k','LineWidth',3);
    set(gca,'XTick',[1 2],'XTickLabel',{'All Days','Novel Days'},'FontSize',xfont,'Fontweight','normal');
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
            if ((length(allcorr_day{an})>=fd) && (~isempty(allcorr_day{an})))
                famcorr{an} = [famcorr{an}, allcorr_day{an}{fd}];
            end
        end
        
        famcorrvec = [famcorrvec; famcorr{an}(:)];
    end
end
    
    if isempty(novcorrvec)
        bar(2,nanmean(famcorrvec),'r'); errorbar(2,nanmean(famcorrvec),nansem(famcorrvec),'k','LineWidth',3);
        set(gca,'XTick',[1 2],'XTickLabel',{'All Days','Fam Days'},'FontSize',xfont,'Fontweight','normal');
        axis([0 3 0 1.0])
    else
        
        bar(3,nanmean(famcorrvec),'r'); errorbar(3,nanmean(famcorrvec),nansem(famcorrvec),'k','LineWidth',3);
        set(gca,'XTick',[1 2 3],'XTickLabel',{'All Days','Novel Days','Fam Days'},'FontSize',xfont,'Fontweight','normal');
        axis([0 4 0 1.0])
    end

title(['Place Field Stability' ],'FontSize',tfont,'Fontweight','normal');
ylabel('Place Field Correlation across epochs','FontSize',yfont,'Fontweight','normal');

if savefig1==1,
    figfile = [figdir,'PlaceFieldStab_ExpAndCon'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end


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

figure; hold on; 
if forppr==1
    redimscreen_figforppr1;
else
    redimscreen_figforppt1;
end
errorbar(startday, cumdaycorrmean, cumdaycorrerr,'ro-','Linewidth',2,'MarkerSize',12);
title('Place Field Stability Across Days','FontSize',tfont,'Fontweight','normal');
ylabel('Place Field Correlation across epochs','FontSize',yfont,'Fontweight','normal');
xlabel('Day','FontSize',xfont,'Fontweight','normal');
axis([startday(1)-0.5 startday(end)+0.5 0 1.0])

if savefig1==1,
    figfile = [figdir,'PlaceFieldStabvsDays_ConandExp'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end


if saveprocdata == 1;
    eval([saveproctag,'_cumdaycorr = cumdaycorr;']);
    eval([saveproctag,'_allcorr = allcorrvec;']);
    eval([saveproctag,'_allcorr_nov = novcorrvec;']);
    eval([saveproctag,'_allcorr_fam = famcorrvec;']); 
    
    % Save file
    %eval(['save(saveprocdatafile,',saveproctag,'_cumdaycorr,',saveproctag,'_allcorr,'...
    %    ,saveproctag,'_allcorr_nov,',saveproctag,'_allcorr_fam);']);
end

keyboard;



% -------------------------------------------------------------------------
% --------------------------- SAVING --------------------------------------
% -------------------------------------------------------------------------


% %% Save con* and exp* data  and then can do following comparisons such as stats, etc.
% % Save Exp first, then run for Con grp and stop code here. Load the already
% % saved file and then save con* in the same file
% %save(saveprocdatafile,'con*'); %save(saveprocdatafile,'exp*')
% 
% 
% %% Statstics on Con and Exp data in PlaceFieldStabData
% % All corr
% [hall,pall] = ttest2(con_allcorr,exp_allcorr); %pall=0.5606
% 

% Save for Anovan
Allcon=[]; Allexp=[]; Condays=[]; Expdays=[];
% % Daywise corr
for i=1:8 
    concorr = con_cumdaycorr{i};
    expcorr = exp_cumdaycorr{i};
    concorrm(i) = mean(con_cumdaycorr{i});
    expcorrm(i) = mean(exp_cumdaycorr{i});
    concorre(i) = sem(con_cumdaycorr{i});
    expcorre(i) = sem(exp_cumdaycorr{i});
    
    [hday(i),pday(i)] = ttest2(concorr,expcorr); 
    [pdayr(i),hdayr(i)] = ranksum(concorr,expcorr); 
    [hdayk(i),pdayk(i)] = kstest2(concorr,expcorr); 
    
    eval(['concorr',num2str(i),'=concorr;']);
    eval(['expcorr',num2str(i),'=expcorr;']);
    
    % For anovan
    Allcon = [Allcon, concorr];
    Allexp = [Allexp, expcorr];
    Condays = [Condays, i*ones(size(concorr))];
    Expdays = [Expdays, i*ones(size(expcorr))];
        
end

startday=1:8;
figure; hold on; 
redimscreen_figforppr1;
%redimscreen_figforppt1;
errorbar(startday, expcorrm, expcorre,'ro-','Linewidth',2,'MarkerSize',12);
errorbar(startday, concorrm, concorre,'bd-','Linewidth',2,'MarkerSize',12);
title('Place Field Stability Across Days','FontSize',16,'Fontweight','normal');
ylabel('Place Field Correlation across epochs','FontSize',14,'Fontweight','normal');
xlabel('Day','FontSize',14,'Fontweight','normal');
axis([startday(1)-0.5 startday(end)+0.5 0 1.0])

if savefig1==1,
    figfile = [figdir,'PlaceFieldStabvsDays_ConandExp'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end


% % t-test2, kstest2, ranksum
% p1=0.5374, >0.41, >0.44
% p2=0.71, 0.47, 0.64
% p3=0.9997, 0.43, 0.50
% p4=0.8113, 0.78, 0.93
% p5=0.1873, 0.47, 0.23
% p6=0.6227, 0.63, 0.96
% p7=0.0708, 0.13, 0.0538
% p8=0.3242, 0.25, 0.42

% % Can also do anovan to compare curves folowed by multcompare
Allcomb = [Allexp'; Allcon'];
Allgrp = [11*ones(size(Allexp')); 12*ones(size(Allcon'))];
Alldays = [Expdays';Condays'];

[pa,table,stats,terms] = anovan(Allcomb,{Allgrp Alldays},'model','interaction','varnames',{'Grp';'Day'},'display','off'); 
% Grp. F=1.95, p=0.1647
% Day F=2.4, p=0.0229
% Grp*Day, F=0.71, p=0.665

% Multcompare following anova
[c,m,h,nms] = multcompare(stats,'dim',[1 2],'display','off'); % Uses tukey-kramer = hsd by default
% None of the within day comparisons across Exp and Con is significant at 0.05 level
% Only Expd1 and Expd7 for dim=[1 2]
% For dim=2 only, days 3,5,7 are diff from d1

[c,m,h,nms] = multcompare(stats,'dim',[1 2],'ctype','bonferroni','display','off'); % Uses bonferroni = lsd to compensate for multiple comparisons
% Same results




% --------------------------- Oth Plots --------------------------------------





% % % 3) Corrln vs pairs of days
% 
% % check and get days with first animal
% if isempty(usedays)
%     days = unique(psf(1).epochs{1}(:,1));
% else
%     days = usedays;
% end
% cumdaycorr = cell(length(days),1);
% 
% if length(days)>=3
%     
%     for anidx = 1:length(useanim)
%         an = useanim(anidx);
%         if isempty(usedays)
%             days = unique(psf(an).epochs{1}(:,1));
%         else
%             days = usedays;
%         end
%         
%         for nd = 1:length(days)-1            
%             d = days(nd);
%             currdays = [d, d+1];
%             cumdaycorr{d} = [cumdaycorr{d}, allcorr_day{an}{d}, allcorr_day{an}{d+1}];
%         end
%     end
%     
%     % check and get days with first animal
%     if isempty(usedays)
%         days = unique(psf(1).epochs{1}(:,1));
%     else
%         days = usedays;
%     end
%     startday = [];
% 
%     for i=1:length(days)
%         d = days(i);
%         cumdaycorrmean(d) = nanmean(cumdaycorr{d});
%         cumdaycorrerr(d) = nansem(cumdaycorr{d});
%         startday = [startday, d];
%     end
%     
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
% 
% % % 4) Corrln vs triplets of days
% 
% 
% % check and get days with first animal
% if isempty(usedays)
%     days = unique(psf(1).epochs{1}(:,1));
% else
%     days = usedays;
% end
% cumdaycorr = cell(length(days),1);
% 
% if length(days)>=4
%     
%     for anidx = 1:length(useanim)
%         an = useanim(anidx);
%         if isempty(usedays)
%             days = unique(psf(an).epochs{1}(:,1));
%         else
%             days = usedays;
%         end
%         
%         for nd = 1:length(days)-2            
%             d = days(nd);
%             currdays = [d, d+1, d+2];
%             cumdaycorr{d} = [cumdaycorr{d}, allcorr_day{an}{d}, allcorr_day{an}{d+1}, allcorr_day{an}{d+2}];
%         end
%     end
%     
%     % check and get days with first animal
%     if isempty(usedays)
%         days = unique(psf(1).epochs{1}(:,1));
%     else
%         days = usedays;
%     end
%     startday = [];
% 
%     for i=1:length(days)
%         d = days(i);
%         cumdaycorrmean(d) = nanmean(cumdaycorr{d});
%         cumdaycorrerr(d) = nansem(cumdaycorr{d});
%         startday = [startday, d];
%     end
% 
%     figure; hold on; redimscreen_figforppt1;
%     errorbar(startday, cumdaycorrmean, cumdaycorrerr,'ko-','Linewidth',2,'MarkerSize',8);
%     title('Place Field Stability Across Days','FontSize',24,'Fontweight','normal');
%     ylabel('Place Field Corrln across epochs','FontSize',20,'Fontweight','normal');
%     xlabel('Day: Triplets','FontSize',16,'Fontweight','normal');
%     axis([startday(1)-0.5 startday(end)+0.5 0 1.0])
%     
% end




