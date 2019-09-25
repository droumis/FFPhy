
% Changing how I plot openfield rate from version 1. I call twodoccupanct instead of openfieldrate
% Open field is now separated by trajectory and a separate plot is made.

% Plotting of Open Field Rate and Linearized
% Place Fields if desired

clear;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1=0; % Figure Options

savedir = '/datatmp/kkay/ProcessedData/';
savefile = [savedir 'PlaceFieldsTrajMaps'];

minabsvel = 3;  % cm/sec - Most conservative for runs and place fields
minlinvel = 5;

plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days
plotanimidx = [1]; % To pick animals for plotting


% If runscript, run Datafilter and save data
if runscript == 1
    
    
    %Animal selection
    %-----------------------------------------------------
    animals = {'Chapati'};
    
    %Filter creation
    %--------------------------------------------------------
    
    % Epoch filter
    % -------------
    dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    epochfilter = 'isequal($type, ''run'')';
    
    % Cell filter
    % -----------
    %placecellfilter = '( strcmp($tag, ''PyrSR'') || strcmp($tag, ''PyrSRe'') || strcmp($tag, ''PyrR'') )';
    %placecellfilter = 'strcmp($tag, ''PyrSR'')';
    ca1cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 100))';
    ca2cellfilter = '(isequal($area, ''CA2'') && ($meanrate < 100))';
    ca3cellfilter = '(isequal($area, ''CA3'') && ($meanrate < 100))';

    
    % Time filter
    % -----------
    
    % Either use tetlist for riple detection for ripfilter,
    % or use ripple tet filter based on tag in tetinfo marking elec as ripple tet
    
    tetlist = [5,6,11,13,14]; % Chapati "pure" CA1 tetrodes % [5,6,11,13,14]
    
    %riptetfilter = '(isequal($area, ''riptet''))';
    %riptetfilter = '(isequal($descrip, ''riptet''))';
    
    % Abs linear velocity(thrs 5)/ velocity(thrs 3) time filter
    %timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getvelpos', '(($absvel >= 3))'}, {'DFTFsj_getriptimes','($nripples == 0)',tetlist,'minthresh',2} }
    %timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getvelpos', '(($absvel >= 3))'}, {'DFTFsj_getriptimes','($nripples == 0)',[],'tetfilter',riptetfilter,'minthresh',2} }
    
    % Linear velocity time filter implemented in getlinstate
    timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getriptimes','($nripples == 0)','tetrodes',tetlist,'minthresh',2} };
    %timefilter = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6}, {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',2} };
    
    % Iterator
    % --------
    iterator = 'singlecellanal';
    
    % Filter creation
    % ----------------
    ca1f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca1cellfilter,'excludetime', timefilter,'iterator', iterator);
    ca2f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca2cellfilter,'excludetime', timefilter,'iterator', iterator);
    ca3f = createfilter('animal',animals,'epochs',epochfilter,'cells',ca3cellfilter,'excludetime', timefilter,'iterator', iterator);
    % psf = createfilter('animal',animals,'days',dayfilter,'epochs',epochfilter,'cells',placecellfilter,'iterator', iterator);
    
    % Set analysis function
    % ----------------------

    ca1lf = setfilterfunction(ca1f, 'DFAsj_filtercalclinfields_tf', {'spikes', 'linpos'}, 'binsize', 2);
    %psf = setfilterfunction(psf, 'DFAsj_filtercalclinfields', {'spikes', 'linpos', 'pos'}, 'binsize', 2, 'minabsvel', 3);
    ca1mf = setfilterfunction(ca1f, 'DFAsj_twodoccupancy', {'spikes', 'linpos','pos'}, 'binsize', 1, 'std', 3);
    %pmf = setfilterfunction(psf, 'DFAsj_openfieldrate_tf', {'spikes', 'linpos', 'pos'}, 'binsize', 1, 'minabsvel', 3, 'std', 3);
    
    ca2lf = setfilterfunction(ca2f, 'DFAsj_filtercalclinfields_tf', {'spikes', 'linpos'}, 'binsize', 2);
    ca2mf = setfilterfunction(ca2f, 'DFAsj_twodoccupancy', {'spikes', 'linpos','pos'}, 'binsize', 1, 'std', 3);
    ca3lf = setfilterfunction(ca2f, 'DFAsj_filtercalclinfields_tf', {'spikes', 'linpos'}, 'binsize', 2);
    ca3mf = setfilterfunction(ca2f, 'DFAsj_twodoccupancy', {'spikes', 'linpos','pos'}, 'binsize', 1, 'std', 3);

    
    % Run analysis
    % ------------
    ca1lf = runfilter(ca1lf);  % Place Field Stability
    ca1mf = runfilter(ca1mf);  % Place Field Map
    
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear figopt1 runscript plotdays plotanimidx savedata
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

% if savedata == 1
%     figopt1 = 0; plotdays = [];
% end
% if exist('runscript') & exist('figopt1')
%     if (runscript == 0) && (figopt1==1)
%         days = plotdays;
%     end
% end


% ---------------------------------


% Set up trajdata, index, and allmap data
trajdata = []; index=[]; allmapdata=[];
for anidx = 1:length(useanim)
    an = useanim(anidx);
    for i=1:length(ca1lf(an).output{1}),
        
        index{an}(i,:)=ca1lf(an).output{1}(i).index;
        alltrajdata{an}{i}=ca1lf(an).output{1}(i).trajdata;
        allmapdata{an}{i}=ca1mf(an).output{1}(i);
        
    end
end

% Initialize
clr = {'b','r','g','m','c','y','k','r'};
novel=1:3;
fam=4:8;


for anidx = 1:length(useanim)
    an = useanim(anidx);
    % Get days and epochs for current animal
    if isempty(usedays) 
        days = unique(ca1lf(an).epochs{1}(:,1));
    else
        days = usedays;
    end
    
    dayepoch_pairs=ca1lf(an).epochs{1};
              %%  allepochs = unique(ca1lf(an).epochs{1}(:,2));  old code from shantanu -- bad if RUN epochs differ from day to day // kk 3.26.13 
    %days = unique(index(:,1));
    %allepochs = unique(index(:,2));
    
    %Initialize
    cnt_allneu=0;  % Count neurons across days
    %%%allcorr=[];
    %%%allcorr_day=[];
    cntneu=0; % These are also across days
    cntneu_ep2=0; % These are also across days
    
    for d = 1:length(days)
        
        day = days(d)
        epoch_indices = find(dayepoch_pairs(:,1)==day)
        cnt_dayneu = 0;
        
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
                
                for ep = 1:length(epoch_indices)
                    epoch = dayepoch_pairs(ep,2);
                    curridx = find( index{an}(:,1)==day & index{an}(:,2)==epoch & index{an}(:,3)==currtet & index{an}(:,4)==currcell);
                    disp(sprintf('%d %d %d %d',day,epoch,currtet,currcell))
                    if isempty(curridx)
                        continue
                    end
                    trajdata = alltrajdata{an}{curridx};        %% added {an} after alltrajdata -- not sure why it was otherwise
                    mapdata = allmapdata{an}{curridx}.smoothedspikerate;
                    maxrate=[];
                    % Get traj with maximum peak fir rate
                    for i=1:length(trajdata),
                        maxrate(i) = max(trajdata{i}(5:end-5,5));
                    end
                    
                    if ~exist('figopt1'), figopt1=0; end
                    if figopt1==1,
                        
                        % Linearized traj
                        figure(str2num([num2str(day) num2str(currtet) num2str(currcell) num2str(0)]));
                        redimscreen_2horsubplots;
                        hold on;
                        subplot(1,length(allepochs),ep); hold on;
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
                        % 1) Each traj spearately - 1 fig for each epoch
                        %                     figure(str2num([num2str(day) num2str(currtet) num2str(currcell) num2str(ep)]));
                        %                     titlestr = ['Out RealLeft';'In  RealLeft';'OutRealRight';'In RealRight'];
                        %                     for tr = 1:length(mapdata)
                        %                         subplot(2,length(allepochs),tr); hold on;
                        %                         imagesc(flipud(mapdata{tr})); colorbar
                        %                         title(['Epoch' num2str(epoch) ': ' titlestr(tr,:)], 'FontSize',14,'Fontweight','bold');
                        %                         if ep==1
                        %                             xlabel ('X-position (cm)','FontSize',14,'Fontweight','bold');
                        %                             ylabel ('Y-position (cm)','FontSize',14,'Fontweight','bold');
                        %                         else
                        %                             xlabel ('X-position (cm)','FontSize',14,'Fontweight','bold');
                        %                         end
                        %                     end
                        
                        % 2) Combine Outs and Ins / Or alternatively, Rights and Lefts. Just 1 fig
                        
                        % Rightmap
                        if isempty(mapdata{3}),
                            rightmap = mapdata{4};
                        elseif isempty(mapdata{4}),
                            rightmap = mapdata{3};
                        else
                            a=min(size(mapdata{3},1),size(mapdata{4},1));
                            b=min(size(mapdata{3},2),size(mapdata{4},2));
                            rightmap = (mapdata{3}(1:a,1:b) + mapdata{4}(1:a,1:b))./2;
                        end
                        
                        % Leftmap
                        if isempty(mapdata{1}),
                            leftmap = mapdata{2};
                        elseif isempty(mapdata{2}),
                            leftmap = mapdata{1};
                        else
                            a=min(size(mapdata{1},1),size(mapdata{2},1));
                            b=min(size(mapdata{1},2),size(mapdata{2},2));
                            leftmap = (mapdata{1}(1:a,1:b) + mapdata{2}(1:a,1:b))./2;
                        end
                        
                        figure(str2num([num2str(day) num2str(currtet) num2str(currcell) num2str(1)]));
                        redimscreen
                        subplot(2,length(allepochs),2*ep-1); hold on;
                        imagesc(flipud(rightmap)); colorbar
                        title(['Epoch' num2str(epoch) ': Rightmap'], 'FontSize',14,'Fontweight','bold');
                        xlabel ('X-position (cm)','FontSize',14,'Fontweight','bold');
                        ylabel ('Y-position (cm)','FontSize',14,'Fontweight','bold');
                        subplot(2,length(allepochs),2*ep); hold on;
                        imagesc(flipud(leftmap)); colorbar
                        title(['Epoch' num2str(epoch) ': Leftmap'], 'FontSize',14,'Fontweight','bold');
                        xlabel ('X-position (cm)','FontSize',14,'Fontweight','bold');
                        ylabel ('Y-position (cm)','FontSize',14,'Fontweight','bold');
                        
                        
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
                
                if figopt1==1,
                    keyboard;  % Pause and return control after neurons plot - 2 or 3 figs in all
                end
                
            end % end neuron
            
        end % end tet
        
        % FOR DAY
        % Get correlation data for current day neurons
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
            
            % Get correlations between place fields
            single_traj=0; double_traj=0;
            
            %         % 2d corr
            sizex = min([length(allplacemaps_ep1{an}{n}.yticks),length(allplacemaps_ep2{an}{n}.yticks)]);
            sizey = min([length(allplacemaps_ep1{an}{n}.xticks),length(allplacemaps_ep2{an}{n}.xticks)]);
            corr2d = corr2(allplacemaps_ep1{an}{n}.smoothedspikerate(1:sizex,1:sizey), allplacemaps_ep2{an}{n}.smoothedspikerate(1:sizex,1:sizey));
            
            % Method 0: Get correlation across all trajectories
            place_ep1all=[]; place_ep2all=[]; % Reset for each neuron
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
            
%             % Method1: Use 1st epoch maximal trajectory only
%             % ----------------------------------------------   
%             %             
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
%             
%               
%             % Meth 2: Original: Avg of maximal trajectories if double traj  
%             % -------------------------------------------------------------
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
        
        
    end % end day
    
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