
% 2D Place Field Map by Trajs

% warning('off','all');
close all
runFilterFramework = 1;
saveFilterOutput = runFilterFramework;
loadFilterOutput = 0;
plotOccNormFields = 1;
savefigs= 1;
pausefigs = 0;
investInfo = animaldef(lower('Demetris'));
%% ---------------- plotting params --------------------------
colorSet = 'DR1';
fonttype = 'Arial';
titlesize = 16;
axissize = 16;
arrowsize = 12;
usecolormap = 'hot';
trajname = {'outbound A', 'inbound A', 'outbound B', 'inbound B'};
SpacingHoriz = 0.02;
SpacingVert = 0.02;
Padding = 0.00;
position = [.1 .1 .8 .5];
Spacing = 0.00;
Padding = 0.00;
MarginLeft = 0.05;
MarginRight = 0.05;
MarginTop = 0.15;
MarginBottom =  0.1;
%% ---------------- Data Filters --------------------------
animals = {'JZ1'};
days = [5];
epochType = {'run', 'run'};
epochEnvironment = {'wtrack', 'openfield'}; %wtrack

ntAreas = {'ca1', 'mec', 'por', 'v2l', 'sub'};
filtfunction = 'occNormFiring';
eventtype = 'rippleskons';
eventSourceArea = 'ca1';
consensus_numtets = 1;   % minimum # of tets for consensus event detection
minstdthresh = 5;        % STD. how big your ripples are
exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
minripvelocity = 0;
maxripvelocity = 4;
run2D = 1;
runLin = 1;
%% ---------------- Paths and Title strings ---------------------------------------------------
filtOutputDirectory = sprintf('%s%s/', investInfo{2}, filtfunction);
figdirectory = sprintf('%s%s/', investInfo{4}, filtfunction);
filenamesave = sprintf('%s_D%s', strjoin(epochEnvironment,'-'), strjoin(arrayfun(@(x) num2str(x),days,'UniformOutput',false),'-')); %add more naming stuff
filename = sprintf('%s_%s.mat', filtfunction, filenamesave);
filenameTitle = strrep(filename,'_', '\_');
%% ---------------- Run FIlter ---------------------------------------------------
if runFilterFramework == 1
    eptypeEnv = [epochType; epochEnvironment];
    epochfilter = sprintf('((isequal($type, ''%s'')) && (isequal($environment, ''%s''))) || ',eptypeEnv{:});
    yuck = strfind(epochfilter, '||');
    epochfilter = epochfilter(1:yuck(end)-1);  %cut off the trailing '||'
    cellfilter = '($numspikes > 10)';
    %     cellfilter = sprintf('((isequal($area,''%s'')) && ($clustermetrics.num_events > 100)) || ', ntAreas{:});
    %     yuck = strfind(cellfilter, '||');
    %     cellfilter = cellfilter(1:yuck(end)-1); %cut off the trailing '||'
    %     tetfilter = '(isequal($descrip, ''riptet''))';
    tetfilter = sprintf('(isequal($area,''%s'')) || ', ntAreas{:});
    yuck = strfind(tetfilter, '||');
    tetfilter = tetfilter(1:yuck(end)-1); %cut off the trailing '||'
    %     timefilter{1} = {'dff_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 3))', 6};
    timefilter{1} = {'get2dstate', '(abs($velocity) >= 3)'};
    %     timefilter{1} = {sprintf('get2dstate','($velocity>=%s)', minvel)};
    timefilter{2} = {'getconstimes', '($cons == 0)', [eventSourceArea eventtype],1,'consensus_numtets',consensus_numtets,...
        'minstdthresh',minstdthresh,'exclusion_dur',exclusion_dur,'minvelocity',minripvelocity,'maxvelocity',maxripvelocity};
    iterator = 'singlecellanal';
    iF.datafilter = whos; %save all filtering parameters in workspace into struct
    if run2D
        F = createfilter('animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);
        F = setfilterfunction(F, ['dfa_' filtfunction], {'spikes', 'linpos', 'pos', 'task'});
        F = runfilter(F);
    end
    if runLin
        epochfilter = '((isequal($type, ''run'')) && (isequal($environment, ''wtrack'')))';
        lF = createfilter('animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);
        lF = setfilterfunction(lF, 'dfa_filtercalclinfields', {'spikes', 'linpos', 'cellinfo'});
        lF = runfilter(lF);
    end
end
%% ---------------- Save Filter Output ---------------------------------------------------
if saveFilterOutput == 1;
    if ~isdir(filtOutputDirectory);
        mkdir(filtOutputDirectory);
    end
    %save the entire workspace 6for filter provenance
    save(sprintf('%s/%s',filtOutputDirectory, filename), 'F', 'lF', '-v7.3');
    disp(sprintf('filteroutput saved to %s/%s',filtOutputDirectory, filename))
end
%% ---------------- Load Filter Output ---------------------------------------------------
if loadFilterOutput == 1;
    load(sprintf('%s/%s',filtOutputDirectory, filename))
    disp(sprintf('filteroutput loaded: %s/%s',filtOutputDirectory, filename))
end
%% plot
if plotOccNormFields
    %--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    % make a subplot with the background grey in the upper 4 panels and the trajectory in the lower 4 panels... in ai just overlay
    for iAn = 1:length(F); %per anim
        animalinfo = animaldef(lower(animals{iAn}));
        animalID = animalinfo{1,3}; %use anim prefix for name
        FFanimdir =  sprintf('%s',animalinfo{1,2});
        %% ---------- Get ntrode tags --------------------------------
        load([FFanimdir, animalID, 'tetinfo']);
        tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
        %find the unique day/tet outputs
        matInds = cell2mat({F.output{1}.index}');
        LINmatInds = cell2mat({lF.output{1}.index}');
        [daytetcells, daytetInds, daytetInds2 ] = unique(matInds(:,[1 3 4]), 'rows', 'stable');
        [~, tagIndMap] = ismember(matInds(:,[1 3]),tetinfoAll.index(:,[1 3]), 'rows');
        ntrodeTags = tetinfoAll.values(tagIndMap);
        try
            numsumSupAreas = cellfun(@(x) sum(uint16(x.suparea)), ntrodeTags, 'UniformOutput', false);
            numsumAreas = cellfun(@(x) sum(uint16(x.area)), ntrodeTags, 'UniformOutput', false);
            numsumSubAreas = cellfun(@(x) sum(uint16(x.subarea)), ntrodeTags, 'UniformOutput', false);
            strSupAreas = cellfun(@(x) x.suparea, ntrodeTags, 'UniformOutput', false);
            strAreas = cellfun(@(x) x.area, ntrodeTags, 'UniformOutput', false);
            strSubAreas = cellfun(@(x) x.subarea, ntrodeTags, 'UniformOutput', false);
        catch
            error('all ntrodes need to have a suparea, subarea, and area tag, even if blank')
        end
        icolors = colorPicker(colorSet, strAreas, strSubAreas);
        %% ---------- for each cell, get all epochs--------------------------------
        for icell = 1:length(daytetcells); %for each allepoch-unique cell
            cellID = daytetcells(icell, 3);
            icellFoutInds = find(icell == daytetInds2);
            epochInds = matInds(icellFoutInds,:);
            numeps = size(epochInds,1);
            %         for icell = 1:length(F(iAn).output{1}); %per cell
            %             iInd = F(iAn).output{1}(icell).index;
            iInd = epochInds(1,:);
            day = iInd(1);
            epochF = iInd(2);
            ntrode = iInd(3);
            idayEpTet = [day epochF ntrode]; %F{iAn}.tetout{day}{ntrode}.indices;
            [~, iIndInds] = ismember(idayEpTet(:,[1 3]),matInds(:,[1 3]),'rows');
            iarea = char(cell2mat(strAreas(iIndInds(1)))); %grab the area tags
            isubarea = num2str(cell2mat(strSubAreas(iIndInds(1))));
            %% ---------- plot per cell --------------------------------
            if savefigs && ~pausefigs; %if you want to save figs but not pause them to take a look.. just keep them invisible. i think this makes it faster and returns keyboard/mouse control during saving
                ifig = figure('Visible','off','units','normalized','position',position);
            else
                ifig = figure('units','normalized','position',position);
            end
            set(gcf,'color','white')
            %load pos, cellinfo, task structs
            load(sprintf('%s%s%s%02d.mat',FFanimdir, animalID, 'pos', day));
            load(sprintf('%s%s%s.mat',FFanimdir, animalID, 'cellinfo'));
            sfA = 0;
            load(sprintf('%s%s%s%02d.mat',FFanimdir, animalID, 'task', day));
            %% ---------- get unique epoch environment types  --------------------------------
            eptypes = [];
            for iepoch = 1:numeps;
                epoch = epochInds(iepoch,2);
                eptypes{iepoch,1} = task{day}{epoch}.environment;
            end
            [envtypes,IA,IC] = unique(eptypes, 'stable');
            %% ---------- for each epoch env type, concat the data --------------------------------
            for ienv = 1:length(envtypes) %numeps;
                ienvFInds = icellFoutInds(ienv == IC);
                numtrajs = length(F(iAn).output{1}(ienvFInds(1)).smoothedspikerate);
                ienvSpikeRates = []; posData.data = []; posData.fields = ''; ienvxticks = []; ienvyticks = [];
                for iepienv = 1:length(ienvFInds)
                    smSR = F(iAn).output{1}(ienvFInds(iepienv)).smoothedspikerate;
                    for it = 1:length(smSR)
                        smSR{it}(smSR{it} < 0) = nan;
%                         smSR{it}(isnan(smSR{it})) = 0;
                    end
                    ienvSpikeRates = cat(1,ienvSpikeRates, smSR);
                    ienvxticks = cat(1,ienvxticks, F(iAn).output{1}(ienvFInds(iepienv)).xticks);
                    ienvyticks = cat(1,ienvyticks, F(iAn).output{1}(ienvFInds(iepienv)).yticks);
                    epoch = F(iAn).output{1}(ienvFInds(iepienv)).index(2);
                    posData.fields = pos{day}{epoch}.fields;
                    posData.data = cat(1,posData.data, pos{day}{epoch}.data); % vert cat pos data across epochs
                end
                %% ---------- trim the sizes of each traj-epoch to the smallest of the epochs --------------------------------
                    % first, set any empty epoch trajs to size inf 
                rowsizes = cell2mat(cellfun(@(x) size(x,1), ienvSpikeRates, 'un', 0));
                colsizes = cell2mat(cellfun(@(x) size(x,2), ienvSpikeRates, 'un', 0));
                rowsizes(rowsizes == 0) = inf;
                colsizes(colsizes == 0) = inf;
                [min_size_row, min_index_row] = nanmin(rowsizes,[],1);
                [min_size_col, min_index_col] = nanmin(colsizes,[],1);
                idayienv = [];
                for itrep = 1:numtrajs
                    itrepSpRa = ienvSpikeRates(:,itrep);
                    %lineup trimmed matrices and average across epochs
                    try
                    itreps(1,1,:) = cellfun(@(x) x(1:min_size_row(itrep),1:min_size_col(itrep)), itrepSpRa,'UniformOutput', false);
                    catch %will error/continue on empty traj/epoch
                        continue
                    end
                    itrepsmat = cell2mat(itreps);
                    itrday = nanmean(itrepsmat,3);
                    idayienv.data{itrep} = itrday;
                    %                     itreps(1,1,:) = cellfun(@(x) x(1:min_size_row(itrep)), ienvxticks(min_index_row,itrep),'UniformOutput', false);
                    idayienv.xticks{itrep} = ienvxticks{min_index_row(itrep),itrep};
                    idayienv.yticks{itrep} = ienvyticks{min_index_col(itrep),itrep};
                    %% ---------- get mean, peak for this day/cell/eptype --------------------------------
                    idayienv.peakFR(itrep) = nanmax(nanmax(itrday));
                    idayienv.meanFR(itrep) = nanmean(nanmean(itrday));
                end
                taskEnv = envtypes{ienv};
                Xstring = 'x-sm';
                Xcol = find(cell2mat(cellfun(@(x) strcmp(x,Xstring), strsplit(posData.fields, ' '), 'UniformOutput', false)));
                xData = posData.data(:,Xcol);
                Ystring = 'y-sm';
                Ycol = find(cell2mat(cellfun(@(x) strcmp(x,Ystring), strsplit(posData.fields, ' '), 'UniformOutput', false)));
                yData = posData.data(:,Ycol);
                xlimx = [min(xData) max(xData)];
                ylimy = [min(yData) max(yData)];
                %% --------- plot each trajectory --------------------------------
                for itraj = 1:numtrajs;
                    currtrajFR = []; currtrajOcc = []; currtrajOccFR = [];
                    try
                        currtrajFR = idayienv.data{itraj};
                    catch
                        continue
                    end
                    currtrajFR(currtrajFR == -1) = 0;
                    currtrajFR(isnan(currtrajFR)) = 0;
                    sfA = sfA + 1;
                    sfB = sfA;
                    if strcmp(taskEnv, 'openfield')
                        % double the open field subplot size in x and y
                        opensf = subaxis(2,6,[sfA sfA+1 sfA+6 sfA+7], 'Spacing', Spacing, 'SpacingVert', SpacingVert, 'SpacingHoriz', SpacingHoriz, 'Padding', Padding, 'MarginLeft', MarginLeft, 'MarginRight', MarginRight,...
                            'MarginTop', MarginTop, 'MarginBottom', MarginBottom);
                        alltrajPeakFRopen = max(idayienv.peakFR);
                        alltrajMeanFRopen = mean(idayienv.meanFR);
                    elseif strcmp(taskEnv, 'wtrack')
                        wtracksf = subaxis(2,6,sfA, 'Spacing', Spacing, 'SpacingVert', SpacingVert, 'SpacingHoriz', SpacingHoriz, 'Padding', Padding, 'MarginLeft', MarginLeft, 'MarginRight', MarginRight,...
                            'MarginTop', MarginTop, 'MarginBottom', MarginBottom);
                        alltrajPeakFRwtrack = max(idayienv.peakFR);
                        alltrajMeanFRwtrack = mean(idayienv.meanFR);
                    end
                    currtrajFRxy = []; %put the traj data into the x y pos
                    try
                        currtrajFRxy(idayienv.xticks{itraj}, idayienv.yticks{itraj}) = currtrajFR;
                    catch
                        disp(sprintf('no data found for icell #%d', ienvFInds(1)))
                        continue
                    end
                    alltrajPeakFR = max(idayienv.peakFR);
                    currtrajFRxy = currtrajFRxy'./alltrajPeakFR; %normalize by peak FR per epoch type
                    % plot all pos as light grey lines 
                    plot(xData, yData,'Color', [.9 .9 .9], 'LineWidth', 2); 
                    hold on;
                    cmap = flipud(colormap(usecolormap));
                    colormap(cmap)
                    mask1 = currtrajFRxy > std2(currtrajFRxy); %keep everything above 1 std opaque
                    %plot occupancy normalized firing ratemap
                    s = pcolor(currtrajFRxy);
%                     caxis([0 alltrajPeakFR]) %max(max((I)))-std2(I)]) % saturate top standard dev
                    s.AlphaData = mask1; % make all zero bins transparent
                    s.EdgeColor = 'none';
                    s.FaceAlpha = 'interp';
                    s.FaceColor = 'interp';
                    %                     lh = plot(xData, yData,'LineWidth', 2); %overlay semi-transparent pos lines
                    %                     lh.Color=[.9 .9 .9 .05];
                    caxis([0 1])
                    xlim(xlimx)
                    ylim(ylimy)
                    %                     subTit = title(sprintf('%d spikes %d seconds ',sum(sum(F(iAn).output{1}(ienvFInds(ienv)).spikes{itraj})), round(sum(sum(F(iAn).output{1}(icellFoutInds(ienv)).occupancy{itraj})))));
                    %                     set(subTit, 'Color', 'k', 'FontName', 'Arial', 'FontSize',6, 'FontWeight', 'normal');%, 'horizontalAlignment', 'left', 'verticalAlignment', 'bottom');
                    if itraj == 1 && strcmp(taskEnv, 'wtrack')
                        sylab = [{['\color[rgb]{.5 .5 .5} \fontsize{12}', '2D-W FR all epoch']};...
                            {['\color[rgb]{.5 .5 .5} \fontsize{8}', 'cm']}]; % ...
                        subylab = text(-.2, .5,sylab, 'FontSize',12,'FontWeight','bold','FontName', 'Arial', 'rotation', 90,...
                            'Units', 'normalized', 'horizontalAlignment', 'center');
                        set(gca,'xtick',[], 'xticklabel',[])
                        set(gca,'YColor', [.5 .5 .5], 'FontSize',6,'FontWeight','bold','FontName', 'Arial')
                    elseif itraj == 1 && strcmp(taskEnv, 'openfield')
                        set(gca,'XColor', [.5 .5 .5], 'FontSize',6,'FontWeight','bold','FontName', 'Arial')
                        set(gca,'ytick',[], 'yticklabel',[])
                        xlabel('cm','FontSize',8,'FontWeight','bold','FontName', 'Arial')
                        title(sprintf('%s FR all epochs', taskEnv), 'Color', [.5 .5 .5], 'FontSize',12,'FontWeight','bold','FontName', 'Arial');                        
                    else
                        set(gca,'xtick',[], 'xticklabel',[],'ytick',[], 'yticklabel',[])
                    end
                    % for w track, plot 2D and 1D subplots for each
                    % trajectory
                    if strcmp(taskEnv, 'wtrack')
                        try %will error on empty epoch/traj.. just continue through
                        wellStartEnd = F(iAn).output{1}(ienvFInds(1)).wellCoords{itraj};
                        text(wellStartEnd(1,1),wellStartEnd(1,2),'\wedge','FontSize',arrowsize,'FontWeight','bold', 'Color', [.4 .8 .5], 'horizontalAlignment', 'center')
                        text(wellStartEnd(2,1),wellStartEnd(2,2),'\vee','Fontsize',arrowsize,'FontWeight','bold', 'Color', [.8 .5 .8], 'horizontalAlignment', 'center')
                        text(wellStartEnd(1,1),wellStartEnd(1,2),'\wedge','FontSize',arrowsize,'FontWeight','bold', 'Color', [.4 .8 .5], 'horizontalAlignment', 'center')
                        text(wellStartEnd(2,1),wellStartEnd(2,2),'\vee','Fontsize',arrowsize,'FontWeight','bold', 'Color', [.8 .5 .8], 'horizontalAlignment', 'center')
                        catch
                            disp('no data for this ep traj cell')
                        end
                        title(sprintf('%s', trajname{itraj}), 'Color', [.8 .8 .8], 'FontSize',10,'FontWeight','bold','FontName', 'Arial');
                        %  gather the linear firing rates from the lF output struct and plot
                        dayntcell = iInd(:,[1 3 4]);
                        LINmatdayntcell = LINmatInds(:,[1 3 4]);
                        LINIndMap = [];
                        LINIndMap = find(ismember(LINmatdayntcell, dayntcell, 'rows'));
                        if ~isempty(LINIndMap)
                            sfB = sfA + 6;
                            subaxis(2,6,sfB, 'Spacing', Spacing, 'SpacingVert', SpacingVert, 'SpacingHoriz', SpacingHoriz, 'Padding', Padding, 'MarginLeft', MarginLeft, 'MarginRight', MarginRight,...
                                'MarginTop', MarginTop, 'MarginBottom', MarginBottom);
                            for iLINInd = 1:length(LINIndMap)
                                trajdatafields = lF(iAn).output{1}(LINIndMap(iLINInd)).trajdatafields;
                                trajdata = lF(iAn).output{1}(LINIndMap(iLINInd)).trajdata;
                                linFRdatastring = 'sm-occnormfiringrate';
                                linlocstring = 'locationbin';
                                linFRcol = find(cell2mat(cellfun(@(x) strcmp(x,linFRdatastring), strsplit(trajdatafields, ' '), 'UniformOutput', false)));
                                linloccol = find(cell2mat(cellfun(@(x) strcmp(x,linlocstring), strsplit(trajdatafields, ' '), 'UniformOutput', false)));
                                ieplinFRdata = trajdata{itraj}(:,linFRcol);
                                ieplinlocdata = trajdata{itraj}(:,linloccol);
                                if itraj == 2 || itraj == 4;
                                    ieplinFRdata = flipud(ieplinFRdata);
                                end
                                plot(ieplinlocdata,ieplinFRdata, 'LineWidth', 2)
                                hold on;
                            end
                            set(gca,'YColor', [.5 .5 .5], 'FontSize',6,'FontWeight','bold','FontName', 'Arial')
                            xlimlin = [ieplinlocdata(1,1) ieplinlocdata(end,1)];
                            xlim(xlimlin)
                            set(gca, 'YLim', [0, get(gca, 'YLim') * [0; 1]])
                            text(xlimlin(1),0,'\wedge','FontSize',arrowsize,'FontWeight','bold', 'Color', [.4 .8 .5], 'horizontalAlignment', 'center', 'verticalAlignment', 'cap')
                            text(xlimlin(1),0,'\wedge','FontSize',arrowsize,'FontWeight','bold', 'Color', [.4 .8 .5], 'horizontalAlignment', 'center', 'verticalAlignment', 'cap')
                            %                                 text(xlimlin(1),0,'\wedge','FontSize',10,'FontWeight','bold', 'Color', [.4 .8 .5], 'horizontalAlignment', 'center', 'verticalAlignment', 'cap')
                            text(xlimlin(2),0,'\vee','Fontsize',arrowsize,'FontWeight','bold', 'Color', [.8 .5 .8], 'horizontalAlignment', 'center', 'verticalAlignment', 'cap')
                            text(xlimlin(2),0,'\vee','Fontsize',arrowsize,'FontWeight','bold', 'Color', [.8 .5 .8], 'horizontalAlignment', 'center', 'verticalAlignment', 'cap')
                            %                                 text(xlimlin(2),0,'\vee','Fontsize',10,'FontWeight','bold', 'Color', [.8 .5 .8], 'horizontalAlignment', 'center', 'verticalAlignment', 'cap')
                            if itraj == 1
                                %                                     ylabel({{sprintf('linFR /ep /trj', taskEnv), 'Color', [.8 .8 .8], 'FontSize',12,'FontWeight','bold','FontName', 'Arial'},...
                                %                                         {sprintf('Hz', taskEnv), 'Color', [.5 .5 .5], 'FontSize',10,'FontWeight','bold','FontName', 'Arial'}});
                                sylab = [{['\color[rgb]{.5 .5 .5} \fontsize{12}', '1D-W FR per epoch']};...
                                    {['\color[rgb]{.5 .5 .5} \fontsize{8}', 'Hz']}];
                                subylab = text(-.2, .5,sylab, 'FontSize',12,'FontWeight','bold','FontName', 'Arial', 'rotation', 90,...
                                    'Units', 'normalized', 'horizontalAlignment', 'center');
                                set(gca,'XColor', [.5 .5 .5], 'FontSize',6,'FontWeight','bold','FontName', 'Arial')
                                xlabel('cm','FontSize',8,'FontWeight','bold','FontName', 'Arial')
                            else
                                set(gca,'xtick',[], 'xticklabel',[])%,'ytick',[], 'yticklabel',[])
                            end
                            
                        end
                    end
                    
                end
            end
            %% ---- super title and colorbar----
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            %             xlab = 'pos x';
            %             ylab = 'pos y';
            %             supxlabel = text(.5, .05, xlab, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'Parent', sprtitleax,...
            %                 'Units', 'normalized', 'horizontalAlignment', 'center');
            %             supylabel = text(.02, .5, ylab, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'rotation', 90, ...
            %                 'Parent', sprtitleax, 'Units', 'normalized', 'horizontalAlignment', 'center');
            %             caxis([0 maxrate]);
%             opensfpos = get(opensf,'position');
            posx1=get(gca,'position');
            clrbar = colorbar('location','eastoutside', 'FontSize',6,'FontName', 'Arial', 'ylim',[0 1], 'ytick', []);%, 'FontWeight','bold');   
            posx=get(clrbar,'Position');
            posx(1)= 1-MarginRight+.01; 
            posx(2)= MarginBottom;
            posx(3)= .01;
            posx(4)= 1 - MarginTop - MarginBottom;
            set(clrbar,'Position',posx)
            %             currytick = get(clrbar,'position');
            %             set(clrbar,'YTicklabel', newytick);
            try
                openFR = axes('position',posx,'yaxislocation','right','color','none','ylim',[0 alltrajPeakFRopen],'FontSize',6,'FontName', 'Arial', 'xtick',[]);
                waxFR = axes('position',posx,'yaxislocation','left','color','none','ylim',[0 alltrajPeakFRwtrack],'FontSize',6,'FontName', 'Arial', 'xtick',[]);
            catch
            end
            sprtit = sprintf('%s %s %s D%d nT%d C%d', filtfunction, strjoin(epochEnvironment,'-'), animalID, day, ntrode, cellID);
            iclr = icolors(iIndInds(1),:);
            sprTags = sprintf('%s %s', iarea, isubarea);
%             set(gca,'position',posx1)
            iStitle = text(.5, .93, [{sprtit};...
                {['\fontsize{12} \color[rgb]' sprintf('{%d %d %d} %s ', iclr, sprTags)]};...
                {['\color[rgb]{.5 .5 .5} \fontsize{8}', sprintf(' {%s}', filenameTitle)]}], 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'horizontalAlignment', 'center');
            clrbartit = text(posx(1)+posx(3)/2, posx(2)-MarginBottom/3, 'Hz', 'FontSize',8,'FontWeight','bold','Color',[.5 .5 .5], 'FontName', 'Arial','horizontalAlignment', 'center', 'Parent', sprtitleax,'Units', 'normalized');
            clrbartit2 = text(posx(1)+posx(3)/2, posx(2)+posx(4), 'W     O', 'FontSize',8,'FontWeight','bold','Color',[.5 .5 .5], 'FontName', 'Arial','horizontalAlignment', 'center', 'Parent', sprtitleax,'Units', 'normalized');
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                if ~isdir(figdirectory);
                    mkdir(figdirectory);
                end
                if ~isdir([figdirectory filenamesave]);
                    mkdir([figdirectory filenamesave]);
                end
                sprtitsave = strrep(sprtit,' ', '_'); %put '_' character in place of whitespace for filename
                set(gcf, 'PaperPositionMode', 'auto'); %saves the png in the dimensions defined for the fig
                currfigfile = sprintf('%s%s/%s',figdirectory, filenamesave, sprtitsave);
                print(currfigfile,'-dpng', '-r0')
                disp(sprintf('plot %s saved', sprtit))
            end
            close all
        end
    end
    
end