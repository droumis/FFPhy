
% 2D Place Field Map by Trajs

% warning('off','all');
close all
runFilterFramework = 0;
saveFilterOutput = runFilterFramework;
loadFilterOutput = 0;
plotOccNormFields = 1;
savefigs=0;
pausefigs = 1;
investInfo = animaldef(lower('Demetris'));
%% ---------------- plotting params --------------------------
colorSet = 'DR1';
fonttype = 'Arial';
titlesize = 16;
axissize = 16;
usecolormap = 'hot';
SpacingHoriz = 0.03;
SpacingVert = 0.02;
% Spacing = 0.00;
% Padding = 0.0;
% MarginLeft = 0.04;
% MarginRight = 0.04;
% MarginTop = 0.09;
% MarginBottom =  0.08;
position = [.1 .1 .5 .7];
% SpacingHorizontal = 0.00;
% SpacingVertical = 0.00;
Spacing = 0.00;
Padding = 0.00;
MarginLeft = 0.05;
MarginRight = 0.05;
MarginTop = 0.15;
MarginBottom =  0.1;
%% ---------------- Data Filters --------------------------
animals = {'JZ1'};
days = [4]; %leave blank to take all days from HP animals and Ndl
epochType = {'run', 'run'};
epochEnvironment = {'wtrack', 'openfield'}; %wtrack

ntAreas = {'ca1', 'mec'};
filtfunction = 'occNormFiring';
eventtype = 'rippleskons';
eventSourceArea = 'ca1';
consensus_numtets = 1;   % minimum # of tets for consensus event detection
minstdthresh = 5;        % STD. how big your ripples are
exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
minripvelocity = 0;
maxripvelocity = 4;
%% ---------------- Paths and Title strings ---------------------------------------------------
filtOutputDirectory = sprintf('%s%s/', investInfo{2}, filtfunction);
figdirectory = sprintf('%s%s/', investInfo{4}, filtfunction);
filenamesave = sprintf('%s_D%s', strjoin(epochEnvironment,'-'), strjoin(arrayfun(@(x) num2str(x),days,'UniformOutput',false),'-')); %add more naming stuff
filename = sprintf('%s_%s.mat', filtfunction, filenamesave);
filenameTitle = strrep(filename,'_', '\_');
%% ---------------- Run FIlter ---------------------------------------------------
if runFilterFramework == 1
    %     epochfilter = 'isequal($type, ''run'') && ~isequal($environment, ''lin'')';
    eptypeEnv = [epochType; epochEnvironment];
    epochfilter = sprintf('((isequal($type, ''%s'')) && (isequal($environment, ''%s''))) || ',eptypeEnv{:});
    yuck = strfind(epochfilter, '||');
    epochfilter = epochfilter(1:yuck(end)-1);  %cut off the trailing '||'
    cellfilter = '($numspikes > 100)';
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
    F = createfilter('animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);
    F = setfilterfunction(F, ['dfa_' filtfunction], {'spikes', 'linpos', 'pos', 'task'});
    F = runfilter(F);
end
%% ---------------- Save Filter Output ---------------------------------------------------
if saveFilterOutput == 1;
    if ~isdir(filtOutputDirectory);
        mkdir(filtOutputDirectory);
    end
    %save the entire workspace for filter provenance
    save(sprintf('%s/%s',filtOutputDirectory, filename), 'F','-v7.3');
    disp(sprintf('filteroutput saved to %s/%s',filtOutputDirectory, filename))
end
%% ---------------- Load Filter Output ---------------------------------------------------
if loadFilterOutput == 1;
    load(sprintf('%s/%s',filtOutputDirectory, filename))
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
        [daytets, daytetInds, daytetInds2 ] = unique(matInds(:,[1 3 4]), 'rows', 'stable');
%         indices = cellfetch(F(iAn).output, 'index');
%         indices = 
%         indices = cell2mat(cellfun(@(x) x(:,[1 3]),indices.values, 'UniformOutput', false)); %scrub the cell id and epoch
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

        
        for icell = 1:length(daytets); %for each xepoch-unique cell
            cellID = daytets(icell, 3);
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
            for iepoch = 1:numeps;
%                 icellFout = F(iAn).output{1}(icellFoutInds(iepoch));
                epoch = epochInds(iepoch,2);
                taskEnv = task{day}{epoch}.environment;
                posData = pos{day}{epoch};
                Xstring = 'x-sm';
                Xcol = find(cell2mat(cellfun(@(x) strcmp(x,Xstring), strsplit(posData.fields, ' '), 'UniformOutput', false)));
                xData = posData.data(:,Xcol);
                Ystring = 'y-sm';
                Ycol = find(cell2mat(cellfun(@(x) strcmp(x,Ystring), strsplit(posData.fields, ' '), 'UniformOutput', false)));
                yData = posData.data(:,Ycol);
                numtrajs = length(F(iAn).output{1}(icellFoutInds(iepoch)).smoothedspikerate);
                xlimx = [min(xData) max(xData)];
                ylimy = [min(yData) max(yData)];
                maxrate = max(cell2mat(cellfun(@(x) max(max(x)), F(iAn).output{1}(icellFoutInds(iepoch)).smoothedspikerate,'UniformOutput', false)));
                for itraj = 1:numtrajs;
                    currtrajFR = []; currtrajOcc = []; currtrajOccFR = [];
                    currtrajFR = F(iAn).output{1}(icellFoutInds(iepoch)).smoothedspikerate{itraj};
                    currtrajOcc = F(iAn).output{1}(icellFoutInds(iepoch)).smoothedoccupancy{itraj};
                    currtrajFR(currtrajFR == -1) = 0;
%                     currtrajOccFR = currtrajFR./currtrajOcc;
%                     currtrajOccFR(isnan(currtrajOccFR)) = 0;
%                     currtrajFR = (currtrajFR)./maxrate; % normalize 0-1 across epochs
                    fronttraj = [];
                    fronttraj(F(iAn).output{1}(icellFoutInds(iepoch)).xticks{itraj}, F(iAn).output{1}(icellFoutInds(iepoch)).yticks{itraj}) = currtrajFR; % currtrajOccFR
                    fronttraj = fronttraj';
                    %plot all pos as light grey lines
                    % double the open field subplot size in x and y
                    sfA = sfA + 1;
                    sfB = sfA;
                    sfC = sfA;
                    sfD = sfA;
                    if strcmp(taskEnv, 'openfield')
                        sfB = sfB + 1;
                        sfC = sfC + 4;
                        sfD = sfD + 5;
                    end
                    subaxis(numeps,4,[sfA sfB sfC sfD], 'Spacing', Spacing, 'SpacingVert', SpacingVert, 'SpacingHoriz', SpacingHoriz, 'Padding', Padding, 'MarginLeft', MarginLeft, 'MarginRight', MarginRight,...
                        'MarginTop', MarginTop, 'MarginBottom', MarginBottom);
                    if strcmp(taskEnv, 'openfield')
                        sfA = sfA + 1;
                    end
                    plot(xData, yData,'Color', [.9 .9 .9], 'LineWidth', 2);
                    hold on;
                    cmap = flipud(colormap(usecolormap));
                    colormap(cmap)
                    mask1 = fronttraj > std2(fronttraj); %keep everything above 1 std opaque
%                     fronttraj = fronttraj;
                    %plot occupancy normalized firing ratemap
                    s = pcolor(fronttraj);
%                     caxis([0 1]) %max(max((I)))-std2(I)]) % saturate top standard dev
                    s.AlphaData = mask1; % make all zero bins transparent
                    s.EdgeColor = 'none';
                    s.FaceAlpha = 'interp';
                    s.FaceColor = 'interp';
%                     lh = plot(xData, yData,'LineWidth', 2); %overlay semi-transparent pos lines
%                     lh.Color=[.9 .9 .9 .05];
                    xlim(xlimx)
                    ylim(ylimy)
                    %                 statematTrajInd = find(linpos{day}{epoch}.statematrix.traj == itraj, 1);
                    
                    %                 wellExitEnter = linpos{day}{epoch}.statematrix.wellExitEnter(statematTrajInd, :);
                    %                 x = mean([linpos{1, 4}{1, 2}.wellSegmentInfo.wellCoord(wellExitEnter(1),1) linpos{1, 4}{1, 2}.wellSegmentInfo.wellCoord(wellExitEnter(2),1)]);
                    %                 y = mean([linpos{1, 4}{1, 2}.wellSegmentInfo.wellCoord(1,2) linpos{1, 4}{1, 2}.wellSegmentInfo.wellCoord(2,2)]);
                    %                 x = mean(wellExitEnter(:,1));
                    % %                 y = mean(wellExitEnter(:,2));
                    % %                 draw a circle for start end well
                    %                 r = 1;
                    %                 ang=0:0.01:2*pi;
                    %                 xp=r*cos(ang);
                    %                 yp=r*sin(ang);
                    %                 plot(x+xp,y+yp, 'Color', [.4 .6 .8]);
                    %                 rectangle('Position',[wellStartEnd(1,1) wellStartEnd(1,2) 2 2],'Curvature',[1 1], 'FaceColor',[0 .5 .5],'EdgeColor','none',...
                    %                     'alpha', .5)
                    %                 plot(x+xp,y+yp, 'Color', [.4 .6 .8]);
                    %                 text(x,y,'\rightarrow','FontSize',15,'FontWeight','bold', 'Color', 'k', 'horizontalAlignment', 'center')
                    if strcmp(taskEnv, 'wtrack')
                        wellStartEnd = F(iAn).output{1}(icellFoutInds(iepoch)).wellCoords{itraj};
                        text(wellStartEnd(1,1),wellStartEnd(1,2),'\wedge','FontSize',10,'FontWeight','bold', 'Color', [.4 .5 .8], 'horizontalAlignment', 'center')
                        text(wellStartEnd(2,1),wellStartEnd(2,2),'\vee','Fontsize',10,'FontWeight','bold', 'Color', [.8 .5 .8], 'horizontalAlignment', 'center')
                    end
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                    set(gca,'ytick',[])
                    set(gca,'yticklabel',[])
                    subTit = title(sprintf('%d spikes %d seconds ',sum(sum(F(iAn).output{1}(icellFoutInds(iepoch)).spikes{itraj})), round(sum(sum(F(iAn).output{1}(icellFoutInds(iepoch)).occupancy{itraj})))));
                    set(subTit, 'Color', 'k', 'FontName', 'Arial', 'FontSize',6, 'FontWeight', 'normal');%, 'horizontalAlignment', 'left', 'verticalAlignment', 'bottom');
                    if itraj == 1
                       ylabel(sprintf('epoch %d', epoch), 'Color', [.8 .8 .8], 'FontSize',12,'FontWeight','bold','FontName', 'Arial');
                    end
                end
%                  subTit = title(sprintf('epoch %d', epoch), 'Color', [.8 .8 .8]);
            end
            %% ---- super title and colorbar----
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            xlab = 'pos x';
            ylab = 'pos y';
            supxlabel = text(.5, .05, xlab, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'Parent', sprtitleax,...
                'Units', 'normalized', 'horizontalAlignment', 'center');
            supylabel = text(.02, .5, ylab, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'rotation', 90, ...
                'Parent', sprtitleax, 'Units', 'normalized', 'horizontalAlignment', 'center');
            caxis([0 maxrate]);
            clrbar = colorbar('location','eastoutside', 'FontSize',6,'FontName', 'Arial');%, 'FontWeight','bold');
            posx1=get(gca,'position');
            posx=get(clrbar,'Position');
            posx(1)= 1-MarginRight+.01;
            posx(2)= MarginBottom;
            posx(3)= .01;
            posx(4)= 1 - MarginTop - MarginBottom;
            set(clrbar,'Position',posx)
            set(gca,'position',posx1)
            sprtit = sprintf('%s %s %s D%d nT%d C%d', filtfunction, strjoin(epochEnvironment,'-'), animalID, day, ntrode, icell);
            iclr = icolors(iIndInds(1),:);
            sprTags = sprintf('%s %s', iarea, isubarea);
            %             iStitle = text(.5, .95, [{sprtit}; ['\color[rgb]{.8 .8 .8} \fontsize{8}', sprintf(' {%s}', filenameTitle)]], ...
            %                 'Parent', sprtitleax, 'Units', 'normalized');
            iStitle = text(.5, .93, [{sprtit};...
                {['\fontsize{12} \color[rgb]' sprintf('{%d %d %d} %s ', iclr, sprTags)]};...
                {['\color[rgb]{.5 .5 .5} \fontsize{8}', sprintf(' {%s}', filenameTitle)]}], 'Parent', sprtitleax, 'Units', 'normalized');
            %             end
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'horizontalAlignment', 'center');
            clrbartit = text(posx(1)+posx(3)/2, posx(2)-MarginBottom/2, 'Hz', 'FontSize',10,'FontWeight','bold','Color','k', 'FontName', 'Arial','horizontalAlignment', 'center');
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