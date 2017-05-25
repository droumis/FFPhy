%% Run/Plot event-triggered spike raster and histograms
% optional: calculate significance of peri-event modulation for each cell

%% ---------------- Dashboard --------------------------
close all 
runFilterFramework = 1;
saveFilterOutput = runFilterFramework;
loadFilterOutput = 0;
postprocessFilterOutput = 1;
savePostprocessing = postprocessFilterOutput;
loadPostprocessing = 0;
plotfigs = 1;
runAntiAlias = 0;
savefigs= 0;
pausefigs = 1;
investInfo = animaldef(lower('Demetris'));
%% ---------------- plotting params --------------------------
colorSet = 'DR1';
subtitlecolor = [.6 .6 .6];
MarkerFaceColor = [.7 .7 .7];
MarkerFaceAlpha = .7;
markersize = 10;
ripMarkerEdgeColor = ([130, 214, 130])./255;
ripMarkerEdgeAlpha = 1;
lineColor = ([117, 87, 130])./255;
lineAlpha = .8;
siglineColor = ripMarkerEdgeColor;
siglineAlpha = .5;
lineWidth = 3;
areaFaceColor = lineColor;
areaAlpha = lineAlpha;
SpacingHoriz = 0.02;
SpacingVert = 0.0;
Padding = 0.00;
position = [.1 .1 .8 .9];
Spacing = 0.00;
Padding = 0.00;
MarginLeft = 0.05;
MarginRight = 0.01;
MarginTop = 0.1;
MarginBottom =  0.05;
%% ---------------- Data Filters --------------------------
animals = {'JZ1'};
days = [4];
epochType = {'sleep', 'run', 'run'};
epochEnvironment = {'sleep', 'wtrack', 'openfield'}; %wtrack
ntAreas = {'ca1', 'mec', 'por', 'v2l', 'sub'};
filtfunction = 'riptrigspiking';
TF = 1; %'(isequal($validripple, 1))';       % i think this is just grabbing time/info??
eventtype = 'rippleskons';
eventSourceArea = 'ca1';
welldist = [];
exclude_ripples = 0;
consensus_numtets = 1;   % minimum # of tets for consensus event detection
minstdthresh = 2;        % STD. how big your ripples are
exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
minvelocity = 0;
maxvelocity = 4;

window = [0.5 0.5];         % size of psth window (in sec)
binsize = .001;             % size of bins (in sec)
frbinsize = 0.02;
time = (-window(1)-0.5*binsize):binsize:(window(2)+0.5*binsize);

calcstatsPEM = 0;  % calculate significance of each cell's peri-event-modulation

%%% SPARSIFY %%%
% Option to randomly delete a certain proportion of spikes for each cell
% to match a control firing rate, during significance calculation
sparsify = 0;

%% ---------------- Paths and Title strings ---------------------------------------------------
filtOutputDirectory = sprintf('%s%s/', investInfo{2}, filtfunction);
postprocDirectory = sprintf('%s%s/', investInfo{3}, filtfunction);
figdirectory = sprintf('%s%s/', investInfo{4}, filtfunction);
filenamesave = sprintf('%s_D%s', strjoin(epochEnvironment,'-'), strjoin(arrayfun(@(x) num2str(x),days,'UniformOutput',false),'-')); %add more naming stuff
filename = sprintf('%s_%s.mat', filtfunction, filenamesave);
filenameTitle = strrep(filename,'_', '\_');
%% ---------------- Run FIlter ---------------------------------------------------
if runFilterFramework == 1
    timefilter{1} = {'getconstimes', '($cons == 1)', [eventSourceArea eventtype],1,'consensus_numtets',consensus_numtets,...
        'minstdthresh',minstdthresh,'exclusion_dur',exclusion_dur,'minvelocity',minvelocity,'maxvelocity',maxvelocity};
    eptypeEnv = [epochType; epochEnvironment];
    epochfilter = sprintf('((isequal($type, ''%s'')) && (isequal($environment, ''%s''))) || ',eptypeEnv{:});
    yuck = strfind(epochfilter, '||');
    epochfilter = epochfilter(1:yuck(end)-1);  %cut off the trailing '||'
    cellfilter = '($numspikes > 100)';
    iterator = 'singlecellanal';
    F = createfilter('animal', animals,'days',days,'epochs', epochfilter, 'cells', cellfilter, 'excludetime', timefilter,...
        'iterator', iterator);
    F = setfilterfunction(F, 'dfa_geteventtrigspiking', {'spikes',[eventSourceArea eventtype],'pos','task'},...
        'TF',TF,'window',window,'binsize',binsize,'frbinsize',frbinsize,'minthresh',minstdthresh,...
        'maxvelocity',maxvelocity,'minvelocity', minvelocity,'consensus_numtets',consensus_numtets,'welldist',welldist);
    F = runfilter(F);
    %     f.datafilter = whos; %save all filtering parameters in workspace into struct
    % append runscript parameters to the raw output structs
    datafilter_params = paramsstruct(eventSourceArea,eventtype,timefilter,window,binsize,time,consensus_numtets,minstdthresh,exclusion_dur);
    F.datafilter_params = datafilter_params;
end
%% ---------------- Save Filter Output ---------------------------------------------------
if saveFilterOutput == 1;
    if ~isdir(filtOutputDirectory);
        mkdir(filtOutputDirectory);
    end
    save(sprintf('%s/%s',filtOutputDirectory, filename), 'F', '-v7.3');
    disp(sprintf('filteroutput saved to %s/%s',filtOutputDirectory, filename))
end
%% ---------------- Load Filter Output ---------------------------------------------------
if loadFilterOutput == 1;
    load(sprintf('%s/%s',filtOutputDirectory, filename))
    disp(sprintf('filteroutput loaded: %s/%s',filtOutputDirectory, filename))
end
%% ---------------- Load Postproc Output ---------------------------------------------------
if loadPostprocessing == 1;
    load(sprintf('%s/%s',postprocDirectory, filename))
    disp(sprintf('postproc output loaded: %s/%s',postprocDirectory, filename))
end

%% ---------------- Postproc ---------------------------------------------------
if postprocessFilterOutput
    tic
    for iAn = 1;%:length(animals)
        animalinfo = animaldef(lower(animals{iAn}));
        animalID = animalinfo{1,3}; %use anim prefix for name
        FFanimdir =  sprintf('%s',animalinfo{1,2});
        load([FFanimdir, animalID, 'tetinfo']);
        tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
        %find the unique day/tet outputs
        matInds = cell2mat({F.output{1}.index}');
        %         LINmatInds = cell2mat({lF.output{1}.index}');
        [daytetcells, daytetInds, daytetInds2 ] = unique(matInds(:,[1 3 4]), 'rows', 'stable');
        %         indices = cellfetch(F(iAn).output, 'index');
        %         indices =
        %         indices = cell2mat(cellfun(@(x) x(:,[1 3]),indices.values, 'UniformOutput', false)); %scrub the cell id and epoch
        
        %% ---------- for each cell, get all epochs--------------------------------
        % initialize output
        ppF = struct;
        for ic = 1:size(daytetcells,1); %for each allepoch-unique cell
            cellID = daytetcells(ic, 3);
            icellFoutInds = find(ic == daytetInds2);
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
%             iarea = char(cell2mat(strAreas(iIndInds(1)))); %grab the area tags
%             isubarea = num2str(cell2mat(strSubAreas(iIndInds(1))));
            
            load(sprintf('%s%s%s%02d.mat',FFanimdir, animalID, 'task', day));
            %% ---------- get unique epoch environment types  --------------------------------
            eptypes = [];
            for iepoch = 1:numeps;
                epoch = epochInds(iepoch,2);
                eptypes{iepoch,1} = task{day}{epoch}.environment;
            end
            [envtypes,IndA,IndC] = unique(eptypes, 'stable');
            %% ---------- for each epoch env type, concat the data --------------------------------
            % initialize output .fields
            dtc = daytetcells;
            ppF(ic).dtc = dtc(ic,:);
            ppF(ic).epoch_types = [];            % indices here correspond to epoch #s
            ppF(ic).epoch_envs = [];             % indices here correspond to epoch #s
            % all detected events, regardless of epoch or state
            ppF(ic).psth = [];
            ppF(ic).frhist = [];
            ppF(ic).instantFR = [];
            ppF(ic).psthsum = [];
            ppF(ic).instantFRmean = [];
            ppF(ic).posteventmatrix = [];
            ppF(ic).eventduration = [];
            ppF(ic).eventtags = [];
            ppF(ic).eventtags_descript = '[ epochnum eventtime epochtype]';
            ppF(ic).nospikes = 0;
            ppF(ic).noevents = 0;                  % number of events reported for the epochs in which the unit was clustered ("events experienced")
            ppF(ic).epochs = [];
            ppF(ic).run_epochs = [];
            ppF(ic).run_nospikes = nan(1,numeps);
            ppF(ic).run_noevents = nan(1,numeps);
            ppF(ic).sleep_epochs = [];
            ppF(ic).sleep_nospikes = nan(1,numeps);
            ppF(ic).sleep_noevents = nan(1,numeps);
            ppF(ic).time = [];
            for ienv = 1:size(envtypes,1) %for each environment type
                
                ienvTypeInds = find(ienv == IndC);
                ienvFInds = icellFoutInds(ienvTypeInds);
                %                 numtrajs = length(F(iAn).output{1}(ienvFInds(1)).smoothedspikerate);
                %                 ienvPSTH = [];
                %                 posData.data = []; posData.fields = ''; ienvxticks = []; ienvyticks = [];
                for iepienv = 1:length(ienvFInds)
                    epID = epochInds(ienvTypeInds(iepienv),2);
                    %                     iPSTH = F(iAn).output{1}(ienvFInds(iepienv)).psth;
                    iout = F(iAn).output{1}(ienvFInds(iepienv));
                    %                     ienvPSTH = cat(1,ienvPSTH, iPSTH);
                    ppF(ic).time = iout.time;   % (all time vectors are the same..)
                    ppF(ic).frtime = iout.frtime;
                    ppF(ic).psth = [ppF(ic).psth ; iout.psth];
                    ppF(ic).frhist = [ppF(ic).frhist; iout.frhist];
                    ppF(ic).instantFR = [ppF(ic).instantFR ; iout.instantFR];
                    ppF(ic).posteventmatrix = [ppF(ic).posteventmatrix ; iout.posteventmatrix];
                    ppF(ic).eventduration = [ppF(ic).eventduration ; iout.eventduration];
                    ppF(ic).eventtags = [ppF(ic).eventtags ; iout.eventtags];
                    ppF(ic).nospikes = ppF(ic).nospikes + iout.nospikes;
                    ppF(ic).noevents = ppF(ic).noevents + iout.noevents;
                    ppF(ic).epochs = [ppF(ic).epochs epID];
                    % epoch-by-epoch outputs
                    ppF(ic).epoch_types{epID} = iout.epoch_type;
                    ppF(ic).epoch_envs{epID} = iout.epoch_environment;
                    ppF(ic).epoch_nospikes(epID) = sum(iout.nospikes);
                    ppF(ic).epoch_noevents(epID) = sum(iout.noevents);
                    ppF(ic).epoch_noeventspikes(epID) = sum(sum(iout.psth));
                end
                % if a no data for this cell, ignore
                if isempty(ppF(ic).psth)
                    continue
                    disp('psth is empty')
                end
                % Convert .psth field to sparse matrix for smaller saved variable
%                 ppF(ic).psth = sparse(ppF(ic).psth);
                %Sum psth
                ppF(ic).psthsum = sum(ppF(ic).psth,1);
                ppF(ic).instantFRmean = mean(ppF(ic).instantFR,1);
                
                
            end
            disp(sprintf('post processed cell %d of %d',ic, size(daytetcells,1)))
        end
        toc
    end
end
%% ---------------- Save Postproc Output ---------------------------------------------------
if savePostprocessing == 1;
    if ~isdir(postprocDirectory);
        mkdir(postprocDirectory);
    end
    save(sprintf('%s/%s',postprocDirectory, filename), 'ppF', '-v7.3');
    disp(sprintf('postproc saved to %s/%s',postprocDirectory, filename))
end
%% RASTER PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%
if plotfigs
    % create smoothing kernel
    smoothing_length = 10;   % std of gaussian (in ms) used to smooth rip psth
    smoothing_width = round(smoothing_length*.001/binsize);   % smoothing width in number of bins
    kernel = gaussian(smoothing_width,smoothing_width*8);
    for iAn = 1;%:length(animals);
        animalinfo = animaldef(lower(animals{iAn}));
        animalID = animalinfo{1,3}; %use anim prefix for name
        FFanimdir =  sprintf('%s',animalinfo{1,2});
        load([FFanimdir, animalID, 'tetinfo']);
        tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
        dtcInds = cell2mat({ppF(:).dtc}');
        [~, tagIndMap] = ismember(dtcInds(:,[1 2]),tetinfoAll.index(:,[1 3]), 'rows');
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
        %% ---------- plot per cell --------------------------------
        for icell = 1:length(ppF);
            % if you want to save figs but not pause them to take a look.. just keep them invisible.
            % i think this makes it faster and returns keyboard/mouse control during saving
            if savefigs && ~pausefigs;
                ifig = figure('Visible','off','units','normalized','position',position);
            else
                ifig = figure('units','normalized','position',position);
            end
            set(gcf,'color','white')
            load(sprintf('%s%s%s.mat',FFanimdir, animalID, 'cellinfo'));
            % Plotting options
            %             mark_confirmed_sleep_ripples = 0;  % only use if "sleep" data structure was passed into the DFA
            %             timecourse_colormap = 0;
            %             plot_variance = 0;
            %             region = 'ca1';
            idata = ppF(icell);
            try
            [iepenvs, ~, epenvsInds] = unique(idata.epoch_envs, 'stable');
            catch
                disp(sprintf('cell not defined for all epochs: %d %d %d',idata.dtc))
                close all
                continue
            end
            %% ---------- subplot per epoch type --------------------------------
            for iepenv = 1:size(iepenvs,2); %for each epoch type
                epienvIDnums = find(iepenv == epenvsInds);
                [idataienvInds, idataienvbounds] = ismember(idata.eventtags(:,1),epienvIDnums);
                [~, epbounds, ~] = unique(idataienvbounds(idataienvbounds > 0)); %get indices into the raster of the first event for each epoch
                ienvipsth = full(idata.psth(idataienvInds,:));
                sfA = iepenv;
                %                 sfB = iepenv + size(iepenvs,2);
                sfC = iepenv + size(iepenvs,2)*2;
                subaxis(4,size(iepenvs,2),[sfA sfC], 'Spacing', Spacing, 'SpacingVert', SpacingVert, 'SpacingHoriz', SpacingHoriz, 'Padding', Padding, 'MarginLeft', MarginLeft, 'MarginRight', MarginRight,...
                    'MarginTop', MarginTop, 'MarginBottom', MarginBottom);
                epbounds = epbounds(epbounds ~= 1); %get rid of ep bnd at row 1
                for iepbnd = 1:size(epbounds,1)
                    lh = line([idata.time(1)*1000 idata.time(end)*1000], [epbounds(iepbnd) epbounds(iepbnd)],'LineWidth',lineWidth);
                    lh.Color= [lineColor, lineAlpha];
                end
                hold on
                %% plot raster   %%%%%%%%%%%%%%%%%%%
                timevec = idata.time*1000;  % convert to ms
                % obtain bincenter form
                bincenters = timevec(1:(end-1)) + (timevec(2)-timevec(1)) / 2 ;
                % also truncate last bin of histc psth output
                psth = ienvipsth(:,1:(end-1));
                xpoints = psth.*repmat(bincenters,size(psth,1),1);
                xpoints(xpoints==0)=[];
                xpoints = xpoints(:);
                levelvec = (1:size(psth,1))';
                levelmat = repmat(levelvec,1,length(bincenters));
                ypoints = psth.*levelmat;
                ypoints(ypoints==0)=[];
                ypoints = ypoints(:);
                scatter(xpoints,ypoints,markersize,'MarkerEdgeColor','none', 'MarkerFaceColor',MarkerFaceColor,'MarkerFaceAlpha',MarkerFaceAlpha)
                hold on
                postripplematrix = idata.posteventmatrix(idataienvInds,:);
                ripmatrix = [ zeros(size(postripplematrix))  ones(size(postripplematrix,1),1)  postripplematrix ];
                psth = psth.*ripmatrix;
                ypoints2 = ripmatrix.*levelmat;
                plotpsth = psth.*ripmatrix;
                xpoints = plotpsth.*repmat(bincenters,size(plotpsth,1),1);
                xpoints(xpoints==0)=[];
                xpoints = xpoints(:);
                levelvec = (1:size(plotpsth,1))';
                levelmat = repmat(levelvec,1,length(bincenters));
                ypoints = plotpsth.*levelmat;
                ypoints(ypoints==0)=[];
                ypoints = ypoints(:);
                scatter(xpoints,ypoints,markersize,'MarkerEdgeColor', ripMarkerEdgeColor,'MarkerEdgeAlpha',ripMarkerEdgeAlpha,'MarkerFaceColor','none')    % highlight ripple-associated spikes
                %                 dff_plotraster(idata.time,ienvipsth,'markersize', 70, 'color', [.8 .8 .8], 'burstisi',6,'postripplematrix',idata.posteventmatrix(idataienvInds,:))
                hold on
                % axes
                xlim([bincenters(1) bincenters(end)])
                ylim([0 size(ienvipsth,1)])
                set(gca,'YDir','reverse')
                % x ticks
                %set(gca,'XTick',mean([timevec(1) timevec(2)]):500:mean([timevec(end-1) timevec(end)]))
                set(gca,'XTick',[])
                % label axes
                set(gca,'FontSize',8,'FontWeight','bold')
                %xlabel('time (ms)','FontSize',16,'FontWeight','bold')
                %                 ylabel('ripple #','FontSize',10,'FontWeight','bold')
                if iepenv == 1;
                    sylab = {['\color[rgb]{.5 .5 .5} \fontsize{12}', 'ripple #']};
                    ylabel(sylab)
%                     subylab = text(0, .5,sylab, 'FontSize',12,'FontWeight','bold','FontName', 'Arial', 'rotation', 90,...
%                         'Units', 'normalized', 'horizontalAlignment', 'center');
                else
                    set(gca,'ylabel',[])
                end
                %                 set(gcf, 'renderer', 'zbuffer')
                subtit = title(sprintf('%s', iepenvs{iepenv}), 'Color', subtitlecolor, 'FontSize',12,'FontWeight','bold','FontName', 'Arial');
                ylimy = get(gca, 'ylim');
                slh = line([0 0], ylimy, 'LineStyle', '-', 'LineWidth',1.5);
                slh.Color= [[0.5 0.5 0.5], .3];
                %% Create smoothed PSTH of mean firing rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                smoothedpsth = smoothvect(sum(ienvipsth,1)./(binsize*idata.noevents),kernel);
                %                 instFR = idata.instantFR(idataienvInds,:);
                %                 instFRstd = std(instFR,0,1);
                %                 instFRmean = mean(instFR,1);
                sfD = iepenv + size(iepenvs,2)*3;
                subaxis(4,size(iepenvs,2),sfD, 'Spacing', Spacing, 'SpacingVert', SpacingVert, 'SpacingHoriz', SpacingHoriz, 'Padding', Padding, 'MarginLeft', MarginLeft, 'MarginRight', MarginRight,...
                    'MarginTop', MarginTop, 'MarginBottom', MarginBottom);
                ah = area(idata.time,zscore(smoothedpsth),'facecolor',areaFaceColor);
                slh = line([idata.time(1) idata.time(end)], [2 2], 'LineStyle', '--', 'LineWidth',2);
                slh.Color= [siglineColor, siglineAlpha];
                slh = line([idata.time(1) idata.time(end)], [-2 -2], 'LineStyle', '--', 'LineWidth',2);
                slh.Color= [siglineColor, siglineAlpha];
                slh = line([0 0], [-3 3], 'LineStyle', '-', 'LineWidth',1.5);
                slh.Color= [[0.5 0.5 0.5], .3];
                %                 plot(idata.time,smoothvect(instFRmean,kernel),'color',areaFaceColor);
                %                 hold on
                %                 plot(idata.time,smoothvect(instFRmean,kernel)+instFRstd,'color',areaFaceColor);
                %                 hold on
                %                 errfill = fill([idata.time fliplr(idata.time)],[instFRmean'+instFRstd; instFRmean'-flipud(instFRstd)],[0 0 1],'linestyle','none');
                set(ah, 'FaceAlpha', areaAlpha)
                axis tight
                ylim([-3 3])
                set(gca,'FontSize',8,'FontWeight','bold')
                if iepenv == 1;
                sylab = {['\color[rgb]{.5 .5 .5} \fontsize{12}', 'zscore spikerate (std)']};
                ylabel(sylab)
%                 subylab = text(-.2, .5,sylab, 'FontSize',12,'FontWeight','bold','FontName', 'Arial', 'rotation', 90,...
%                     'Units', 'normalized', 'horizontalAlignment', 'center');
                else
                    set(gca,'YTick',[])
                end
                hold on
                
%                 % (Optional) plot variance - should look about the same as psth outline
%                 %smoothing_width = round(smoothing_length*.001/binsize);   % smoothing width in number of bins
%                 %kernel = gaussian(smoothing_width,smoothing_width*8);
%                 if plot_variance
%                     subplot(8,2,[14 15])
%                     psth_variance_smoothed = smoothvect(var(full(idata.psth),1),kernel);
%                     var_mean = mean(psth_variance_smoothed);
%                     y_shift = ymaximum/2;
%                     y_scale = 0.5 * ymaximum/range(psth_variance_smoothed);
%                     psth_variance_toplot = (psth_variance_smoothed - var_mean) * y_scale + y_shift;
%                     h = plot(idata.time,psth_variance_smoothed,'r','linewidth',3);
%                     %         set(gca,'Color',bgcolor)
%                     %         axis tight
%                 end
            end
            
            %% ---- super title and colorbar----
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            xlab = 'Time (s)';
            supxlabel = text(.5, .02, xlab, 'FontSize',12,'FontWeight','bold','Color',[.5 .5 .5], 'FontName', 'Arial', 'Parent', sprtitleax,...
                'Units', 'normalized', 'horizontalAlignment', 'center');
            sprtit = sprintf('%s %s D%d nT%d C%d', filtfunction, animalID, idata.dtc(1), idata.dtc(2), idata.dtc(3));
            iclr = icolors(icell,:);
            iarea = char(cell2mat(strAreas(icell))); %grab the area tags
            isubarea = num2str(cell2mat(strSubAreas(icell)));
            sprTags = sprintf('%s %s', iarea, isubarea);
            %             set(gca,'position',posx1)
            iStitle = text(.5, .95, [{sprtit};...
                {['\fontsize{12} \color[rgb]' sprintf('{%d %d %d} %s ', iclr, sprTags)]};...
                {['\color[rgb]{.5 .5 .5} \fontsize{8}', sprintf(' {%s}', filenameTitle)]}], 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'horizontalAlignment', 'center');
            %% ---- pause, save figs ----
            if runAntiAlias
                myaa; %will drastically improve jagginess but slow things down quite a bit. only use for print figure making
                close(ifig)
            end
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


%%  Calculate significance of modulation

close all

if calcstatsPEM
    
    % option to save ripplemod struct, below
    save_ripplemod = 0;
    daylist = unique(dtc(:,1))';
    tic
    nshuffles = 1000;
    % set up ripplemod struct for storing p values - currently saves all cells and days in one struct
    ripplemod = struct;
    ripplemod.days = daylist;
    ripplemod.event = eventconsname;
    ripplemod.celltype = cellfilter;
    ripplemod.epochtype = epochtype;
    
    % previously had an option to save this for each day, a la "eliripplemod04"
    % ^^ do we want this??
    %       for day = daylist
    %       tmpdtc =  dtc(dtc(:,1)==day,:);
    
    % Modvalues = [day, tet, cell, direction(pos=+1, neg=-1), pvalue, depth-of-modulation, mean-FR-in-window]
    ripplemod.modvalues = [dtc repmat(zeros(size(dtc,1),1),1,4)];
    
    for cellind = 1:size(dtc,1)      % (!!) if doing this for cherry-picked cells, enter cellind here
        %             if dtc(cellind,1) == day
        analyzing_cell = dtc(cellind,:);
        disp(['analyzing cell: ' num2str(analyzing_cell)])
        
        spikedata = full(ppF(cellind).psth);
        numspikes = sum(full(ppF(cellind).psthsum),2);
        
        % mean FR in the plotted window
        windowFR =  numspikes/((length(ppF(cellind).psthsum)*binsize)*size(ppF(cellind).psth,1));
        
        % (Optional) Randomly delete a certain proportion of spikes for each cell
        % to match a control firing rate, specified above.  "sparsify" must equal 1
        control = 0.5; %control FR to match
        if windowFR > control && sparsify
            spikeinds = find(spikedata);
            flip = rand(size(spikeinds));
            deletespikes = flip<=0.5;
            spikedata(spikeinds(deletespikes))=0;
            windowFR = sum(sum(spikedata,1),2)/((length(ppF(cellind).psthsum)*binsize)*size(ppF(cellind).psth,1))
            sparsified = 1;
            disp('sparsified')
        end
        
        % time range to analyze significance... could also use the average event duration, on either side of 0
        range = [0 0.2]; %mean(A(cellind).eventduration);
        % varRange = index for range(1) to range(2) in time series
        varRange = [lookup(range(1),ppF(cellind).time):lookup(range(2),ppF(cellind).time)];
        
        % set windowisrange to 1 to only shuffle spikes within the analysis range, not the whole -0.5 to 0.5
        windowisrange = 0;
        if windowisrange
            tmpspikedata=spikedata(:,varRange); %for shuffling only within the analysis range
            spikedata = tmpspikedata;
        end
        
        binwindow = size(spikedata,2);
        halfbinwindow = size(spikedata,2)/2;
        shufStep = [-halfbinwindow halfbinwindow]; % maximum time to shuffle spikes
        
        allShufPsth = zeros(nshuffles,binwindow);
        
        % smoothing should be the same as the plotted real data
        smoothing_length = 10;   % std of gaussian (in ms) used to smooth rip psth
        smoothing_width = round(smoothing_length*.001/F.runscript_params.binsize);   % smoothing width in number of bins
        kernel = gaussian(smoothing_width,smoothing_width*8);
        
        % circularly shuffle individual ripple-triggered spiketrains by r # of bins,
        for n = 1:nshuffles
            currshuffle = zeros(size(spikedata));
            for rip = 1:size(spikedata,1) % number of rips depends on number of clustered epochs for the cell
                %create an empty "spiketrain"
                dummyspikes = zeros(1,size(spikedata,2));
                %generate random value from the uniform distribution -0.5 to .5 (peri-ripple window in bins)
                r = round(shufStep(1) + (shufStep(2)-shufStep(1)).*rand(1,1));
                %find indices of spiketimes in current train
                spikeinds = find(spikedata(rip,:));
                shiftedspikeinds = spikeinds+r;
                % if times fall off the beginning or end, bring them around circularly
                circindend = find(shiftedspikeinds>size(spikedata,2));
                shiftedspikeinds(circindend)=shiftedspikeinds(circindend)-size(spikedata,2);
                circindstart = find(shiftedspikeinds<1);
                shiftedspikeinds(circindstart)=size(spikedata,2)-(1-shiftedspikeinds(circindstart));
                %fill in dummy spike train
                dummyspikes(shiftedspikeinds) = 1;
                currshuffle(rip,:) = dummyspikes;
            end
            % smoothed the shuffled psth
            smoothedcurrpsth = smoothvect(sum(currshuffle,1)./(binsize*size(spikedata,1)),kernel);
            %      figure
            %      h = bar(1:length(smoothedcurrpsth),smoothedcurrpsth,'facecolor',[0 0 0]);
            allShufPsth(n,:) = smoothedcurrpsth;
        end
        
        %take mean of shuffles
        meanShuf = mean(allShufPsth,1);
        
        %if shuffle was only done within the analysis range, use the whole thing
        if windowisrange
            meanShufRange = meanShuf;
        else
            % otherwise we analyze only the shuffled psth in "range"
            meanShufRange = meanShuf(varRange);
        end
        % figure
        % h = bar(1:length(meanShuf),meanShuf,'facecolor',[0 0 0]);
        
        %find summed squared distance of each shuffle
        meanvarShuf = zeros(nshuffles,1);
        for n=1:nshuffles
            if windowisrange
                allShufRange = allShufPsth(n,:);
            else
                allShufRange = allShufPsth(n,varRange);
            end
            meanvarShuf(n,:) = sum(((allShufRange-meanShufRange).^2),2);
        end
        
        %find summed squared distance of the real data psth from the mean of the shuffles, within the analysis range
        RealPsth = smoothvect(sum(spikedata,1)./(binsize*size(spikedata,1)),kernel);
        if windowisrange
            RealPsthRange = RealPsth;
        else
            RealPsthRange = RealPsth(varRange);
        end
        meanvarReal = sum((RealPsthRange-meanShufRange).^2);
        
        if sparsified
            figure
            h = bar(1:length(RealPsth),RealPsth,'facecolor',[0.5 0.5 0.5]);
            % figure
            % h = bar(varRange,RealPsthRange,'facecolor','b');
            % figure
            % h1 = bar(varRange,meanShufRange,'facecolor','k');
        end
        
        % p value is 1- (what fraction of shuffles have a variance that exceeds the real data)
        p = 1 - sum(meanvarShuf<meanvarReal)/nshuffles
        
        % direction of modulation
        if  sum(RealPsthRange-meanShufRange,2)<0
            direction = -1;
        else
            direction = 1;
        end
        
        ripplemod.modvalues(cellind,4) = direction;
        ripplemod.modvalues(cellind,5) = p;
        ripplemod.modvalues(cellind,6) = meanvarReal; %depth of modulation
        ripplemod.modvalues(cellind,7) = windowFR;
        
        %save ripplemod variable
        if save_ripplemod
            modname = sprintf('%sripplemoddata/%sripplemod_%s_%s_%s.mat',filtOutputDirectory,anim,epochtype,eventconsname,date);
            save(modname,'ripplemod','-mat')
        end
        
    end
    
    %         end
    ripplemod.modvalues
    %     end
    toc
end
