

function [F, ppF] = dfs_riptrigspiking(Fp, varargin)
%% Run event-triggered spike raster and histograms
% optional: calculate significance of peri-event modulation for each cell

%% ---------------- Dashboard --------------------------
% F = []; ppF = [];
runFilterFramework = 0;
saveFilterOutput = runFilterFramework;
loadFilterOutput = 0;

combineEpochs = 0;
saveCombinedEpochs = 0;
loadCombinedEpochs = 0;

plotfigs = 0;
runAntiAlias = 0;
savefigs= 0;
pausefigs = 0;

calcstatsPEM = 0;  % calculate significance of each cell's peri-event-modulation
%%% SPARSIFY %%%
% Option to randomly delete a certain proportion of spikes for each cell
% to match a control firing rate, during significance calculation
sparsify = 0;

% ---------------- Data Filters -----------------------------------------------

if ~isempty(varargin)
    assign(varargin{:})
end
% ---------------- Paths ---------------------------------------------------
paths = make_paths(Fp);

%% ---------------- Get Data---------------------------------------------------
% ---------------- Run FIlter -------------------------------------------------
if runFilterFramework == 1    
    F = createfilter('animal',Fp.animals,'days',Fp.days,'epochs', ... 
        Fp.epochfilter, 'cells',Fp.cellfilter, 'excludetime', Fp.timefilter, ...
        'iterator', Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, {'spikes', ...
        Fp.eventDataLabel ,'pos','task'}, 'TF',Fp.TF,'window',Fp.window, ...
        'binsize',Fp.binsize,'frbinsize',Fp.frbinsize,'minthresh', ...
        Fp.minstdthresh,'maxvelocity',Fp.maxvelocity,'minvelocity', ...
        Fp.minvelocity, 'consensus_numtets',Fp.consensus_numtets,'welldist', ...
        Fp.welldist);
    F = runfilter(F);
    for a = 1:length(F)
        F(a).datafilter_params = Fp;
    end
end

% ---------------- Save Filter Output ----------------------------------------
if saveFilterOutput == 1;
    save_filter_output(F, paths.filtOutputDirectory, paths.filenamesave)
end
% ---------------- Load Filter Output ----------------------------------------
if loadFilterOutput == 1;
    F = load_filter_output(paths.filtOutputDirectory, paths.filenamesave);
end
% ---------------- combine epochs ---------------------------------------------
if combineEpochs
    if exist('F', 'var')
        ppF = combine_epochs(F, Fp, saveCombinedEpochs, paths);
    else
        error('create or load data filter output to combine epochs \n')
    end
end
% ---------------- Load combined epochs ---------------------------------------------------
if loadCombinedEpochs == 1;
    ppF = load_filter_output(paths.resultsDirectory, paths.filenamesave, ...
        'filetail', '_combEps');
end



end
% 
% 
% %% ---------------- Plot---------------------------------------------------
% if plotfigs
%     fprintf('u sure? \n')
%     pause
%     Pp = load_plotting_params('riptrigspiking');
%     % create smoothing kernel
%     smoothing_length = 10;   % std of gaussian (in ms) used to smooth rip psth
%     smoothing_width = round(smoothing_length*.001/Fp.binsize);   % smoothing width in number of bins
%     kernel = gaussian(smoothing_width,smoothing_width*8);
%     for iAn = 1;%:length(animals);
%         animalinfo = animaldef(lower(Fp.animals{iAn}));
%         animalID = animalinfo{1,3}; %use anim prefix for name
%         FFanimdir =  sprintf('%s',animalinfo{1,2});
%         load([FFanimdir, animalID, 'tetinfo']);
%         tetinfoAll = cellfetch(tetinfo, '', 'alltags', 1);
%         dtcInds = cell2mat({ppF(:).dtc}');
%         [~, tagIndMap] = ismember(dtcInds(:,[1 2]),tetinfoAll.index(:,[1 3]), 'rows');
% %         ntrodeTags = tetinfoAll.values(tagIndMap);
% %         try
% %             numsumSupAreas = cellfun(@(x) sum(uint16(x.suparea)), ntrodeTags, 'UniformOutput', false);
% %             numsumAreas = cellfun(@(x) sum(uint16(x.area)), ntrodeTags, 'UniformOutput', false);
% %             numsumSubAreas = cellfun(@(x) sum(uint16(x.subarea)), ntrodeTags, 'UniformOutput', false);
% %             strSupAreas = cellfun(@(x) x.suparea, ntrodeTags, 'UniformOutput', false);
% %             strAreas = cellfun(@(x) x.area, ntrodeTags, 'UniformOutput', false);
% %             strSubAreas = cellfun(@(x) x.subarea, ntrodeTags, 'UniformOutput', false);
% %         catch
% %             error('all ntrodes need to have a suparea, subarea, and area tag, even if blank')
% %         end
% %         icolors = colorPicker(colorSet, strAreas, strSubAreas);
%         %% ---------- plot per cell --------------------------------
%         for icell = 1:length(ppF);
%             % if you want to save figs but not pause them to take a look.. just keep them invisible.
%             % i think this makes it faster and returns keyboard/mouse control during saving
%             set(0,'DefaultFigureWindowStyle','docked')
%             if savefigs && ~pausefigs;
%                 ifig = figure('Visible','off','units','normalized');%,'position',Pp.position);
%             else
%                 ifig = figure('units','normalized');%,'position',Pp.position);
%             end
%             set(gcf,'color','white')
%             load(sprintf('%s%s%s.mat',FFanimdir, animalID, 'cellinfo'));
%             % Plotting options
%             %             mark_confirmed_sleep_ripples = 0;  % only use if "sleep" data structure was passed into the DFA
%             %             timecourse_colormap = 0;
%             %             plot_variance = 0;
%                         region = 'ca1';
%             idata = ppF(icell);
% %             try
%                 empt_num = find(cellfun('isempty', idata.epoch_envs));
%                 idata.epoch_envs(empt_num) = {''};                
%                 [iepenvs, ~, epenvsInds] = unique(idata.epoch_envs, ...
%                     'stable');
% 
% %             catch
% %                 fprintf('cell not defined for all epochs: %d %d %d. skipping \n',...
% %                     idata.dtc);
% %                 close all
% %                 continue
% %             end
%             %% ---------- subplot per epoch type --------------------------------
%             n = 0;
%             for iepenv = 1:size(iepenvs,2); %for each epoch type
%                 if isempty(iepenvs{iepenv})
% %                     fprintf('empty eps \n')
%                     continue
%                 end
%                 n = n+1;
%                 epienvIDnums = find(iepenv == epenvsInds);
%                 [idataienvInds, idataienvbounds] = ismember( ...
%                     idata.eventtags(:,1),epienvIDnums);
%                 %get indices into the raster of the first event for each epoch
%                 [~, epbounds, ~] = unique(idataienvbounds(idataienvbounds > 0));
%                 ienvipsth = full(idata.psth(idataienvInds,:));
%                 sfA = n;
%                 %                 sfB = iepenv + size(iepenvs,2);
%                                 % get rid of empty eptype
%                 num_nonempty_eptypes = numel(find(~cellfun('isempty',iepenvs)));
%                 sfC = n + num_nonempty_eptypes*2;
%                 subaxis(4,num_nonempty_eptypes,[sfA sfC], 'Spacing',...
%                     Pp.Spacing, ...
%                     'SpacingVert', Pp.SpacingVert, 'SpacingHoriz', ...
%                     Pp.SpacingHoriz, 'Padding', Pp.Padding, 'MarginLeft', ...
%                     Pp.MarginLeft, 'MarginRight', Pp.MarginRight,...
%                     'MarginTop', Pp.MarginTop, 'MarginBottom', Pp.MarginBottom);
%                 epbounds = epbounds(epbounds ~= 1); %get rid of ep bnd at row 1
%                 for iepbnd = 1:size(epbounds,1)
%                     lh = line([idata.time(1)*1000 idata.time(end)*1000],...
%                         [epbounds(iepbnd) epbounds(iepbnd)],'LineWidth', ...
%                         Pp.lineWidth);
%                     lh.Color= [Pp.lineColor, Pp.lineAlpha];
%                 end
%                 hold on
%                 %% plot raster   %%%%%%%%%%%%%%%%%%%
%                 timevec = idata.time*1000;  % convert to ms
%                 % obtain bincenter form
%                 bincenters = timevec(1:(end-1)) + (timevec(2)-timevec(1)) / 2 ;
%                 % also truncate last bin of histc psth output
%                 psth = ienvipsth(:,1:(end-1));
%                 xpoints = psth.*repmat(bincenters,size(psth,1),1);
%                 xpoints(xpoints==0)=[];
%                 xpoints = xpoints(:);
%                 levelvec = (1:size(psth,1))';
%                 levelmat = repmat(levelvec,1,length(bincenters));
%                 ypoints = psth.*levelmat;
%                 ypoints(ypoints==0)=[];
%                 ypoints = ypoints(:);
%                 scatter(xpoints,ypoints,Pp.markersize,'MarkerEdgeColor','none',...
%                     'MarkerFaceColor',Pp.MarkerFaceColor,'MarkerFaceAlpha', ...
%                     Pp.MarkerFaceAlpha)
%                 hold on
%                 postripplematrix = idata.posteventmatrix(idataienvInds,:);
%                 ripmatrix = [ zeros(size(postripplematrix))  ones(size(postripplematrix,1),1)  postripplematrix ];
%                 psth = psth.*ripmatrix;
%                 ypoints2 = ripmatrix.*levelmat;
%                 plotpsth = psth.*ripmatrix;
%                 xpoints = plotpsth.*repmat(bincenters,size(plotpsth,1),1);
%                 xpoints(xpoints==0)=[];
%                 xpoints = xpoints(:);
%                 levelvec = (1:size(plotpsth,1))';
%                 levelmat = repmat(levelvec,1,length(bincenters));
%                 ypoints = plotpsth.*levelmat;
%                 ypoints(ypoints==0)=[];
%                 ypoints = ypoints(:);
%                 scatter(xpoints,ypoints,Pp.markersize,'MarkerEdgeColor', ...
%                     Pp.ripMarkerEdgeColor,'MarkerEdgeAlpha', ...
%                 Pp.ripMarkerEdgeAlpha,'MarkerFaceColor','none')    % highlight ripple-associated spikes
%                 %                 dff_plotraster(idata.time,ienvipsth,'markersize', 70, 'color', [.8 .8 .8], 'burstisi',6,'postripplematrix',idata.posteventmatrix(idataienvInds,:))
%                 hold on
%                 % axes
%                 xlim([bincenters(1) bincenters(end)])
%                 ylim([0 size(ienvipsth,1)])
%                 set(gca,'YDir','reverse')
%                 % x ticks
%                 %set(gca,'XTick',mean([timevec(1) timevec(2)]):500:mean([timevec(end-1) timevec(end)]))
%                 set(gca,'XTick',[])
%                 % label axes
%                 set(gca,'FontSize',8,'FontWeight','bold')
%                 %xlabel('time (ms)','FontSize',16,'FontWeight','bold')
%                 %                 ylabel('ripple #','FontSize',10,'FontWeight','bold')
%                 if n == 1;
%                     sylab = {['\color[rgb]{.5 .5 .5} \fontsize{12}', 'ripple #']};
%                     ylabel(sylab)
% %                     subylab = text(0, .5,sylab, 'FontSize',12,'FontWeight','bold','FontName', 'Arial', 'rotation', 90,...
% %                         'Units', 'normalized', 'horizontalAlignment', 'center');
%                 else
%                     set(gca,'ylabel',[])
%                 end
%                 %                 set(gcf, 'renderer', 'zbuffer')
%                 subtit = title(sprintf('%s', iepenvs{iepenv}), 'Color', ...
%                     Pp.subtitlecolor, 'FontSize',12,'FontWeight','bold',...
%                     'FontName', 'Arial');
%                 ylimy = get(gca, 'ylim');
%                 slh = line([0 0], ylimy, 'LineStyle', '-', 'LineWidth',1.5);
%                 slh.Color= [[0.5 0.5 0.5], .3];
%                 %% Create smoothed PSTH of mean firing rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 smoothedpsth = smoothvect(sum(ienvipsth,1)./(Fp.binsize*idata.noevents),kernel);
%                 %                 instFR = idata.instantFR(idataienvInds,:);
%                 %                 instFRstd = std(instFR,0,1);
%                 %                 instFRmean = mean(instFR,1);
%                 sfD = n + num_nonempty_eptypes*3;
%                 subaxis(4,num_nonempty_eptypes,sfD, 'Spacing', Pp.Spacing,...
%                     'SpacingVert', Pp.SpacingVert, 'SpacingHoriz', ...
%                     Pp.SpacingHoriz, 'Padding', Pp.Padding, 'MarginLeft', ...
%                     Pp.MarginLeft, 'MarginRight', Pp.MarginRight,...
%                     'MarginTop', Pp.MarginTop, 'MarginBottom', Pp.MarginBottom);
%                 ah = area(idata.time,zscore(smoothedpsth),'facecolor',...
%                     Pp.areaFaceColor);
%                 slh = line([idata.time(1) idata.time(end)], [2 2], ...
%                     'LineStyle', '--', 'LineWidth',2);
%                 slh.Color= [Pp.siglineColor, Pp.siglineAlpha];
%                 slh = line([idata.time(1) idata.time(end)], [-2 -2], ...
%                     'LineStyle', '--', 'LineWidth',2);
%                 slh.Color= [Pp.siglineColor, Pp.siglineAlpha];
%                 slh = line([0 0], [-3 3], 'LineStyle', '-', 'LineWidth',1.5);
%                 slh.Color= [[0.5 0.5 0.5], .3];
%                 %                 plot(idata.time,smoothvect(instFRmean,kernel),'color',areaFaceColor);
%                 %                 hold on
%                 %                 plot(idata.time,smoothvect(instFRmean,kernel)+instFRstd,'color',areaFaceColor);
%                 %                 hold on
%                 %                 errfill = fill([idata.time fliplr(idata.time)],[instFRmean'+instFRstd; instFRmean'-flipud(instFRstd)],[0 0 1],'linestyle','none');
%                 set(ah, 'FaceAlpha', Pp.areaAlpha)
%                 axis tight
%                 ylim([-3 3])
%                 set(gca,'FontSize',8,'FontWeight','bold')
%                 if n == 1;
%                 sylab = {['\color[rgb]{.5 .5 .5} \fontsize{12}', ...
%                     'zscore spikerate (std)']};
%                 ylabel(sylab)
% %                 subylab = text(-.2, .5,sylab, 'FontSize',12,'FontWeight','bold','FontName', 'Arial', 'rotation', 90,...
% %                     'Units', 'normalized', 'horizontalAlignment', 'center');
%                 else
%                     set(gca,'YTick',[])
%                 end
%                 hold on
%                 
% %                 % (Optional) plot variance - should look about the same as psth outline
% %                 %smoothing_width = round(smoothing_length*.001/binsize);   % smoothing width in number of bins
% %                 %kernel = gaussian(smoothing_width,smoothing_width*8);
% %                 if plot_variance
% %                     subplot(8,2,[14 15])
% %                     psth_variance_smoothed = smoothvect(var(full(idata.psth),1),kernel);
% %                     var_mean = mean(psth_variance_smoothed);
% %                     y_shift = ymaximum/2;
% %                     y_scale = 0.5 * ymaximum/range(psth_variance_smoothed);
% %                     psth_variance_toplot = (psth_variance_smoothed - var_mean) * y_scale + y_shift;
% %                     h = plot(idata.time,psth_variance_smoothed,'r','linewidth',3);
% %                     %         set(gca,'Color',bgcolor)
% %                     %         axis tight
% %                 end
%             end
%             
%             %% ---- super title and colorbar----
%             sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
%             sprtit = sprintf('%s - %d %d %d', animalID, idata.dtc(1),...
%                 idata.dtc(2), idata.dtc(3));
%             iStitle = text(.5, .95, {sprtit}, 'Parent', sprtitleax, 'Units', ...
%                 'normalized');
%             set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
%                 'horizontalAlignment', 'center');
% %             sprtit = sprintf('%s %s D%d nT%d C%d', Fp.filtfunction, animalID, ...
% %                 idata.dtc(1), idata.dtc(2), idata.dtc(3));
% %             iclr = icolors(icell,:);
% %             iarea = char(cell2mat(strAreas(icell))); %grab the area tags
% %             isubarea = num2str(cell2mat(strSubAreas(icell)));
% %             sprTags = sprintf('%s %s', iarea, isubarea);
% %             %             set(gca,'position',posx1)
% %             iStitle = text(.5, .95, [{sprtit};...
% %                 {['\fontsize{12} \color[rgb]' sprintf('{%d %d %d} %s ', iclr, sprTags)]};...
% %                 {['\color[rgb]{.5 .5 .5} \fontsize{8}', sprintf(' {%s}', filenameTitle)]}], 'Parent', sprtitleax, 'Units', 'normalized');
% %             set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'horizontalAlignment', 'center');
%             %% ---- pause, save figs ----
%             if runAntiAlias
%                 myaa; %will drastically improve jagginess but slow things down quite a bit. only use for print figure making
%                 close(ifig)
%             end
%             if pausefigs
%                 pause
%             end
%             if savefigs
%                 save_figure(paths.figdirectory, paths.filenamesave, sprtit)
%             end
%             
%             close all
%         end
%     end
% end
% 
% 
% %%  Calculate significance of modulation
% 
% close all
% 
% if calcstatsPEM
%     
%     % option to save ripplemod struct, below
%     save_ripplemod = 0;
%     daylist = unique(dtc(:,1))';
%     tic
%     nshuffles = 1000;
%     % set up ripplemod struct for storing p values - currently saves all cells and days in one struct
%     ripplemod = struct;
%     ripplemod.days = daylist;
%     ripplemod.event = eventconsname;
%     ripplemod.celltype = cellfilter;
%     ripplemod.epochtype = epochtype;
%     
%     % previously had an option to save this for each day, a la "eliripplemod04"
%     % ^^ do we want this??
%     %       for day = daylist
%     %       tmpdtc =  dtc(dtc(:,1)==day,:);
%     
%     % Modvalues = [day, tet, cell, direction(pos=+1, neg=-1), pvalue, depth-of-modulation, mean-FR-in-window]
%     ripplemod.modvalues = [dtc repmat(zeros(size(dtc,1),1),1,4)];
%     
%     for cellind = 1:size(dtc,1)      % (!!) if doing this for cherry-picked cells, enter cellind here
%         %             if dtc(cellind,1) == day
%         analyzing_cell = dtc(cellind,:);
%         disp(['analyzing cell: ' num2str(analyzing_cell)])
%         
%         spikedata = full(ppF(cellind).psth);
%         numspikes = sum(full(ppF(cellind).psthsum),2);
%         
%         % mean FR in the plotted window
%         windowFR =  numspikes/((length(ppF(cellind).psthsum)*binsize)*size(ppF(cellind).psth,1));
%         
%         % (Optional) Randomly delete a certain proportion of spikes for each cell
%         % to match a control firing rate, specified above.  "sparsify" must equal 1
%         control = 0.5; %control FR to match
%         if windowFR > control && sparsify
%             spikeinds = find(spikedata);
%             flip = rand(size(spikeinds));
%             deletespikes = flip<=0.5;
%             spikedata(spikeinds(deletespikes))=0;
%             windowFR = sum(sum(spikedata,1),2)/((length(ppF(cellind).psthsum)*binsize)*size(ppF(cellind).psth,1))
%             sparsified = 1;
%             disp('sparsified')
%         end
%         
%         % time range to analyze significance... could also use the average event duration, on either side of 0
%         range = [0 0.2]; %mean(A(cellind).eventduration);
%         % varRange = index for range(1) to range(2) in time series
%         varRange = [lookup(range(1),ppF(cellind).time):lookup(range(2),ppF(cellind).time)];
%         
%         % set windowisrange to 1 to only shuffle spikes within the analysis range, not the whole -0.5 to 0.5
%         windowisrange = 0;
%         if windowisrange
%             tmpspikedata=spikedata(:,varRange); %for shuffling only within the analysis range
%             spikedata = tmpspikedata;
%         end
%         
%         binwindow = size(spikedata,2);
%         halfbinwindow = size(spikedata,2)/2;
%         shufStep = [-halfbinwindow halfbinwindow]; % maximum time to shuffle spikes
%         
%         allShufPsth = zeros(nshuffles,binwindow);
%         
%         % smoothing should be the same as the plotted real data
%         smoothing_length = 10;   % std of gaussian (in ms) used to smooth rip psth
%         smoothing_width = round(smoothing_length*.001/F.runscript_params.binsize);   % smoothing width in number of bins
%         kernel = gaussian(smoothing_width,smoothing_width*8);
%         
%         % circularly shuffle individual ripple-triggered spiketrains by r # of bins,
%         for n = 1:nshuffles
%             currshuffle = zeros(size(spikedata));
%             for rip = 1:size(spikedata,1) % number of rips depends on number of clustered epochs for the cell
%                 %create an empty "spiketrain"
%                 dummyspikes = zeros(1,size(spikedata,2));
%                 %generate random value from the uniform distribution -0.5 to .5 (peri-ripple window in bins)
%                 r = round(shufStep(1) + (shufStep(2)-shufStep(1)).*rand(1,1));
%                 %find indices of spiketimes in current train
%                 spikeinds = find(spikedata(rip,:));
%                 shiftedspikeinds = spikeinds+r;
%                 % if times fall off the beginning or end, bring them around circularly
%                 circindend = find(shiftedspikeinds>size(spikedata,2));
%                 shiftedspikeinds(circindend)=shiftedspikeinds(circindend)-size(spikedata,2);
%                 circindstart = find(shiftedspikeinds<1);
%                 shiftedspikeinds(circindstart)=size(spikedata,2)-(1-shiftedspikeinds(circindstart));
%                 %fill in dummy spike train
%                 dummyspikes(shiftedspikeinds) = 1;
%                 currshuffle(rip,:) = dummyspikes;
%             end
%             % smoothed the shuffled psth
%             smoothedcurrpsth = smoothvect(sum(currshuffle,1)./(binsize*size(spikedata,1)),kernel);
%             %      figure
%             %      h = bar(1:length(smoothedcurrpsth),smoothedcurrpsth,'facecolor',[0 0 0]);
%             allShufPsth(n,:) = smoothedcurrpsth;
%         end
%         
%         %take mean of shuffles
%         meanShuf = mean(allShufPsth,1);
%         
%         %if shuffle was only done within the analysis range, use the whole thing
%         if windowisrange
%             meanShufRange = meanShuf;
%         else
%             % otherwise we analyze only the shuffled psth in "range"
%             meanShufRange = meanShuf(varRange);
%         end
%         % figure
%         % h = bar(1:length(meanShuf),meanShuf,'facecolor',[0 0 0]);
%         
%         %find summed squared distance of each shuffle
%         meanvarShuf = zeros(nshuffles,1);
%         for n=1:nshuffles
%             if windowisrange
%                 allShufRange = allShufPsth(n,:);
%             else
%                 allShufRange = allShufPsth(n,varRange);
%             end
%             meanvarShuf(n,:) = sum(((allShufRange-meanShufRange).^2),2);
%         end
%         
%         %find summed squared distance of the real data psth from the mean of the shuffles, within the analysis range
%         RealPsth = smoothvect(sum(spikedata,1)./(binsize*size(spikedata,1)),kernel);
%         if windowisrange
%             RealPsthRange = RealPsth;
%         else
%             RealPsthRange = RealPsth(varRange);
%         end
%         meanvarReal = sum((RealPsthRange-meanShufRange).^2);
%         
%         if sparsified
%             figure
%             h = bar(1:length(RealPsth),RealPsth,'facecolor',[0.5 0.5 0.5]);
%             % figure
%             % h = bar(varRange,RealPsthRange,'facecolor','b');
%             % figure
%             % h1 = bar(varRange,meanShufRange,'facecolor','k');
%         end
%         
%         % p value is 1- (what fraction of shuffles have a variance that exceeds the real data)
%         p = 1 - sum(meanvarShuf<meanvarReal)/nshuffles
%         
%         % direction of modulation
%         if  sum(RealPsthRange-meanShufRange,2)<0
%             direction = -1;
%         else
%             direction = 1;
%         end
%         
%         ripplemod.modvalues(cellind,4) = direction;
%         ripplemod.modvalues(cellind,5) = p;
%         ripplemod.modvalues(cellind,6) = meanvarReal; %depth of modulation
%         ripplemod.modvalues(cellind,7) = windowFR;
%         
%         %save ripplemod variable
%         if save_ripplemod
%             modname = sprintf('%sripplemoddata/%sripplemod_%s_%s_%s.mat',filtOutputDirectory,anim,epochtype,eventconsname,date);
%             save(modname,'ripplemod','-mat')
%         end
%         
%     end
%     
%     %         end
%     ripplemod.modvalues
%     %     end
%     toc
% end
% end
