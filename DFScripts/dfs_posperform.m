
close all

% loadBehaveStructs = 1;
% calculateStateSpace = 1;  
% saveStateSpaceResults = 1;
plotPosPerform = 1;
savefigs = 1;
pausefigs = 0;
runAntiAlias = 0;

% savePerformResults = 0;
% plotPerform = 0;


%     MarginRight = 0.04;
%     MarginTop = 0.09;
%     MarginBottom =  0.08;

%% ---------------- plotting params --------------------------
colorSet = 'DR1';
% clims = [0 1]; %[0 .7]
% position = [.1 .1 .9 .8];
SpacingHorizontal = 0.1;
SpacingVertical = 0.0;
position = [.1 .1 .8 .6];
SpacingHorizontal = 0.00;
SpacingVertical = 0.05;
Spacing = 0.01;
Padding = 0.00;
MarginLeft = 0.04;
MarginRight = 0.04;
MarginTop = 0.09;
MarginBottom =  0.08;
usecolormap = 'jet';
win = [.5 .5]; %in seconds
% indwin = win*1500;
%% ---------------- Data Filters --------------------------
animal = 'D10';
% animals = {'JZ1', 'D13'};
days = [1]; %1:10
behavestruct = 'BehaveState';
plottype = 'posperform';
eventtype = 'rippleskons';
epochEnvironment = 'wtrack';% 'wtrack'; %wtrack, wtrackrotated, openfield, sleep
eventSourceArea = 'ca1';
ripAreas = {'ca1'}; %, 'mec', 'por'
% ntAreas = {'ca1', 'sub', 'mec', 'por', 'v2l', 'ref'};
consensus_numtets = 1;   % minimum # of tets for consensus event detection
minstdthresh = 2;        % STD. how big your ripples are
% exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
% minvelocity = 0;  
maxvelocity = 4;
% outputDirectory = '/typhoon/droumis/analysis';
investInfo = animaldef(lower('Demetris'));
figdirectory = sprintf('%s%s/', investInfo{4}, plottype);
animalinfo = animaldef(lower(animal));
animalID = animalinfo{1,3}; %use anim prefix for name
FFanimdir =  sprintf('%s',animalinfo{1,2});
%% ---------------- Paths and Title strings ---------------------------------------------------
currfigdirectory = sprintf('%s%s/',figdirectory);
% filenamesave = sprintf('%s%s_%s_%s_%s', eventSourceArea, eventtype, epochEnvironment, cell2mat(animals));
% filename = sprintf('%s_%s.mat', behavestruct, filenamesave);
% resultsOutputDirectory = sprintf('%s%s/', investInfo{3}, behavestruct);
filename = sprintf('%s_%s_%s_%s', plottype, eventtype, epochEnvironment, animal);
filenameTitle = strrep(filename,'_', ' ');
% %% ---------------- Load BehaveStruct---------------------------------------------------
% if loadBehaveStructs;
%
% end

%% ---------------- plot posperform---------------------------------------------------
if plotPosPerform
    
    load(sprintf('%s%s%s.mat',FFanimdir, animalID, behavestruct));
    wTrackBehave = cellfetch(BehaveState.statechanges, 'wtrack', 'alltags', 1);
    allepsMat = []; eplengths = [];
    %     days = unique(wTrackBehave.index(:,1),'stable');
    for day = days;
        eps = unique(wTrackBehave.index(wTrackBehave.index(:,1)==day,2),'stable');
        if savefigs && ~pausefigs; %if you want to save figs but not pause them to take a look.. just keep them invisible. i think this makes it faster and returns keyboard/mouse control during saving
            ifig = figure('Visible','off','units','normalized', 'position',position);
        else
            ifig = figure('units','normalized', 'position',position);
            %             set(0,'DefaultFigureWindowStyle','docked') %dock the fig
            %             set(0,'DefaultFigureWindowStyle','normal')
        end
        set(gcf,'color','white');
        %         load(sprintf('%s%s%s%02d.mat',FFanimdir, animalID, 'task', day));
            load(sprintf('%s%s%s%02d.mat',FFanimdir, animalID, 'linpos', day));
        %         load(sprintf('%s%s%s%02d.mat',FFanimdir, animalID, 'ca1rippleskons', day));
        %         load(sprintf('%s%s%s%02d.mat',FFanimdir, animalID, 'mecrippleskons', day));
        %         load(sprintf('%s%s%s%02d.mat',FFanimdir, animalID, 'porrippleskons', day));
        for iep = 1:length(eps);
            ep = eps(iep);
            isfig = subaxis(length(eps),1,iep,'SpacingVertical', SpacingVertical, 'SpacingHorizontal', SpacingHorizontal, ...
                'Padding', Padding, 'MarginLeft', MarginLeft, 'MarginRight', MarginRight, 'MarginTop', MarginTop, ...
                'MarginBottom', MarginBottom);
            %% get the ripples from each area
            clear ripout
            strSubAreas = {'d' '5' '6'}; %hack for dorsal ca1, mec layer 5, and por layer 6 colors
            icolors = colorPicker(colorSet, ripAreas, strSubAreas);
            cmap = jet(101);
            for iarea = 1:length(ripAreas)
                ripout{iarea} = DR_kk_getconstimes(animalinfo{2},animalinfo{3},[day ep],sprintf('%s%s',ripAreas{iarea}, eventtype),1,...
                    'consensus_numtets', consensus_numtets, 'maxvelocity',maxvelocity, 'minstdthresh', minstdthresh);
                % kk_getconstimes(animaldir,animalprefix, epochs, eventconsname, tetfilter, varargin)
                % ripout = kk_getconstimes(animalinfo{2},animalinfo{3},[day ep], 'ripplesdons',1,'consensus_numtets',3);
                %	maxthresh - the largest threshold in stdev units at which this ripple would still be detected.
                maxstd_valsVec = list2vec([ripout{iarea}{day}{ep}.starttime ripout{iarea}{day}{ep}.endtime], ripout{iarea}{day}{ep}.time, 'marker', ripout{iarea}{day}{ep}.maxthresh);
                [periodtimes_rip2 periodripsinds] = dr_vec2list(ripout{iarea}{day}{ep}.cons,ripout{iarea}{day}{ep}.time);
                maxstd_vals = maxstd_valsVec(periodripsinds(:,1));
                maxstd_valsNormed = round((maxstd_vals-min(maxstd_vals))./(max(maxstd_vals) - min(maxstd_vals))*100);
                maxstd_valsNormed = maxstd_valsNormed+1;
                ripout{iarea}{day}{ep}.periodtimes_rip2 = periodtimes_rip2;
                ripout{iarea}{day}{ep}.periodripsinds = periodripsinds;
                irips = ripout{iarea}{day}{ep}.periodtimes_rip2;
%                 plot(repmat(irips(:,1),1,2)', repmat([-iarea*5;(-iarea-1)*5],1, size(irips,1)),'-', 'Color', icolors(iarea,:), 'LineWidth',1); hold on;
%                 set(isfig,'defaultAxesColorOrder',cmap(maxstd_valsNormed,:))
                %the colors are normalized max threshold std per area/epoch
                set(isfig,'ColorOrder',cmap(maxstd_valsNormed,:))
                line(repmat(irips(:,1),1,2)', repmat([-iarea*5;(-iarea-1)*5],1, size(irips,1)), 'LineWidth',2); hold on;
                
            end
            %% plot lindist and overlay all the rip times and trial out corr/mistake patches
            
            postime = linpos{day}{ep}.statematrix.time;
            traj = linpos{day}{ep}.statematrix.traj;
            segmentnum = linpos{day}{ep}.statematrix.segmentIndex;
            % plot animal's travel lindist from center well and colorize by segment
            lindist = linpos{day}{ep}.statematrix.linearDistanceToWells(:,1);
            plot(postime,lindist,'linewidth',1.5,'color',[.8 .8 .8]);
            hold on
            %         clr1 = [.5 .3 .8; .2 .1 .1; .4 .1 .6; .2 .9 .6];
            clr1 = [.6 .5 .2; .8 .5 .8; .6 0 .6; 0 0.4980 0.3255; 0 0.8 0.5216]; %.2 .9 .6];
            %             %plot ripple times from all areas
            for iseg = 1:5
                inds = [];
                inds = double(segmentnum == iseg);
                inds(find(inds==0))=nan;
                plot(postime.*inds,lindist.*inds,'-','linewidth',2,'color',clr1(iseg,:));
                hold on;
            end
            
            %plot correct/incorrect outbound patches
            clr2 = ([.5 .5 1; .5 .7 .7; .3 .8 .3; .5 .5 1; .5 .7 .7; .3 .8 .3;]);
            trialIO = BehaveState.statechanges{day}{ep}.statechangeseq;
            trialIOfields = BehaveState.statechanges{day}{ep}.fields;
            corrcol = find(cell2mat(cellfun(@(x) strcmp(x,'correct'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
            timeportoutcol = find(cell2mat(cellfun(@(x) strcmp(x,'timeportout'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
            outBcol = find(cell2mat(cellfun(@(x) strcmp(x,'outbound'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
            inBcol = find(cell2mat(cellfun(@(x) strcmp(x,'inbound'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
            lastTimecol = find(cell2mat(cellfun(@(x) strcmp(x,'lasttime'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
            currTimecol = find(cell2mat(cellfun(@(x) strcmp(x,'currenttime'), strsplit(trialIOfields, ' '), 'UniformOutput', false)));
            outBCorr = trialIO((trialIO(:,outBcol)==1 & trialIO(:,corrcol)==1),[lastTimecol currTimecol corrcol timeportoutcol]);
            outBMist = trialIO((trialIO(:,outBcol)==1 & trialIO(:,corrcol)==0),[lastTimecol currTimecol corrcol timeportoutcol]);
            if ~isempty(find(outBMist(:,3:4))) || ~isempty(find(~outBCorr(:,3:4)))
                error('behavestate struct indicates incorrect trial scoring')
            end
            yl = ylim;
            xl = xlim;
            for itrial = 1:length(outBCorr(:,1)); %patch correct outbound trials
                patch([outBCorr(itrial,1) outBCorr(itrial,2) outBCorr(itrial,2) outBCorr(itrial,1)],  [yl(2) yl(2) yl(1) yl(1)] , 'y', 'edgecolor','none', 'FaceAlpha',.2);
            end
            for itrial = 1:length(outBMist(:,1)); %patch incorrect outbound trials
                patch([outBMist(itrial,1) outBMist(itrial,2) outBMist(itrial,2) outBMist(itrial,1)],  [yl(2) yl(2) yl(1) yl(1)] , 'k', 'edgecolor','none', 'FaceAlpha',.1);
            end
            for itrial = 1:length(outBCorr(:,1)); %is a portOUT
                plot(outBCorr(itrial,4), yl(2)-1, 'd', 'linewidth', 1, 'MarkerEdgeColor', 'k','MarkerFaceColor', 'y'); hold on;
            end
            
            axis tight
            set(gca,'TickDir','out')
            a = gca;
            % set box property to off and remove background color
            set(a,'box','off','color','none')
            % create new, empty axes with box but without ticks
            b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
            % set original axes as active
            axes(a) %not sure why but this is negating the invisibility when saving figs
            % link axes in case of zooming
            linkaxes([a b])
            a = get(gca,'XTickLabel');
            set(gca,'XTickLabel',a,'FontName','Arial','fontsize',8)
            a = get(gca,'XTick');
            set(gca,'XTick',a,'FontName','Arial','fontsize',8)
            a = get(gca,'YTickLabel');
            set(gca,'YTickLabel',a,'FontName','Arial','fontsize',8)
            set(gca,'TickLength',[0.005,0.005]);
            %             fcclr = ([1 1 1 ; 0.2 .6 .9; 0 0 0]);
            %             edgclr = ([1 0 0; 0 0 1; 0 0 0]);
            %             lintyp = (['o' 'd' 'x']); % o * x
            %             linwid = ([2 3 1.5]);
            %             marksz = ([7 4 7]);
            %             for iarea = 1:length(ripAreas)
            %                 riptimesind = knnsearch(postime, ripout{iarea}{day}{ep}.periodtimes_rip2(:,1));
            %                 riptimesLIN = postime(riptimesind);
            %                 ripYLIN = lindist(riptimesind);
            %                 plot(riptimesLIN, ripYLIN,lintyp(iarea), 'MarkerSize',marksz(iarea), 'LineWidth',linwid(iarea), 'MarkerEdgeColor',edgclr(iarea,:), 'MarkerFaceColor',fcclr(iarea,:)); hold on;%'MarkerFaceColor',[.49 1 .63] 'LineWidth',1,
            %             end
%             if iep == 1
% %                 set(gca, 'XTick', []);
%             end
            title(sprintf('epoch %d', ep),'FontSize',8,'FontWeight','bold','Color',[.5 .5 .5], 'FontName', 'Arial')
        end
        
        sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
%             

                    clrbar = colorbar('location','eastoutside', 'FontSize',6,'FontName', 'Arial');%, 'FontWeight','bold');
                    colormap(cmap)
                    caxis([min(maxstd_vals) max(maxstd_vals)])
            posx1=get(gca,'position');
            posx=get(clrbar,'Position');
            posx(1)= 1-MarginRight+.01;
            posx(2)= MarginBottom;
            posx(3)= .01;
            posx(4)= 1 - MarginTop - MarginBottom;
            set(clrbar,'Position',posx)
            set(gca,'position',posx1)
            
        suptitclr = [0 0 0];
        suptittxt = sprintf('%s Day%d', filenameTitle, day);
        iStitle = text(.5, .98, ['\fontsize{12} \color[rgb]' sprintf('{%d %d %d} %s ', suptitclr, suptittxt)]);
        set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'horizontalAlignment', 'center');
        xlab = 'time (s)';
        ylab = 'linear distance from centerwell (cm)';
        supxlabel = text(.5, .02, xlab, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
            'Parent', sprtitleax, 'Units', 'normalized', 'horizontalAlignment', 'center');
        supylabel = text(.01, .5, ylab, 'FontSize',12,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
            'rotation', 90, 'Parent', sprtitleax, 'Units', 'normalized', 'horizontalAlignment', 'center');
%         for iarea = 1:length(ripAreas)
%                     arealab = text(.01, (1-(iarea/(length(ripAreas)+1)))*.1, ripAreas{iarea}, 'FontSize',8,'FontWeight','bold','Color',icolors(iarea,:), 'FontName', 'Arial', ...
%             'Parent', sprtitleax, 'Units', 'normalized', 'horizontalAlignment', 'center');
%         end
        clrbartit = text(posx(1)+posx(3)/2, posx(2)-MarginBottom/2, 'swr-std', 'FontSize',10,'FontWeight','bold','Color','k', 'FontName', 'Arial','horizontalAlignment', 'center');
            
        %% ---- pause, save figs ----
        if runAntiAlias
            myaa; %will drastically improve jagginess but slow things down quite a bit. only use for print figure making
            close(ifig)
        end
        if pausefigs
            pause
        end
        if savefigs
            if ~isdir(currfigdirectory);
                mkdir(currfigdirectory);
            end
            sprtitsave = strrep(suptittxt,' ', '_'); %put '_' character in place of whitespace for filename
            set(gcf, 'PaperPositionMode', 'auto'); %saves the png in the dimensions defined for the fig
            currfigfile = sprintf('%s%s',currfigdirectory, sprtitsave);
            print(currfigfile,'-dpng', '-r0')
            disp(sprintf('plot %s saved', suptittxt))
        end
        close all;
    end
end


% %% ---------------- plot perform---------------------------------------------------
% if plotPerform
%
% end
% %% ---------------- save PerformResults---------------------------------------------------
% if savePerformResults
% end
















































































