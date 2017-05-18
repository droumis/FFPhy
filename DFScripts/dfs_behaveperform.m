
close all

loadBehaveStruct = 1;
calculateStateSpace = 1;
saveStateSpaceResults = 1;
plotStateSpace = 0;
savefigs = 0;
pausefigs = 0;

% savePerformResults = 0;
% plotPerform = 0;

%% ---------------- plotting params --------------------------
colorSet = 'DR1';
% clims = [0 1]; %[0 .7]
% position = [.1 .1 .9 .8];
SpacingHorizontal = 0.01;
SpacingVertical = 0.02;
% Spacing = 0.00;
% Padding = 0.0;
% MarginLeft = 0.04;
% MarginRight = 0.04;
% MarginTop = 0.09;
% MarginBottom =  0.08;
position = [.1 .1 .9 .8];
SpacingHorizontal = 0.00;
SpacingVertical = 0.00;
Spacing = 0.00;
Padding = 0.00;
MarginLeft = 0.05;
MarginRight = 0.05;
MarginTop = 0.14;
MarginBottom =  0.08;
usecolormap = 'jet';
win = [.5 .5]; %in seconds
% indwin = win*1500;
%% ---------------- Data Filters --------------------------
animal = 'JZ1';
% animals = {'JZ1', 'D13'};
days = [1];
behavestruct = 'BehaveState';
eventtype = 'rippleskons';
epochEnvironment = 'wtrack';% 'wtrack'; %wtrack, wtrackrotated, openfield, sleep
eventSourceArea = 'ca1';
% ripAreas = {'ca1', 'mec', 'por'};
% ntAreas = {'ca1', 'sub', 'mec', 'por', 'v2l', 'ref'};
% consensus_numtets = 1;   % minimum # of tets for consensus event detection
% minstdthresh = 5;        % STD. how big your ripples are
% exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
% minvelocity = 0;
% maxvelocity = 4;
% outputDirectory = '/typhoon/droumis/analysis';
investInfo = animaldef(lower('Demetris'));
figdirectory = sprintf('%s%s/', investInfo{4}, behavestruct);
animalinfo = animaldef(lower(animal));
animalID = animalinfo{1,3}; %use anim prefix for name
FFanimdir =  sprintf('%s',animalinfo{1,2});
%% ---------------- Paths and Title strings ---------------------------------------------------
currfigdirectory = sprintf('%s%s/',figdirectory);
% filenamesave = sprintf('%s%s_%s_%s_%s', eventSourceArea, eventtype, epochEnvironment, cell2mat(animals));
% filename = sprintf('%s_%s.mat', behavestruct, filenamesave);
resultsOutputDirectory = sprintf('%s%s/', investInfo{3}, behavestruct);
resultfilename = sprintf('%s_%s_%s', behavestruct, epochEnvironment, animal);
filenameTitle = strrep(resultfilename,'_', ' ');
%% ---------------- Load BehaveStruct---------------------------------------------------
if loadBehaveStruct;
    load(sprintf('%s%s%s.mat',FFanimdir, animalID, behavestruct));
end
%% ---------------- calculate StateSpace Model---------------------------------------------------
wTrackBehave = cellfetch(BehaveState.statechanges, 'wtrack', 'alltags', 1);
allepsMat = []; dayepinds = [];
for iep = 1:length(wTrackBehave.index(:,1))
    allepsCell{iep} = BehaveState.statechanges{wTrackBehave.index(iep,1)}{wTrackBehave.index(iep,2)}.statechangeseq;
    allepsCell{iep} = allepsCell{iep}(2:end,:); %get rid of the first state end change time as it doesn't have a an accurate start time
    allepsMat = [allepsMat; allepsCell{iep}];
    dayepinds = [dayepinds; repmat(wTrackBehave.index(iep,:), length(allepsCell{iep}(:,1)), 1)];
end
allepsMatFields = BehaveState.statechanges{wTrackBehave.index(iep,1)}{wTrackBehave.index(iep,2)}.fields;
chance = 0.5;
if calculateStateSpace
    % [probcorrect lastpoint] = getestprobcorrect(behavperform, backgroundprob, initialcond, doplot)
    %   probcorrect
    %       The first column of probcorrect is the mode
    %       The second column of probcorrect is the lower 5% confidence bound
    %       The third column of probcorrect is the upper 5% confidence bound
    %       DR04/30/17  Fourth col is certainty matrix
    %       first row gets exluded
%     tic
    [pcALL] = getestprobcorrect(allepsMat(:,7),chance,0,0);
    if length(pcALL(2:end,1)) ~= length(allepsMat(:,1));
        error('length of state space results do not match trial input length')
    end
%     alltoc = toc;
%     tic
    [pcINB] = getestprobcorrect(allepsMat(allepsMat(:,8)==1,7),chance,0,0);
%     intoc = toc;
%     tic
    [pcOUTB] = getestprobcorrect(allepsMat(allepsMat(:,9)==1,7),chance,0,0);
%     outtoc = toc;
    statespace.allbound = [pcALL(2:end,:) dayepinds];
    statespace.inbound = [pcINB(2:end,:) dayepinds(allepsMat(:,8)==1,:)];
    statespace.outbound = [pcOUTB(2:end,:) dayepinds(allepsMat(:,9)==1,:)];
    statespace.pcFields = 'mode lower5 upper5 certainty day epoch';
    statespace.behavestruct = behavestruct;
    statespace.epochEnvironment = epochEnvironment;
    statespace.allepsMat = allepsMat;
    statespace.allepsMatFields = allepsMatFields; 
end

%% ---------------- save StateSpaceResults into BehaveState struct---------------------------------------------------
if saveStateSpaceResults == 1;
%     if ~isdir(resultsOutputDirectory);
%         mkdir(resultsOutputDirectory);
%     end
    %save the entire workspace for filter provenance
    BehaveState.statespace = statespace;
    save(sprintf('%s%s%s.mat',FFanimdir, animalID, behavestruct), 'BehaveState','-v7.3');
    disp(sprintf('%s%s%s.mat Saved with StateSpace',FFanimdir, animalID, behavestruct));
end
% %% ---------------- load StateSpaceResults---------------------------------------------------
% if loadStateSpaceResults == 1;
%     load(sprintf('%s%sStateSpace_%s.mat',resultsOutputDirectory,animal, epochEnvironment));
%     disp(sprintf('%s%sStateSpace_%s.mat Loaded',resultsOutputDirectory,animal, epochEnvironment))
% end

%% ---------------- plot statespace---------------------------------------------------
if plotStateSpace
    try
        statespace = BehaveState.statespace;
    catch
%         if the statespace wasnt yet saved into the behavstate studt then jst use the one generated in this session
    end
    for iseq = 1:3; % 1 Allbound; 2 Inbound; 3 Outbound
        switch iseq
            case 1
                pc = statespace.allbound;
                corrstring = 'correct';
                corrcol = find(cell2mat(cellfun(@(x) strcmp(x,corrstring), strsplit(statespace.allepsMatFields, ' '), 'UniformOutput', false)));
                behavperform = allepsMat(:,corrcol);
                trialType = 'Allbound';
            case 2
                pc = statespace.inbound;
                corrstring = 'correct';
                corrcol = find(cell2mat(cellfun(@(x) strcmp(x,corrstring), strsplit(statespace.allepsMatFields, ' '), 'UniformOutput', false)));
                trialType = 'Inbound';
                trialtypestring = lower(trialType);
                trialtypecol = find(cell2mat(cellfun(@(x) strcmp(x,trialtypestring), strsplit(statespace.allepsMatFields, ' '), 'UniformOutput', false)));
                behavperform = allepsMat(allepsMat(:,trialtypecol)==1,corrcol);
            case 3
                pc = statespace.outbound;
                corrstring = 'correct';
                corrcol = find(cell2mat(cellfun(@(x) strcmp(x,corrstring), strsplit(statespace.allepsMatFields, ' '), 'UniformOutput', false)));
                trialType = 'Outbound';
                trialtypestring = lower(trialType);
                trialtypecol = find(cell2mat(cellfun(@(x) strcmp(x,trialtypestring), strsplit(statespace.allepsMatFields, ' '), 'UniformOutput', false)));
                behavperform = allepsMat(allepsMat(:,trialtypecol)==1,corrcol);
        end
        modecol = find(cell2mat(cellfun(@(x) strcmp(x,'mode'), strsplit(statespace.pcFields, ' '), 'UniformOutput', false)));
        lower5col = find(cell2mat(cellfun(@(x) strcmp(x,'lower5'), strsplit(statespace.pcFields, ' '), 'UniformOutput', false)));
        upper5col = find(cell2mat(cellfun(@(x) strcmp(x,'upper5'), strsplit(statespace.pcFields, ' '), 'UniformOutput', false)));
        certaintycol = find(cell2mat(cellfun(@(x) strcmp(x,'certainty'), strsplit(statespace.pcFields, ' '), 'UniformOutput', false)));
        %plot statespace results
        tALL=1:size(pc,1)-1;
        if savefigs && ~pausefigs; %if you want to save figs but not pause them to take a look.. just keep them invisible. i think this makes it faster and returns keyboard/mouse control during saving
            ssFig = figure('Visible','off','units','normalized');
        else
            ssFig = figure('units','normalized');
        end
        set(gcf,'color','white');
        subplot(2,1,1);
        hold on;
        line([cumsum(eplengths); cumsum(eplengths)], [0 1], 'Color', [.9 .9 .9], 'LineStyle', '-', 'LineWidth', 1)
        plot(tALL, pc(2:end,modecol),'b-', 'LineWidth', 2); %plot behavior SS score
        errfillAll = fill([tALL fliplr(tALL)],[pc(2:end,lower5col); flipud(pc(2:end,upper5col))],[0 0 1],'linestyle','none');
        set(errfillAll, 'FaceAlpha', .2)
        hold on; [x, y] = find(behavperform > 0);
        h = plot(x,y-0.03,'+'); set(h, 'MarkerFaceColor','none');
        set(h, 'MarkerEdgeColor', [.2 .7 .2]);
        hold on; [x, y] = find(behavperform == 0);
        h = plot(x,y-0.05,'+'); set(h, 'MarkerFaceColor', 'none');
        set(h, 'MarkerEdgeColor', [.5 .5 .5]);
        axis([1 tALL(end)  0 1.05]);
        line([1 tALL(end)], [chance  chance], 'Color', [.5 .5 .5], 'LineStyle', '--');
        xlabel('Trial Number')
        ylabel([{'Probability of a'};{'Correct Response'}])
        
        subplot(2,1,2)
        hold on;
        line([cumsum(eplengths); cumsum(eplengths)], [0 1], 'Color', [.9 .9 .9], 'LineStyle', '-', 'LineWidth', 1)
        plot(tALL,1 - pc(2:end,certaintycol),'k', 'LineWidth', 2)
        line([ 1 tALL(end)],[0.90 0.90], 'LineStyle', '-');
        line([ 1 tALL(end)],[0.99 0.99], 'LineStyle', '--');
        line([ 1 tALL(end)],[0.95 0.95], 'LineStyle', '-.');
        axis([1 tALL(end)  0 1]);
        grid off;
        xlabel('Trial Number')
        ylabel('Certainty')
        
        sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ssFig);
        suptitclr = [0 0 0];
        suptittxt = sprintf('%s %s', filenameTitle, trialType);
        iStitle = text(.5, .95, ['\fontsize{10} \color[rgb]' sprintf('{%d %d %d} %s ', suptitclr, suptittxt)]);
        set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', 'horizontalAlignment', 'center');
        %% ---- pause, save figs ----
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
%             %     load(sprintf('%s%s%s%02d.mat',FFanimdir, animalID, 'task', day));
%         %     load(sprintf('%s%s%s%02d.mat',FFanimdir, animalID, 'linpos', day));
%         %     load(sprintf('%s%s%s%02d.mat',FFanimdir, animalID, 'ripplekonsca1', day));
%         %     load(sprintf('%s%s%s%02d.mat',FFanimdir, animalID, 'ripplekonsmec', day));
%         %     load(sprintf('%s%s%s%02d.mat',FFanimdir, animalID, 'ripplekonspor', day));
% end
% %% ---------------- save PerformResults---------------------------------------------------
% if savePerformResults
% end
