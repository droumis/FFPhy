
close all

loadBehaveStruct = 0;
calculateStateSpace = 1;
filterpureWdays = 0;
saveStateSpaceResults = 0;

plotStateSpace = 1;
savefigs = 1;
pausefigs = 0;

% savePerformResults = 0;
% plotPerform = 0;
%% ---------------- Data Filters --------------------------
% animal = 'D10';
% behavestruct = 'BehaveState';

animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'};
Fp.animals = animals;
Fp.filtfunction = 'behavestate';
Fp = load_filter_params(Fp);

%% ---------------- Load BehaveStruct---------------------------------------------------
if loadBehaveStruct
    F = load_filter_output(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filtfunction', Fp.filtfunction);
    %     load(sprintf('%s%s%s.mat',FFanimdir, animalID, behavestruct));
end
%% ---------------- calculate StateSpace Model---------------------------------------------------
if calculateStateSpace
    %     statespace = struct;
    for ani = 1:length(Fp.animals)
        animal = Fp.animals{ani};
        andef = animaldef(animal);
        F(ani).animal = andef; % add it to F
        
        task(ani).animal = animal;
        task(ani).info = loaddatastruct(andef{2}, animal, 'task');
%         dayeps = evaluatefilter(task(ani).info, 'isequal($environment, ''wtrack'')');
        dayeps = evaluatefilter(task(ani).info, 'strfind($environment, ''track'')');
        envs = cellfetch(task(ani).info, 'environment');
        notwtrack = evaluatefilter(task(ani).info, '(~isequal($environment, ''openfield'') && ~isequal($environment, ''wtrack'') && isequal($type, ''run''))');
        if filterpureWdays
            % save list of ONLY pure wtrack days and epochs (not rotated, 6arm)
            dayeps = dayeps(~ismember(dayeps(:,1), notwtrack(:,1)),:);
        else
            dayeps = dayeps;
        end
        StateChanges = []; dayepinds = []; eplengths = [];
        StateChangesCell = {};
        for iep = 1:length(dayeps(:,1))
            try
                tmp = F(ani).statechanges{dayeps(iep,1)}{dayeps(iep,2)}.statechangeseq;
                %get rid of the first state end change time as it doesn't have a an accurate start time
                StateChangesCell{iep,1} = [tmp(2:end,:) repmat(dayeps(iep,:), length(tmp(2:end,1)), 1)];
                %         eplengths = [eplengths; length(StateChangesCell{iep}(:,1))];
                %             dayepinds = [dayepinds; repmat(wdayeps(iep,:), length(StateChangesCell{iep}(:,1)), 1)];
            catch
%                 pause
                continue
            end
        end
        
        StateChanges = cell2mat(StateChangesCell); % [allepsMat; allepsStateChanges{iep}];
        daybounds = find(diff(StateChanges(:,12)));
        epbounds = find(abs(diff(StateChanges(:,13))));
        
        SCFields = F(ani).statechanges{dayeps(1,1)}{dayeps(1,2)}.fields;
        chance = 0.5;
        %       The first column of probcorrect is the mode
        %       The second column of probcorrect is the lower 5% confidence bound
        %       The third column of probcorrect is the upper 5% confidence bound
        %       DR04/30/17  Fourth col is certainty matrix
        %       first row gets exluded
        [pcALL] = getestprobcorrect(StateChanges(:,7),chance,0,0);
        if length(pcALL(2:end,1)) ~= length(StateChanges(:,1))
            error('length of state space results do not match trial input length')
        end
        
        [pcINB] = getestprobcorrect(StateChanges(StateChanges(:,8)==1,7),chance,0,0);
        [pcOUTB] = getestprobcorrect(StateChanges(StateChanges(:,9)==1,7),chance,0,0);
        F(ani).statespace.allbound = [pcALL(2:end,:) StateChanges(:,[12 13])];
        F(ani).statespace.inbound = [pcINB(2:end,:) StateChanges(StateChanges(:,8)==1,[12 13])];
        F(ani).statespace.outbound = [pcOUTB(2:end,:) StateChanges(StateChanges(:,9)==1,[12 13])];
        F(ani).statespace.pcFields = 'mode lower5 upper5 certainty day epoch';
        %         statespace(ani).behavestruct = behavestruct;
        %         statespace(ani).epochEnvironment = epochEnvironment;
        F(ani).statespace.allepsMat = StateChanges;
        F(ani).statespace.allepsMatFields = SCFields;
        F(ani).statespace.eplengths = epbounds;
        F(ani).statespace.daylengths = daybounds;
    end
end

%% ---------------- save StateSpaceResults into BehaveState struct---------------------------------------------------
if saveStateSpaceResults == 1
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave)
end
%% ---------------- plot statespace---------------------------------------------------
if plotStateSpace
    Pp = load_plotting_params('behaveperform');
    for ani= 1:length(Fp.animals)
        animal = Fp.animals{ani};
        statespace = F(ani).statespace;
        SPmat = F(ani).statespace.allepsMat;
        for iseq = 1:3 % 1 Allbound; 2 Inbound; 3 Outbound
            switch iseq
                case 1
                    pc = statespace.allbound;
                    behavperform = SPmat(:,7);
                    trialType = 'Allbound';
                case 2
                    pc = statespace.inbound;
                    trialType = 'Inbound';
                    behavperform = SPmat(SPmat(:,8)==1,7);
                case 3
                    pc = statespace.outbound;
                    trialType = 'Outbound';
                    behavperform = SPmat(SPmat(:,9)==1,7);
            end
            modecol = 1;
            lower5col = 2;
            upper5col = 3;
            certaintycol = 4;
            
            % plot statespace results
            tALL=1:size(pc,1)-1;
            if savefigs && ~pausefigs
                close all
                ifig = figure('Visible','off','units','normalized','position', ...
                    Pp.position);
            else
                ifig = figure('units','normalized','position',Pp.position);
            end
            set(gcf,'color','white')
            subplot(2,1,1);
            hold on;
            daybounds = [0; find(diff(pc(:,5)))];
            epbounds = [0; find(abs(diff(pc(:,6))))];
            
            line([epbounds'; epbounds'], [0 .3], 'Color', [.5 .5 .5], 'LineStyle', '-', 'LineWidth', 1)
            line([daybounds'; daybounds'], [0 .4], 'Color', [0 0 0], 'LineStyle', '-', 'LineWidth', 2)
            days = unique(pc(:,5));
            for d = 1:length(days)
                text(daybounds(d), .5, sprintf('%d', days(d)), 'Color',[0 0 0], 'FontSize', 8)
            end
            
            plot(tALL, pc(2:end,modecol),'b-', 'LineWidth', 2); %plot behavior SS score
            errfillAll = fill([tALL fliplr(tALL)],[pc(2:end,lower5col); flipud(pc(2:end,upper5col))],[0 0 1],'linestyle','none');
            set(errfillAll, 'FaceAlpha', .1)
            hold on; [x, y] = find(behavperform > 0);
            h = plot(x,y-0.03,'+'); set(h, 'MarkerFaceColor','none');
            set(h, 'MarkerEdgeColor', [.2 .7 .2]);
            hold on; [x, y] = find(behavperform == 0);
            h = plot(x,y-0.05,'+'); set(h, 'MarkerFaceColor', 'none');
            set(h, 'MarkerEdgeColor', [.5 .5 .5]);
            text(0, 1.1, 'rewarded', 'Color',[.2 .7 .2], 'FontSize', 10)
            axis([1 tALL(end)  0 1.05]);
            line([1 tALL(end)], [chance  chance], 'Color', [.5 .5 .5], 'LineStyle', '--');
            xlabel('Trial # (Day:black, Epoch:gray)')
            ylabel([{'Probability of a'};{'Correct Response'}])
            
            subplot(2,1,2)
            hold on;
%             text(0, 1.1, 'Day', 'Color',[0 0 0], 'FontSize', 10)
%             text(daybounds(2), 1.1, 'Epoch', 'Color',[.5 .5 .5], 'FontSize', 10)
%             line([epbounds'; epbounds'], [0 1], 'Color', [.5 .5 .5], 'LineStyle', '-', 'LineWidth', 1)
%             line([daybounds'; daybounds'], [0 1], 'Color', [0 0 0], 'LineStyle', '-', 'LineWidth', 2)
            line([epbounds'; epbounds'], [0 .3], 'Color', [.5 .5 .5], 'LineStyle', '-', 'LineWidth', 1)
            line([daybounds'; daybounds'], [0 .4], 'Color', [0 0 0], 'LineStyle', '-', 'LineWidth', 2)
            for d = 1:length(days)
                text(daybounds(d), .5, sprintf('%d', days(d)), 'Color',[0 0 0], 'FontSize', 8)
            end
            %         line([cumsum(statespace.eplengths); cumsum(statespace.eplengths)], [0 1], 'Color', [.9 .9 .9], 'LineStyle', '-', 'LineWidth', 1)
            plot(tALL,1 - pc(2:end,certaintycol),'k', 'LineWidth', 2)
            line([ 1 tALL(end)],[0.90 0.90], 'LineStyle', '-');
            line([ 1 tALL(end)],[0.99 0.99], 'LineStyle', '--');
            line([ 1 tALL(end)],[0.95 0.95], 'LineStyle', '-.');
            axis([1 tALL(end)  0 1]);
            grid off;
            xlabel('Trial #')
            ylabel('Certainty')
            
            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('%s statespace %s performance d%d-%d%s', animal, trialType, days(1), days(end));
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', 12);
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                save_figure(Fp.paths.figdirectory, Fp.paths.filenamesave, sprtit)
                close all
            end
            close all;
        end
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
