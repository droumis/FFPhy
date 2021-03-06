


run_ff = 0;
savedata = run_ff;

load_stack = 0;
get_ripstate = 0;
calc_AS = 0;

batch_load_AS = 0; % use parfor 30 and 1 animal
calc_PWR = 1; % use parfor 4
load_PWR = 0;
run_permtest = 0;
load_ntinfo = 0;

plotfigs = 0;
savefigs = 1;
pausefigs = 0;

dsamp = 2;

% use_filters = {'firstwell', 'noise'};
Fp.animals = {'D12'};
Fp.add_params = {'wtrack', 'wavelets4-300Hz', 'excludeNoise','excludePriorFirstWell'};
Fp.filtfunction = 'dfa_riptriglfp';
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
me = animaldef('Demetris');

%% load lfpstack
if load_stack
    tic
    fprintf('loading %s\n', Fp.animals{:})
    lfpstack = load_data(Fp.paths.resultsDirectory, 'riptriglfpstack_wtrack', Fp.animals);
    lfpstacktimer = toc
end
%% Calculate and Save Analytic Signal %% this should use 30 workers (on typhoon)
if calc_AS
    tic
    computeAnalyticSignal(lfpstack, 'waveSet', Fp.waveSet, 'overwrite', 1);
    calcAStimer = toc
end
%% get behavioral state for each rip
if get_ripstate
    ripstate = getStateFilters(lfpstack);
end
%%
if batch_load_AS
    tic
    %     for ian = 1:length(Fp.animals)
    ian = 1;
    lfptypes = {'eeg'};
    wp = getWaveParams(Fp.waveSet, []);
    animal = Fp.animals{ian};
    fprintf('animal %s \n', animal);
    andef = animaldef(animal);
    % batch load
    anidx = find(strcmp(animal, {lfpstack.animal}));
    ntrodes = lfpstack(anidx).ntrodes;
    as = cell(1,length(lfptypes));
    for itype = 1:length(lfptypes)
        asnt = cell(1,length(ntrodes));
        parfor nti = 1:length(ntrodes) % should use 30 parfor
            nt = ntrodes(nti);
            asnt{nti} = loadAS(Fp.animals{anidx}, nt, Fp.waveSet, lfptypes{itype});
        end
        as{itype} = asnt;
        clear asnt
        % vectorize all ntrodes
        a = [as{itype}{:}];
        b = {a.as};
        asdatacat{itype} = cat(4,b{:});
    end
    batchloadtimer = toc
    %     end
end
%% Compute and Save Power %% this should use 4 workers (on typhoon)
if calc_PWR
    tic
    getPower(asdatacat, lfpstack, ripstate, Fp, 'lfptypes', {'eeg'}, ...
        'ripstatetypes', {'all'}, 'run_permutation_test', run_permtest, ...
        'savepower', 1);
    calcpwrtimer = toc
end

%% load D10 power
if load_PWR
    pwr = {};
    for ian = 1:length(Fp.animals)
        savestr = sprintf('%s/power/power_%s_%s_waveSet-%s_%s.mat', me{2}, ...
            Fp.animals{ian}, 'wtrack', Fp.waveSet, 'eeg');
        tmp = load(savestr);
        pwr = [tmp.pwrout; pwr];
    end
    clear tmp
end
%% load ntrode info
if load_ntinfo
    ntinfo = struct;
    for ian = 1:length(Fp.animals)
        ntinfo(ian).animal = Fp.animals{ian};
        andef = animaldef(ntinfo(ian).animal);
        ntinfo(ian).ntinfo = loaddatastruct(andef{2}, ntinfo(ian).animal, 'tetinfo');
    end
end
%% plot ftheatmap+zmap - per ntrode - all wtrack rips
% confirm that things look ok and there might be mec area groups
if plotfigs
    Pp = load_plotting_params({'defaults', 'power'});
    wp = getWaveParams('4-300HzJustfreqs', []);
    for ian = 1:numel(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        anidx = find(strcmp({pwr.animal}, animal));
        invalidtets = evaluatefilter(ntinfo(anidx).ntinfo, 'strcmp($valid, ''no'')');
        invalidtets = unique(invalidtets(:,3));
        areas = {'ca1', 'mec'};
        for ar = 1:length(areas) % for each area
            
            areantrodes = evaluatefilter(ntinfo(anidx).ntinfo, ...
                sprintf('isequal($area, ''%s'')', areas{ar}));
            areantrodes = unique(areantrodes(~ismember(areantrodes(:,3), invalidtets),3));
            %             conditions = pwr(ian).ripstate(1).statesetsfields;
            % for each condition
            for co = 1%:length(conditions)
                %                 nlfptypes = pwr(ian).lfptype;
                %                 for ty = 1:nlfptypes % for each eeg type
                if savefigs && ~pausefigs
                    close all
                    ifig =figure('Visible','off','units','normalized','position', ...
                        Pp.position);
                else
                    ifig = figure('units','normalized','position',Pp.position);
                end
                set(gcf,'color','white')
                for nti = 1:length(areantrodes)
                    sf = subaxis(4,5,nti, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                        'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, ...
                        'MarginBottom', Pp.MgBm);
                    nt = areantrodes(nti);
                    
                    idata2plot = squeeze(pwr(anidx).meandbpower{co}(:,:,:,nt))';
                    idata2plot = trim2win(idata2plot, Fp.srate, Pp.pwin, 'dsamp', dsamp);
                    time = linspace(-Pp.pwin(1), Pp.pwin(2), length(idata2plot(1,:)));
                    contourf(sf, time, wp.frex, idata2plot, Pp.contourRes, ...
                        'linecolor','none');
                    set(gca,'ydir','normal','yscale','log');
                    
                    if strcmp(areas{ar}, 'mec')
                        coloraxis = [-1 1]; %[-5 50];
                    elseif strcmp(areas{ar}, 'ca1')
                        coloraxis = [-2 5]; %[-5 50];
                    end
                    caxis(coloraxis)
                    colormap(Pp.usecolormap)
                    
                    % now overlap p<.05 thresholded zmap as contour
%                     MCmax = squeeze(pwr(anidx).permutes_max{co}(:,:,:,nt))';
%                     MCin = squeeze(pwr(anidx).permutes_min{co}(:,:,:,nt))';

%                     eval(sprintf('irawdata = squeeze(ixpc.%sout{ian}{nt.numsumSortInds(int)}(:,iDT,:))'';', calcfunction));
                    
                    % plot MC minmax thresh contour
%                     irawdata = trim2win(irawdata, srate, plotwin);
%                     hold on;
%                     MCthresh = MCminmax(ceil(length(MCminmax)*(1-pval)));
%                     MCthreshmap = idata2plot;
%                     MCthreshmap(abs(idata2plot)<MCthresh) = 0;
%                     [~,mc] = contour(intfig,plottimeWin,frex,logical(MCthreshmap),1);
%                     mc.LineColor = fig.mcLineColor;
                    
                    hold on
                    % thresholded single pix zmask
                    zmask2plot = squeeze(pwr(anidx).threshmean{co}(:,:,nt))'; % ntrode is in 3rd dim here
                    zmask2plot = trim2win(zmask2plot, Fp.srate, Pp.pwin, 'dsamp', dsamp);
                    [~,h] = contour(sf, time, wp.frex, logical(zmask2plot), 1);
                    h.LineColor = 'black';
                    
                    %
                    hold on;
                    title(sprintf('%s nt%d',areas{ar},nt), 'FontSize',Pp.FontS,'FontWeight',Pp.FontW, 'FontName', ...
                        Pp.FontNm)
                    if mod(nti+4,5) == 0
                        xlabel('time s', 'FontSize',Pp.FontS,'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
                        ylabel('freq Hz','FontSize',Pp.FontS, 'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
                        ytickskip = 2:4:wp.numfrex;
                        set(gca,'ytick', round(wp.frex(ytickskip)),'yticklabel',round(wp.frex(ytickskip)))
                    else
                        set(gca, 'xtick', [])
                        set(gca, 'ytick', [])
                    end
                    yl = ylim;
                    line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
                end
                %% super
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('%s %s mean dbPower p01zmap %s %s %s %s', animal, areas{ar}, ...
                    Fp.add_params{1}, Fp.add_params{2});
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center','FontSize', 12);
                
                %% ---- pause, save figs ----
                if pausefigs
                    pause
                end
                if savefigs
                    me = animaldef('Demetris');
                    save_figure(sprintf('%s/wavepower/',me{4}), 'wavepower', sprtit)
                    close all
                end
                close all;
                %                 end
            end
        end
    end
end


%% run the other animals to get power. all wtrack rips
% confirm that things look ok and group into mec area

%% plot - per animal - per ntrode - all wtrack rips

%% plot - all animals - per area - all wtrack rips

%% Do all that for phase clustering (this is where i need to be by monday)

%% Conditions
% - circular continuous: head direction
% - linear continuous : estimated % correct, learning rate, velocity
% - binary : last well rewarded, next well rewarded (correct/error), next well triggered 
%        is center well (inbound/outbound), next well is outer left, next
%        well is outer right.
%% plot ftmaps+zmap for power and phase - per condition - per area - per animal

%% plot ftmaps+zmap for power and phase - per condition - per area - all animals

