


% load D12 power per rip..
% plot the 300 Hz traces across time/rips to detect where the noisy rips are
% just use the raw power.. don't db normalize to baseline
pconf = paramconfig;
run_ff = 0;
savedata = 0;
load_lfp = 0;
stack_lfp = 0;

load_stack = 0;
calc_AS = 0;
get_ripstate = 0;
load_ripstate = 0;
calc_PWR = 0;
calc_ITPC = 0;
run_permtest = 0;

load_raw_pwr = 0;
plot_raw_pwr = 1;

load_pwr = 0;
plot_pwr = 0;

load_itpc = 0;
plot_itpc = 0;

pausefigs = 1;
savefigs = 0;

Fp.animals = {'D12'};
Fp.filtfunction = 'dfa_riptriglfp';
Fp.add_params = {'wtrackdays', 'excludeNoise','excludePriorFirstWell', '<4cm/s', ...
    'wavelets4-300Hz'}; %'correcttrials',

Fp = load_filter_params(Fp, 'add_params', Fp.add_params);

Fp.uselfptype = 'eeggnd';
Fp.useripstates = {'onlywdays'}; %, 'rewarded', 'unrewarded', 'inbound' , 'outbound', ...
%     'rewarded_inbound', 'unrewarded_inbound', 'rewarded_outbound', 'unrewarded_outbound'};
plot_frex = [7 100 200 300];
%% run filter/func
if run_ff == 1
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes', ...
        Fp.tetfilter, 'excludetime', Fp.timefilter, 'iterator', Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
    F = runfilter(F);
    for d = 1:length(F); F(d).datafilter_params = Fp; end
end
%% save data
if savedata == 1
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, 'filetail',...
        sprintf('_%s%s', Fp.epochEnvironment, filetail))
end
%% load lfp
if load_lfp
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', sprintf('_%s%s', Fp.epochEnvironment, filetail));
end
%% stack LFP (ntrode x time x ripple)
if stack_lfp
    lfpstack = stack_riptriglfp(F);
    clear F
    save_data(lfpstack, Fp.paths.resultsDirectory, 'riptriglfpstack_wtrack')
end
%% load lfpstack
if load_stack
    lfpstack = load_data(Fp.paths.resultsDirectory, 'riptriglfpstack_wtrack', Fp.animals);
end

%% Analytic Signal
if calc_AS
    computeAnalyticSignal(lfpstack, 'waveSet', Fp.waveSet, 'overwrite', 1, ...
        'saveAnalyticSignal', 1, 'uselfptype', Fp.uselfptype);
end
%% RIPSTATE
if get_ripstate
    ripstate = getStateFilters(lfpstack);
    save_data(ripstate, [pconf.andef{2}, 'ripstate/'], 'ripstate_wtrack');
end
if load_ripstate
    ripstate = load_data([pconf.andef{2}, 'ripstate/'], 'ripstate_wtrack', Fp.animals);
end
%% Compute and Save MEAN Power, ITPC
if calc_PWR
    pwr = getPower(ripstate, Fp, 'uselfptype', Fp.uselfptype, 'ripstatetypes', ...
        Fp.useripstates, 'run_permutation_test', run_permtest, 'savepower', 1);
end
if calc_ITPC
    itpc = getITPC(ripstate, Fp, 'uselfptype', Fp.uselfptype, 'ripstatetypes', ...
        Fp.useripstates, 'run_permutation_test', run_permtest, 'saveresult', 1);
end
%% load power
if load_pwr
    savedir = sprintf('%s/power/', pconf.andef{2});
    savestr=sprintf('/power_waveSet-%s_%s_%s',Fp.waveSet,Fp.uselfptype,Fp.epochEnvironment);
    pwr = load_data(savedir, savestr, Fp.animals);
end
%% load itpc
if load_itpc
    savedir = sprintf('%s/itpc/', pconf.andef{2});
    savestr=sprintf('/itpc_waveSet-%s_%s_%s',Fp.waveSet,Fp.uselfptype,Fp.epochEnvironment);
    itpc = load_data(savedir, savestr, Fp.animals);
end

%% load raw power
if load_raw_pwr
    pwr = struct;
    for ian = 1:length(Fp.animals)
        animal = Fp.animals{ian};
        wp = getWaveParams(Fp.waveSet);
        pwr(ian) = load_data(sprintf('%s/analyticSignal/', me{2}), ...
            sprintf('AS_waveSet-%s_%s_power', wp.waveSet, Fp.uselfptype), animal);
    end
end
%%
if plot_raw_pwr
    Pp = load_plotting_params({'defaults', 'riptriglfp_perstatefrex_allntrodes'});
    for ian = 1:numel(Fp.animals)
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        anidx = find(strcmp({pwr.animal}, animal));
        ntrodes = lfpstack(anidx).ntrodes;
        dayep = [lfpstack(anidx).day lfpstack(anidx).epoch];
        
        % exclude invalid tets
        invalidtets = evaluatefilter(ntinfo, 'isequal($valid, ''no'')');
        invalidtets = unique(invalidtets(:,3));
        
        for fx = 1:length(plot_frex)
            frexidx = knnsearch(pwr.wp.frex', plot_frex(fx));
            knnfrex = round(pwr.wp.frex(frexidx));
            % plot per ripstate
            for rs = 1:length(Fp.useripstates)
                %         use_filts = find(any(cell2mat(cellfun(@(x) strcmp(x, lfpstack(anidx).filterfields), ...
                %             use_filters, 'un', 0)), 2));
                istate = Fp.useripstates{rs};
                istateidx = find(strcmp(istate, ripstate.statesetsfields));
                include_rips = ripstate.statesets(:,istateidx);
                %% ---- init fig----
                if savefigs && ~pausefigs
                    close all
                    ifig = figure('Visible','off','units','normalized','position', ...
                        Pp.position);
                else
                    ifig = figure('units','normalized','position',Pp.position);
                end
                set(gcf,'color','white')
                
                for nti = 1:length(ntrodes)
                    ntrode = ntrodes(nti);
                    if ismember(ntrode, invalidtets)
                        continue
                    end
                    
                    sf = subaxis(2,ceil(max(ntrodes)/2), nti, 'SpacingVert', Pp.SpVt, ...
                        'SpacingHoriz', Pp.SpHz, 'MarginLeft', Pp.MgLt, 'MarginRight', ...
                        Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);

                    exdayep = dayep(find(include_rips),:);
                    de = unique(exdayep, 'rows');
                    daybounds = find(diff(exdayep(:,1)));
                    epbounds = find(abs(diff(exdayep(:,2))));
                    
                    excld_stack = squeeze(double(pwr(anidx).pwr(nti, :, find(include_rips),frexidx)))';
                    idata2plot = trim2win(excld_stack, Fp.srate, Pp.pwin, ...
                        'dsamp', pwr(anidx).wp.dsamp);
                    %                         mididx = ceil(size(excld_stack,2)/2); % right now assumes center is rip start
                    ptime = linspace(-Pp.pwin(1),Pp.pwin(2),size(idata2plot,2));
                    m = nanmean(idata2plot,2);
                    s = nanstd(idata2plot, [], 2);
                    z = (idata2plot-m)./s;
                    imagesc(ptime, 1:size(z,1), z)
                    colormap(parula)
                    caxis(sf, 'auto')
                    
                    line([-Pp.pwin(1) Pp.pwin(2)], [epbounds'; epbounds'], 'color',[.9 .9 .9])
                    line([-Pp.pwin(1) Pp.pwin(2)], [daybounds'; daybounds'], 'color', [0 0 0])
                    line([0 0], [1 size(z,1)], 'color', [0 0 0], 'linestyle', '--')
                    title(sprintf('%d', ntrode), 'FontSize',Pp.FontS, 'FontWeight',Pp.FontW, ...
                        'FontName', Pp.FontNm)
                    %                         caxis([-1 1])
                    xlabel('time s', 'FontSize',Pp.FontS,'FontWeight',Pp.FontW,'FontName', ...
                        Pp.FontNm)
                    ylabel('ripnum (day-b epoch-w)','FontSize',Pp.FontS, ...
                        'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
                    if nti ~= 1
                        xlabel('')
                        ylabel('')
                        set(gca, 'ytick', []);
                        set(gca, 'xtick', []);
                    end
                end
                %% super
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('%s %s %s %dHz %s power', animal, istate, Fp.uselfptype, ...
                    knnfrex, Fp.epochEnvironment);
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center','FontSize', 12);
                spylabel = text(.02, .5, sprintf('ripnum'), 'Parent', sprtitleax, 'Units', ...
                    'normalized', 'Rotation', 90);
                set(spylabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center', 'FontSize', 12);
                
                spxlabel = text(.5, .02, sprintf('time'), 'Parent', sprtitleax, 'Units', ...
                    'normalized');
                set(spxlabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center', 'FontSize', 12);
                %% ---- pause, save figs ----
                if pausefigs
                    pause
                end
                if savefigs
                    save_figure(Fp.paths.figdirectory, Fp.paths.filenamesave, sprtit)
                    close all
                end
            end
        end
    end
end


