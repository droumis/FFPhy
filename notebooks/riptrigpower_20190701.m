
% goal: what is the rip trig lfp spectrum in ca1, mec super, deep

% exclude the pre frist well and the noie ripples
% label the ntrode by area so i can see which are ref
% create AS for unref as well and plot for both the ref and unref
% plot alongside riptrig lfp traces
% plot for each day
% how to test statistically the single condition spects

% save lfp stacked.
% exclude firstwell, noise rips
% i should just exclude the firstwell and noise rips from the beginning..

run_ff = 0;
savedata = run_ff;

load_lfp = 0;
stack_lfp = 0;
load_stack = 1;
get_ripstate = 1;
calc_AS = 1;
calc_PWR = 1;
load_PWR = 0;
run_permtest = 0;

plotfigs = 0;
savefigs = 0;
pausefigs = 0;

% use_filters = {'firstwell', 'noise'};
Fp.animals = {'JZ2', 'JZ3', 'JZ4'};
Fp.add_params = {'wtrack', 'wavelets4-300Hz', 'excludeNoise','excludePriorFirstWell'};
Fp.filtfunction = 'dfa_riptriglfp';
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
me = animaldef('Demetris');

filetail = '';
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
%% stack LFP into a cube per animal (ntrode x time x ripple)
if stack_lfp
    lfpstack = stack_riptriglfp(F);
    clear F
    save_data(lfpstack, Fp.paths.resultsDirectory, 'riptriglfpstack_wtrack')
end
if load_stack
    lfpstack = load_data(Fp.paths.resultsDirectory, 'riptriglfpstack_wtrack', Fp.animals);
end
%% Calculate and Save Analytic Signal
if calc_AS
    computeAnalyticSignal(lfpstack, 'waveSet', Fp.waveSet, 'overwrite', 1);
end

%% get behavioral state for each rip
if get_ripstate
    ripstate = getStateFilters(lfpstack);
end
%%
% %% add waveletFFT to the lfpstack
% if add_cmwFFT
%     for ian = 1:length(lfpstack)
%         eegidx = find(cellfun(@(x) strcmp(x,'eeg'),lfpstack(ian).lfptypes, 'un', 1));
%         wp = getWaveParams(Fp.waveSet, lfpstack(ian).data{eegidx});
%         lfpstack(ian).waveletFFT = createCMWfft(wp);
%     end
%     save_data(lfpstack, Fp.paths.resultsDirectory, 'riptriglfpstack_wtrack')
% end

% use_filts = find(any(cell2mat(cellfun(@(x) strcmp(x, excludefields), use_filters, 'un', 0)), 2));
%     [excludefields, excludestate] = getExclusionFilters(lfpstack);
% exclude_rips = any(excludestate,2);
%% Compute and Save Power from analytic signal for each ntrode, condition
if calc_PWR
    lfptypes = {'eeggnd', 'eeg'};
%     ntrodes = 1:30;
%     as = cell(1,length(lfptypes));
%     for itype = 1:length(lfptypes)
%         asnt = cell(1,length(ntrodes));
%         parfor nti = 1:length(ntrodes)
%             nt = ntrodes(nti);
%             asnt{nti} = loadAS(Fp.animals{1}, nt, Fp.waveSet, lfptypes{itype});
%         end
%         as{itype} = asnt;
%     end
    getPower(as,ripstate, Fp, 'lfptypes', lfptypes);
end

% %% Compute and Save ITPC from analytic signal for each ntrode, condition
% if calc_PWR
%     wp = getWaveParams(Fp.waveSet, lfpstack(ian).data{2});
%     pwr = getPower(ripstate, Fp, 'baseline', wp.baseind);
% end

%% load power
if load_PWR
    lfptypes = {'eeggnd', 'eeg'};
    pwr = {};
    for ian = 1:length(Fp.animals)
        for itype = 1:length(lfptypes)
            lfptype = lfptypes{itype};
            savestr = sprintf('%s/power/power_%s_%s_waveSet-%s_%s.mat', me{2}, ...
                Fp.animals{ian}, 'wtrack', Fp.waveSet, lfptype);
            tmp = load(savestr);
            pwr{end+1} = [tmp.pwrout; pwr];
        end
    end
end

%% runPermutationTest
% if run_permtest
% %     aIdx = find(ripstate(ian).statesets(:,2)); % rewarded (2)
% %     bIdx = find(ripstate(ian).statesets(:,3)); % unrewarded (3)
%     pwr = powerpermtest(pwr);
% end

%% load other stuff 
ntinfo = struct;
for ian = 1:length(Fp.animals)
    ntinfo(ian).animal = Fp.animals{ian};
    andef = animaldef(ntinfo(ian).animal);
    ntinfo(ian).ntinfo = loaddatastruct(andef{2}, ntinfo(ian).animal, 'tetinfo');
end
%% plot

% invalidtets = evaluatefilter(ntinfo{ian}, 'isequal($valid, ''no'')');
%         invalidtets = unique(invalidtets(:,3));

if plotfigs
    Pp = load_plotting_params({'defaults', 'power'});
    for ian = 1:numel(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        anidx = find(strcmp({pwr.animal}, animal));
        invalidtets = evaluatefilter(ntinfo(anidx).ntinfo, 'isequal($valid, ''no'')');
        invalidtets = unique(invalidtets(:,3));
        areas = {'ca1', 'mec'};
        for ar = 1:length(areas) % for each area
            areantrodes = evaluatefilter(ntinfo(anidx).ntinfo, ...
                sprintf('isequal($area, ''%s'')'), areas{ar});
            
            conditions = pwr(ian).ripstate(1).statesetsfields;
            % for each condition
            for co = 1:length(conditions)
                
                nlfptypes = pwr(ian).ripstate;
                for ty = 1:nlfptypes % for each eeg type
                    if savefigs && ~pausefigs
                        close all
                        ifig =figure('Visible','off','units','normalized','position', ...
                            Pp.position);
                    else
                        ifig = figure('units','normalized','position',Pp.position);
                    end
                    set(gcf,'color','white')
                    for nti = 1:length(areantrodes)
                        sf = subaxis(3,5,nti, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                            'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, ...
                            'MarginBottom', Pp.MgBm);
                        nt = areantrodes(nti);
                        idata2plot = squeeze(pwr(anidx).mediandbpower{nt}{co})';
                        idata2plot = trim2win(idata2plot, Fp.srate, Pp.pwin);
                        time = linspace(-Pp.pwin(1), Pp.pwin(2), length(idata2plot(1,:)));
                        contourf(sf, time, wp.frex, idata2plot, Pp.contourRes, ...
                            'linecolor','none');
                        set(gca,'ydir','normal','yscale','log');
                        if areas{ar} == 'mec'
                            coloraxis = [-1 1]; %[-5 50];
                        elseif areas{ar} == 'ca1'
                            coloraxis = [-2 5]; %[-5 50];
                        end
                        caxis(coloraxis)
                        colormap(Pp.usecolormap)
                        
                        hold on;
                        title(sprintf('nt%d',nt), 'FontSize',Pp.FontS,'FontWeight',Pp.FontW, 'FontName', ...
                            Pp.FontNm)
                        if int == 1
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
                    sprtit = sprintf('%s dbPower %s %s %s %s', animal, ...
                        pwr(anidx).ripstate(ian).statesetsfields{iset}, add_params{1}, ...
                        add_params{2}, areas{area});
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
                end
            end
        end
    end
end


