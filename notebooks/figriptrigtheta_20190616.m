
%{
eventually
- split the pre-wsr and the post-swr.. maybe even ignore the first 100-200ms of the swr..
but first just get the mean or median theta power per ripple.
- so load the referenced lfp,
- create the analytic signal
%}


run_ff = 0;
savedata = run_ff;

load_lfp = 1;
stack_lfp = 1;
load_events = 1;
make_filter_vecs = 1;

plotfigs = 0;
% plot_heatrast_traces_perep_allnt = 0;
plot_traces = 0;
plot_heatrasters = 1;

savefigs = 1;
pausefigs = 0;

use_filters = {'firstwell', 'noise'};
animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ4'}; %, 'D13'}; % next run ff jz1, jz3, jz4
add_params = {'wtrack'};

Fp.animals = animals;
Fp.add_params = add_params;
Fp.filtfunction = 'dfa_riptriglfp';
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
filetail = '';

%% run filter/func
if run_ff == 1
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes', ...
        Fp.tetfilter, 'excludetime', Fp.timefilter, 'iterator', Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
    F = runfilter(F);
    for d = 1:length(F)
        F(d).datafilter_params = Fp;
    end
end
%% save data
if savedata == 1
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, 'filetail',...
        sprintf('_%s%s', Fp.epochEnvironment, filetail))
end
%% load lfp
if load_lfp
    F = load_filter_output(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', sprintf('_%s%s', Fp.epochEnvironment, filetail));
end
%% stack LFP
if stack_lfp
    lfpstack = stack_riptriglfp(F);
end
%% load events, infostructs, timefilters
if load_events
    for ian = 1:length(Fp.animals)
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo{ian} = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        rips{ian} = loaddatastruct(aninfo{2}, animal, 'ca1rippleskons');
        noise{ian} = loaddatastruct(aninfo{2}, animal, 'ca1noisekons');
        refrip{ian} = loaddatastruct(aninfo{2}, animal, 'refrippleskons');
        immobile{ian} = getextendedimmobility(aninfo{2},aninfo{3}, ...
            F(1).epochs{1});
        firstwell{ian} = getpriortofirstwell(aninfo{2},aninfo{3}, ...
            F(1).epochs{1});
    end
end
%% Add filter vecs to stack
if make_filter_vecs
    for ian = 1:length(Fp.animals)
        filtervecs = struct;
        % need to go event by event bc the different epochs
        for irip = 1:length(lfpstack(ian).ripStartTime)
            day = lfpstack(ian).day(irip);
            epoch = lfpstack(ian).epoch(irip);
            ripstarttime = lfpstack(ian).ripStartTime(irip);
            
            % std ADD std directly to the lfpstack
            ripind = knnsearch(rips{ian}{day}{epoch}{1}.starttime,ripstarttime);
            ripstd = rips{ian}{day}{epoch}{1}.maxthresh(ripind);
            lfpstack(ian).std(irip, 1) = rips{ian}{day}{epoch}{1}.maxthresh(ripind);
            
            % pre first well
            exclrip = isExcluded(lfpstack(ian).ripStartTime(irip), ...
                firstwell{ian}{day}{epoch}.prefirst_list);
            filtervecs.firstwell(irip, 1) = exclrip;
            
            % immobile
            exclrip = isExcluded(lfpstack(ian).ripStartTime(irip), ...
                immobile{ian}{day}{epoch}.immobile_list);
            filtervecs.immobile(irip, 1) = exclrip;
            
            % near large noise
            noisewin = 1; %s
            noisestd = 15; %std
            near = abs(noise{ian}{day}{epoch}{1}.starttime - ripstarttime) <= noisewin;
            filtervecs.noise(irip, 1) = any(noise{ian}{day}{epoch}{1}.maxthresh(near) >= ...
                noisestd);
            
            % near large refripband events
            refripwin = 1; %s
            refripstd = 15; %std
            near = abs(refrip{ian}{day}{epoch}{1}.starttime - ripstarttime) <= refripwin;
            filtervecs.refrip(irip, 1) = any(refrip{ian}{day}{epoch}{1}.maxthresh(near) >= ...
                refripstd);
            
            % or very close to refrip events larger than ca1rip
            refripsmwin = .2;
            near = abs(refrip{ian}{day}{epoch}{1}.starttime - ripstarttime) <= refripsmwin;
            tmp = any(refrip{ian}{day}{epoch}{1}.maxthresh(near) >= ...
                ripstd);
            filtervecs.refrip(irip, 1) = any([filtervecs.refrip(irip, 1), tmp]);
        end
        lfpstack(ian).filterfields = fieldnames(filtervecs);
        lfpstack(ian).filtervecs = struct2array(filtervecs);
    end
end

%%
if 0
    ian = 1;
    t = find(strcmp(lfpstack(ian).lfptypes, 'eeg'));
    [unqdayeps, ~, ude_idx] = unique([lfpstack(ian).day lfpstack(ian).epoch], 'rows');
    epmean = cell2mat(arrayfun(@(x) mean(lfpstack(ian).data{t}(ude_idx==x,:,:)), ...
        1:length(unqdayeps(:,1)), 'un', 0)');
    nti = 11;
    plot(squeeze(epmean(:,:,nti))')
end

%%

if plotfigs
    if plot_traces
        Pp = load_plotting_params({'defaults', 'riptriglfp_perLFPtype_allntrodes'});
        for ian = 1:numel(Fp.animals)
            animal = Fp.animals{ian};
            anidx = find(strcmp({lfpstack.animal}, animal));
            ntrodes = lfpstack(anidx).ntrodes;
            dayep = [lfpstack(anidx).day lfpstack(anidx).epoch];
            
            use_filts = find(any(cell2mat(cellfun(@(x) strcmp(x, ...
                lfpstack(anidx).filterfields), ...
                use_filters, 'un', 0)), 2));
            exclude_rips = any(lfpstack(anidx).filtervecs(:,use_filts),2);
            exdayep = dayep(~exclude_rips,:);
            exstack = lfpstack(ian).data{t}(~exclude_rips,:,:);
            exstack(isnan(exstack)) = 0;
            for t = 1:length(lfpstack(ian).lfptypes)
                %         t = find(strcmp(lfpstack(ian).lfptypes, 'eeg'));
                [unqdayeps, ~, ude_idx] = unique(exdayep, 'rows');
                epzmeans = cell2mat(arrayfun(@(x) mean(zscore(exstack(ude_idx==x,:,:),[],2)), ...
                    1:length(unqdayeps(:,1)), 'un', 0)');
                
                invalidtets = evaluatefilter(ntinfo{ian}, 'isequal($valid, ''no'')');
                invalidtets = unique(invalidtets(:,3));
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
                    
                    subaxis(2,ceil(max(ntrodes)/2), ntrode, 'SpacingVert', Pp.SpVt, ...
                        'SpacingHoriz', Pp.SpHz, 'MarginLeft', Pp.MgLt, 'MarginRight', ...
                        Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    
                    a = repmat((1:size(epzmeans,1))',1, size(epzmeans,2));
                    b = squeeze(epzmeans(:,:,nti));
                    p = plot(Fp.time, (a+.5*b)', 'LineWidth', .1);
                    yticks([1:size(epzmeans,1)])
                    yticklabels(cellfun(@(x) num2str(x), flipud(num2cell(unqdayeps,2)), 'un', 0)')
                    axis tight
                    
                    xlabel('time s', 'FontSize',Pp.FontS,'FontWeight',Pp.FontW,'FontName', ...
                        Pp.FontNm)
                    ylabel('day-epoch zscore mean','FontSize',Pp.FontS, ...
                        'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
                    ntarea = ntinfo{ian}{unqdayeps(1,1)}{unqdayeps(1,2)}{ntrode}.area;
                    title(sprintf('%d-%s', ntrode, ntarea), 'FontSize',Pp.FontS, 'FontWeight',Pp.FontW, ...
                        'FontName', Pp.FontNm)
                    if nti ~= 1 && nti ~= ceil(max(ntrodes)/2)+1
                        xlabel('')
                        ylabel('')
                        set(gca, 'ytick', []);
                        set(gca, 'xtick', []);
                    end
                end
                %% super
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('%s %s %s %s %s %s squish', animal, lfpstack(anidx).lfptypes{t}, ...
                    Fp.paths.filenamesave(5:end), Fp.epochEnvironment, filetail, datestr(now, 'yyyymmdd'));
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
            end
        end
    end
     if plot_heatrasters
        Pp = load_plotting_params({'defaults', 'riptriglfp_perLFPtype_allntrodes'});
        for ian = 1:numel(Fp.animals)
            animal = Fp.animals{ian};
            anidx = find(strcmp({lfpstack.animal}, animal));
            ntrodes = lfpstack(anidx).ntrodes;
            dayep = [lfpstack(anidx).day lfpstack(anidx).epoch];
            
            use_filts = find(any(cell2mat(cellfun(@(x) strcmp(x, lfpstack(anidx).filterfields), ...
                use_filters, 'un', 0)), 2));
            exclude_rips = any(lfpstack(anidx).filtervecs(:,use_filts),2);
            exdayep = dayep(~exclude_rips,:);
            exstack = lfpstack(ian).data{t}(~exclude_rips,:,:);
%             exstack(isnan(exstack)) = 0;
            for t = 1:length(lfpstack(ian).lfptypes)
                
                invalidtets = evaluatefilter(ntinfo{ian}, 'isequal($valid, ''no'')');
                invalidtets = unique(invalidtets(:,3));
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
                    
                    subaxis(2,ceil(max(ntrodes)/2), ntrode, 'SpacingVert', Pp.SpVt, ...
                        'SpacingHoriz', Pp.SpHz, 'MarginLeft', Pp.MgLt, 'MarginRight', ...
                        Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    
                    a = repmat((1:size(epzmeans,1))',1, size(epzmeans,2));
                    b = squeeze(epzmeans(:,:,nti));
                    p = plot(Fp.time, (a+.5*b)', 'LineWidth', .1);
                    yticks([1:size(epzmeans,1)])
                    yticklabels(cellfun(@(x) num2str(x), flipud(num2cell(unqdayeps,2)), 'un', 0)')
                    axis tight
                    
                    excld_stack = squeeze(double(exstack(:,:,nti)));
                    
                    mididx = ceil(size(excld_stack,2)/2); % right now assumes center is rip start
                    ptime = Fp.time(mididx-(Pp.pwin(1)*Fp.srate):mididx+(Pp.pwin(2)*Fp.srate));
                    m = nanmean(excld_stack,2);
                    s = nanstd(excld_stack, [], 2);
                    z = (excld_stack-m)./s;
                    imagesc(ptime, 1:size(z,1), z)
                    colormap(parula)
                    
                    line([-Pp.pwin(1) Pp.pwin(2)], [epbounds'; epbounds'], 'color',[.9 .9 .9])
                    line([-Pp.pwin(1) Pp.pwin(2)], [daybounds'; daybounds'], 'color', [0 0 0])
                    line([0 0], [1 size(z,1)], 'color', [0 0 0], 'linestyle', '--')
                    ntarea = ntinfo{ian}{unqdayeps(1,1)}{unqdayeps(1,2)}{ntrode}.area;
                    title(sprintf('%d-%s', ntrode, ntarea), 'FontSize',Pp.FontS, 'FontWeight',Pp.FontW, ...
                        'FontName', Pp.FontNm)
                    caxis([-1 1])
                    xlabel('time s', 'FontSize',Pp.FontS,'FontWeight',Pp.FontW,'FontName', ...
                        Pp.FontNm)
                    ylabel('ripnum (day-b epoch-w)','FontSize',Pp.FontS, ...
                        'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
                    if nti ~= 1 && nti ~= ceil(max(ntrodes)/2)+1
                        xlabel('')
                        ylabel('')
                        set(gca, 'ytick', []);
                        set(gca, 'xtick', []);
                    end
                end
                %% super
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('%s %s %s %s %s %s squish heatrast', animal, lfpstack(anidx).lfptypes{t}, ...
                    Fp.paths.filenamesave(5:end), Fp.epochEnvironment, filetail, datestr(now, 'yyyymmdd'));
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
            end
        end
    end
end