
% finish cleaning the rips.. then
% make day/epoch averaged LFP for selected ntrodes and selected frequencies
% probably tettags to mark for exclusion

% conditions to impose:
% DONE **most important!! exclude events that occur before the first well visit
% DONE *next most important: events that occur hasn't moved more than 10 cm in
% the last 60 seconds
% * create refkons (rip-like events detected on ref chan) and use it with noise kons to exclude
% for each animal, load the reference nt ripple lfp. power should be
% there.. run a filter on stf > 3.
% * what is the std of noise kons? try lowering it.. D12 is noisy
% * plot ref lfp, the rest of ca1 nt lfp, ref kons
% rerun the ripple figure printing using the additional exclusion criteria.
% run the waterfall plots and the heatraster plots using the additional
% exclusion criteria.
if 0
    savefigs = 1;
    pausefigs = 0;
    Fp = load_filter_params({'ripcleaning'});
    Fp.animals = {'D10'};
    for ian = 1:length(Fp.animals)
        andef = animaldef(Fp.animals{ian});
        linpos{ian} = loaddatastruct(andef{2},andef{3},'linpos');
    end
    %% exclude long immobility periods
    Pp = load_plotting_params('ripcleaning');
    for ian = 1:length(Fp.animals)
        animal = Fp.animals{ian};
        for day = 1:length(linpos{ian})
            if isempty(linpos{ian}{day})
                continue
            else
                for epoch = 1:length(linpos{ian}{day})
                    if isempty(linpos{ian}{day}{epoch})
                        continue
                    else
                        %% ---- init fig----
                        if savefigs && ~pausefigs
                            close all
                            ifig = figure('Visible','off','units','normalized','position', ...
                                Pp.position);
                        else
                            ifig = figure('units','normalized','position',Pp.position);
                        end
                        set(gcf,'color','white')
                        fprintf('%s %d %d\n',animal, day, epoch);
                        pre_excl_win = 30; %seconds
                        post_excl_win = 5; %seconds
                        ax_srate = 30; %Hz
                        cm_thresh = 5; %cm
                        sig = 600;
                        % linpos{ian}{day}{epoch}
                        
                        lindist = linpos{ian}{day}{epoch}.statematrix.linearDistanceToWells(:,1);
                        linddist_win = conv(diff(lindist),ones(pre_excl_win*ax_srate,1),'same');
                        
                        linddist_win = zeros(length(lindist),1);
                        linddist_win(2:end) = movsum(diff(lindist),[pre_excl_win*ax_srate post_excl_win*ax_srate],'omitnan');
                        
                        g = gaussian(sig,1200);
                        u = smoothvect(abs(linddist_win), g);
                        exl_immobile =  u <= cm_thresh;
                        excl_immobile_list = vec2list(exl_immobile, linpos{ian}{day}{epoch}.statematrix.time);
                        
                        % smoothvect(exl_immobile)
                        % sum(exl_immobile)/length(exl_immobile)
                        
                        p = plot(exl_immobile*100, '+', 'DisplayName','exl_immobile');
                        hold on
                        p.Color = [1 0 1 .01];
                        p.MarkerSize = 10;
                        sm = plot(u, 'k', 'DisplayName','moving lindist'); % smoothed, abs,
                        ld = plot(lindist, 'b', 'DisplayName','lindist');
                        
                        
                        % exlude times before first well visit
                        wellfirstidx = find(diff(linpos{ian}{day}{epoch}.statematrix.wellExitEnter(:,2)),1);
                        centerwell_enteridx = find(lindist < 4, 1);
                        outerwell_enteridx = find(lindist > max(lindist)-4 , 1);
                        firstidx = min([outerwell_enteridx centerwell_enteridx wellfirstidx]);
                        
                        excl_prefirst = (1:length(linpos{ian}{day}{epoch}.statematrix.time) < firstidx)';
                        excl_prefirst_list = vec2list(excl_prefirst, linpos{ian}{day}{epoch}.statematrix.time);
                        
                        pf = plot(excl_prefirst*100, 'y+', 'DisplayName','exl_prefirst');
                        pf.MarkerSize = 8;
                        hold off
                        axis tight
                        xlabel('sample')
                        ylabel('lindist cm')
                        legend([p ld sm pf],{'immobile', 'lindist', 'rolling diff', 'prefirst'})
                        %% super
                        sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                        sprtit = sprintf('lindist exclude times %s %d %d', animal, day, epoch);
                        iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                        set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                            'horizontalAlignment', 'center','FontSize', 12);
                        %% ---- pause, save figs ----
                        if pausefigs
                            pause
                        end
                        if savefigs
                            save_figure(Fp.paths.figdirectory, 'riptrig_all', sprtit)
                            close all
                        end
                    end
                end
            end
        end
    end
    %%
    
    %% plot ref tet ripband trace
    plot(time, zscore(double(reflfp{day}{epoch}{reftet}.data(1:end,[1 3]))), 'LineWidth', 1);
    fields = strsplit(strrep(reflfp{day}{epoch}{reftet}.fields,'_', '-'), ' ');
    legend(fields([1 3]))
    hold on
    
    refstdmask = zscore(double(reflfp{day}{epoch}{reftet}.data(1:end,3))) > refstdthresh;
    time = reflfp{day}{epoch}{reftet}.starttime + (0:length(reflfp{day}{epoch}{reftet}.data(:,1))-1)/reflfp{day}{epoch}{reftet}.samprate;
    
    refriplist = vec2list(refstdmask, time);
    axis tight
    yl = ylim;
    xs = refriplist(:,1);
    xe = refriplist(:,2);
    patch([xs'; xe'; xe'; xs'], repmat([yl(1) yl(1) yl(2) yl(2)]', 1, length(xe)),'k', ...
        'FaceAlpha',.5, 'edgecolor','none');
    hold off
    %% load refripplekons, plot on ref ripband trace
    ian = 1;
    animal = Fp.animals{ian};
    andef = animaldef(animal);
    refrips = loaddatastruct(andef{2}, andef{3}, 'refrippleskons');
    ca1rips = loaddatastruct(andef{2}, andef{3}, 'ca1rippleskons');
    
    day = 1;
    epoch = 2;
    refstdthresh = 5;
    tetinfo = loaddatastruct(andef{2}, andef{3}, 'tetinfo');
    
    reftet = evaluatefilter(tetinfo{day}{epoch}, 'isequal($area, ''ref'')');
    
    refriplfp = loadeegstruct(andef{2}, andef{3}, 'ripple', day, epoch, reftet);
    zrefriplfpmag = zscore(double(refriplfp{day}{epoch}{reftet}.data(:,3)'), [], 2);
    
    reflfp = loadeegstruct(andef{2}, andef{3}, 'eeg', day, epoch, reftet);
    zreflfp = zscore(double(reflfp{day}{epoch}{reftet}.data'), [], 2);
    
    % CA1
    ca1tets = evaluatefilter(tetinfo{day}{epoch}, 'isequal($area, ''ca1'')');
    ca1RIPlfp = loadeegstruct(andef{2}, andef{3}, 'ripple', day, epoch, ca1tets);
    zca1RIPlfp = zscore(double(cell2mat(arrayfun(@(x) ca1RIPlfp{day}{epoch}{x}.data(:,3), ca1tets, 'un', 0)')'), [], 2);
    ca1EEGlfp = loadeegstruct(andef{2}, andef{3}, 'eeg', day, epoch, ca1tets);
    zca1EEGlfp = zscore(double(cell2mat(arrayfun(@(x) ca1EEGlfp{day}{epoch}{x}.data, ca1tets, 'un', 0)')'), [], 2);
    
    time = geteegtimes(ca1EEGlfp{day}{epoch}{ca1tets(1)});
    %
    % imagesc(time, [ca1tets; reftet], [zca1lfpstack(:,1:1500); zreflfp(1:1500)])
    %%
    clf
    ax1 = subplot(2,1,1);
    plot(time, [zreflfp' zrefriplfpmag']);
    hold on
    axis tight
    ylim([-5 5])
    yl = ylim;
    xs = refrips{day}{epoch}{1}.starttime;
    xe = refrips{day}{epoch}{1}.endtime;
    patch([xs'; xe'; xe'; xs'], repmat([yl(1) yl(1) yl(2) yl(2)]', 1, length(xe)),'k', ...
        'FaceAlpha',.5, 'edgecolor','none');
    
    ax2 = subplot(2,1,2);
    plot(time, [mean(zca1EEGlfp)' mean(zca1RIPlfp)']);
    hold on
    axis tight
    yl = ylim;
    xs = ca1rips{day}{epoch}{1}.starttime;
    xe = ca1rips{day}{epoch}{1}.endtime;
    patch([xs'; xe'; xe'; xs'], repmat([yl(1) yl(1) yl(2) yl(2)]', 1, length(xe)),'k', ...
        'FaceAlpha',.5, 'edgecolor','none');
    linkaxes([ax1,ax2],'x')
    
    %%
    %
    %
    %
    % win = 10;
    % r = 100;
    % xlim([xs(r)-win xe(r)+win])
    %
    % timefilter{end+1} = {'getconstimes', '($cons == 1)', ...
    %     'ca1rippleskons', 1,'consensus_numtets',consensus_numtets, ...
    %     'minstdthresh', minstdthresh,'exclusion_dur',exclusion_dur, ...
    %     'minvelocity', minvelocity,'maxvelocity',maxvelocity};
    %
    % [out] = getconstimes(animaldir,animalprefix, epochs, eventconsname, tetfilter, varargin);
    % % filterresults{i} = evaluatefilter2(dataStruct, timefilter{i}{2});
    %
    % excludetimes =
    % goodrips = ~isExcluded(ripstarttimes, excludetimes);
end
%%
if 1
    animals = {'JZ1'};
    add_params = {'wtrack'};
    
    run_lfp_ff = 0;
    save_lfp = run_lfp_ff;
    
    run_spikes_ff = 0;
    save_spikes = run_spikes_ff;
    
    loadstuff = 1;
    load_lfp = loadstuff;
    load_spikes = loadstuff;
    load_pos = loadstuff;
    load_rips = loadstuff;
    
    stackstuff = 1;
    stack_lfp = stackstuff;
    stack_spikes = stackstuff;
    
    plotfigs = 1;
    pausefigs = 0;
    savefigs = 1;
    
    %% run lfp filter/func
    Fp.animals = animals;
    Fp.add_params = add_params;
    Fp.filtfunction = 'dfa_riptriglfp';
    Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
    
    if run_lfp_ff == 1
        F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes', ...
            Fp.tetfilter, 'excludetime', Fp.timefilter, 'iterator', Fp.iterator);
        F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
        F = runfilter(F);
        for d = 1:length(F)
            F(d).datafilter_params = Fp;
        end
    end
    
    %% save lfp data
    if save_lfp == 1
        save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, 'filetail',...
            ['_', Fp.epochEnvironment])
    end
    %% load lfp
    if load_lfp
        F = load_filter_output(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
            Fp.animals, 'filetail', ['_', Fp.epochEnvironment]);
    end
    %% stack LFP
    if stack_lfp
        lfpstack = stack_riptriglfp(F);
    end
    %% run spiking FF
    sFp.animals = animals;
    sFp.add_params = add_params;
    sFp.filtfunction = 'dfa_riptrigspiking';
    sFp = load_filter_params(sFp, 'add_params',sFp.add_params);
    
    if run_spikes_ff == 1
        spikesF = createfilter('animal',sFp.animals, 'epochs', sFp.epochfilter, ...
            'cells',sFp.cellfilter, 'excludetime',sFp.timefilter, 'iterator',sFp.iterator);
        spikesF = setfilterfunction(spikesF,sFp.filtfunction,sFp.datatypes,sFp.options{:});
        spikesF = runfilter(spikesF);
        for d = 1:length(spikesF)
            spikesF(d).datafilter_params = sFp;
        end
    end
    %% save spikes data
    if save_spikes == 1
        save_data(spikesF,sFp.paths.filtOutputDirectory,sFp.paths.filenamesave, 'filetail',...
            ['_',sFp.epochEnvironment])
    end
    %% load spikes
    if load_spikes
        spikesF = load_filter_output(sFp.paths.filtOutputDirectory,sFp.paths.filenamesave, ...
            sFp.animals, 'filetail', ['_',sFp.epochEnvironment]);
    end
    %% stack spikes
    if stack_spikes
        spikestack = stack_riptrigspiking(spikesF);
    end
    %% load pos
    if load_pos
        for ian = 1:length(Fp.animals)
            andef = animaldef(Fp.animals{ian});
            pos{ian} = loaddatastruct(andef{2},andef{3},'pos');
            linpos{ian} = loaddatastruct(andef{2},andef{3},'linpos');
            behave{ian} = load(sprintf('%s/%sBehaveState.mat',andef{2}, andef{3}));
        end
    end
    %% load events, infostructs, timefilters
    if load_rips
        for ian = 1:length(Fp.animals)
            animal = Fp.animals{ian};
            aninfo = animaldef(animal);
            ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
            ntinfoAll = cellfetch(ntinfo, '', 'alltags', 1);
            rips{ian} = loaddatastruct(aninfo{2}, animal, 'ca1rippleskons');
            noise{ian} = loaddatastruct(aninfo{2}, animal, 'ca1noisekons');
            refrip{ian} = loaddatastruct(aninfo{2}, animal, 'refrippleskons');
            immobile{ian} = getextendedimmobility(aninfo{2},aninfo{3}, ...
                spikesF(ian).epochs{1});
            firstwell{ian} = getpriortofirstwell(aninfo{2},aninfo{3}, ...
                spikesF(ian).epochs{1});
        end
    end
    
    %% plot
    if plotfigs
        Pp = load_plotting_params({'riptrigall'});
        for ian = 1:length(Fp.animals)
            animal = Fp.animals{ian};
            hwin = floor(length(lfpstack(ian).data{1}(1,:,1))/2);
            time = [1/Fp.srate*-hwin:1/Fp.srate:0 ...
                1/Fp.srate:1/Fp.srate:1/Fp.srate*hwin];
            for t = 1:length(lfpstack(ian).lfptypes)+1
                gridspacemat(t,:,:) = Pp.tscale(t)* ...
                    permute(lfpstack(ian).ntrodes*ones(1, ...
                    length(time)), [3 2 1]);
            end
            %         iderip = 0;
            andef = animaldef(animal);
            tetinfo = loaddatastruct(andef{2}, andef{3}, 'tetinfo');
            reftet_keys = evaluatefilter(tetinfo, 'isequal($area, ''ref'')');
            ca1tet_keys = evaluatefilter(tetinfo, 'isequal($area, ''ca1'')');
            % make a lookup for day epoch irip irip within dayep
            rel_irip = cell2mat(cellfun(@(x) [1:x]', ...
                num2cell(lfpstack(ian).numrips_perep),'un', 0));
            iripmat = [[1:length(lfpstack(ian).epoch)]' lfpstack(ian).day ...
                lfpstack(ian).epoch rel_irip];
            
            for irip = 1:length(iripmat(:,1))
                day = lfpstack(ian).day(irip);
                epoch = lfpstack(ian).epoch(irip);
                ripstarttime = lfpstack(ian).ripStartTime(irip);
                reftet = reftet_keys(ismember(reftet_keys(:,[1 2]),[day epoch], 'rows'),3);
                ca1tets = ca1tet_keys(ismember(ca1tet_keys(:,[1 2]), [day epoch], 'rows'),3);
                
                                
                % exclude events prior to the first well visit
                if isExcluded(ripstarttime, firstwell{ian}{day}{epoch}.prefirst_list) 
                    fprintf('skipping bc within prefirst period \n')
                    continue
                end
                
                % exclude events within an extended immobility period
                if isExcluded(ripstarttime, immobile{ian}{day}{epoch}.immobile_list) 
                    fprintf('skipping bc within an extended immobility period\n')
                    continue
                end
                                
                % std thresh
                ripind = knnsearch(rips{ian}{day}{epoch}{1}.starttime,ripstarttime);
                ripstd = rips{ian}{day}{epoch}{1}.maxthresh(ripind);
                fprintf('rip%d start:%d\n', irip, ripstarttime)
                if ripstd < Pp.stdthresh
                    fprint('skipping bc outside of stdthresh \n')
                    continue
                end
                
                % exclude events within 1 second of a large noise event
                noiseripmaxstd = 15;
                noiseripnearwin = 1; %s
                noiseind = knnsearch(noise{ian}{day}{epoch}{1}.starttime,ripstarttime);
                nearnoiseripmaxstd = noise{ian}{day}{epoch}{1}.maxthresh(noiseind);
                if nearnoiseripmaxstd > noiseripmaxstd && ...
                        abs(noise{ian}{day}{epoch}{1}.starttime(noiseind) - ...
                        ripstarttime) < noiseripnearwin
                    fprintf('skipping bc near %d std noise event \n', nearnoiseripmaxstd)
                    continue
                end
                
                % exclude events within 1 second of a refrip event
                refripmaxstd = 15;
                refripnearwin = 1; %s
                refripind = knnsearch(refrip{ian}{day}{epoch}{1}.starttime,ripstarttime);
                nearrefripmaxstd = refrip{ian}{day}{epoch}{1}.maxthresh(refripind);
                if nearrefripmaxstd > refripmaxstd && ...
                        abs(refrip{ian}{day}{epoch}{1}.starttime(refripind) - ...
                        ripstarttime) < refripnearwin
                    fprintf('skipping bc near %d std refrip event \n', nearrefripmaxstd)
                    continue
                end
                
                % exclude events with refrip larger than ca1rip
                if nearrefripmaxstd > ripstd && ...
                        abs(refrip{ian}{day}{epoch}{1}.starttime(refripind) - ...
                        ripstarttime) < .5
                    fprintf('skipping bc near %d std refrip >  %d std ca1rip  \n', nearrefripmaxstd, ripstd)
                    continue
                end
                
                %% ---- init fig----
                if savefigs && ~pausefigs
                    close all
                    ifig = figure('Visible','off','units','normalized','position', ...
                        Pp.position);
                else
                    ifig = figure('units','normalized','position',Pp.position);
                end
                set(gcf,'color','white')
                %% PLOT SPIKES
                
                %             ha = tight_subplot(2, length(dstack(ian).lfptypes), [.05 .02],[.1 .1],[.05 .05]);
                subaxis(2,Pp.nrow,1, 'Spacing', Pp.Spacing, ...
                    'SpacingVert', Pp.SpacingVert, 'SpacingHoriz', Pp.SpacingHoriz, 'Padding', ...
                    Pp.Padding, 'MarginLeft', Pp.MarginLeft, 'MarginRight', Pp.MarginRight,...
                    'MarginTop', Pp.MarginTop, 'MarginBottom', Pp.MarginBottom);
                
                irelrip = iripmat(irip, end);
                irip_muspikes = squeeze(spikestack(ian).muspikes{day}{epoch}(irelrip,:,:))';
                fprintf('%s %d %d allrip%d derip%d \n', animal, day, epoch, irip, irelrip);
                iripmuntrodes = spikestack(ian).mu_de_keys{day}{epoch}(:,3);
                
                [sx, sy] = find(irip_muspikes');
                ca1y = find(ismember(iripmuntrodes, ca1tets));
                sz = zeros(length(sy), 3);
                sz(ismember(sy, ca1y), 3) = 1; % make ca1 blue
                refy = find(ismember(iripmuntrodes, reftet));
                sz(sy==refy, 2) = 1; % make ref green
                
                f1 = scatter(sx/1000-1.001, sy, 10, sz,'+'); % could ad 'z' data for coloring
                %             f1.MarkerFaceAlpha = 0.1;
                f1.MarkerEdgeAlpha = 0.05;
                %             f1.MarkerEdgeColor = [0 0 0];
                axis tight
                xlim([-1 1])
                yticks(iripmuntrodes);
                yticklabels({num2str(iripmuntrodes)});
                ylabel('ntrode','FontSize',8,'FontWeight','bold', 'FontName','Arial')
                title('mu', 'FontSize',8,'FontWeight','bold', 'FontName','Arial')
                
                yl = ylim;
                
                %% plot event patches in window
                winstarttime = lfpstack(ian).ripStartTime(irip)-(hwin*1/Fp.srate);
                winendtime = lfpstack(ian).ripStartTime(irip)+(hwin*1/Fp.srate);

                % plot rips in win
                ripindsinwin = find(rips{ian}{day}{epoch}{1}.starttime>winstarttime & ...
                    rips{ian}{day}{epoch}{1}.endtime<winendtime);
                ripsinwinTimes = [rips{ian}{day}{epoch}{1}.starttime(ripindsinwin) ...
                    rips{ian}{day}{epoch}{1}.endtime(ripindsinwin)];
                xs = ripsinwinTimes(:,1)-lfpstack(ian).ripStartTime(irip);
                xe = ripsinwinTimes(:,2)-lfpstack(ian).ripStartTime(irip);
                patch([xs'; xe'; xe'; xs'], repmat([0 0 yl(2) yl(2)]', 1, length(xe)),'y', ...
                    'FaceAlpha',.15, 'edgecolor','none');
                
                
                %% mu FR
                subaxis(2,Pp.nrow,2, 'Spacing', Pp.Spacing, ...
                    'SpacingVert', Pp.SpacingVert, 'SpacingHoriz', Pp.SpacingHoriz, 'Padding', ...
                    Pp.Padding, 'MarginLeft', Pp.MarginLeft, 'MarginRight', Pp.MarginRight,...
                    'MarginTop', Pp.MarginTop, 'MarginBottom', Pp.MarginBottom);
                
                irip_muinstaFR = squeeze(spikestack(ian).muinstantFR{day}{epoch}(irelrip,:,:))';
                imagesc(time, iripmuntrodes, 10*log10(flipud(irip_muinstaFR)));
                yticks(iripmuntrodes);
                yticklabels({num2str(flipud(iripmuntrodes))});
                colormap(jet)
                %             xlabel('time')
                title('mu firing rate', 'FontSize',8,'FontWeight','bold', 'FontName','Arial')
                patch([xs'; xe'; xe'; xs'], repmat([0 0 yl(2) yl(2)]', 1, length(xe)),'k', ...
                    'FaceAlpha',0, 'EdgeAlpha', .5, 'edgecolor','k', 'LineWidth', .05);
                
                %% plot LFP
                for t = 1:length(lfpstack(ian).lfptypes)
                    %                 axes(ha(t));
                    hold off
                    subaxis(2,Pp.nrow,t+2, 'Spacing', Pp.Spacing, ...
                        'SpacingVert', Pp.SpacingVert, 'SpacingHoriz', Pp.SpacingHoriz, 'Padding', ...
                        Pp.Padding, 'MarginLeft', Pp.MarginLeft, 'MarginRight', Pp.MarginRight,...
                        'MarginTop', Pp.MarginTop, 'MarginBottom', Pp.MarginBottom);
                    hold on
                    if t > 2 % zscore only the bandpass filtered traces
                        iripdata = squeeze(gridspacemat(t,:,:) + ...
                            zscore(double(lfpstack(ian).data{t}(irip, :, :)), [], 2))';
                    else
                        try
                            iripdata = squeeze(gridspacemat(t,:,:) + ...
                                lfpstack(ian).data{t}(irip, :, :))';
                        catch
                            iripdata = squeeze(int16(gridspacemat(t,:,:)) + ...
                                lfpstack(ian).data{t}(irip, :, :))';
                        end
                    end
                    p = plot(time, iripdata);
                    
                    ca1y = find(ismember(lfpstack(ian).ntrodes, ca1tets));
                    sz = zeros(length(lfpstack(ian).ntrodes), 3);
                    sz(ismember(lfpstack(ian).ntrodes, ca1y), 3) = 1; % make ca1 blue
                    refy = find(ismember(lfpstack(ian).ntrodes, reftet));
                    sz(lfpstack(ian).ntrodes==refy, 2) = 1; % make ref green
                    set(p, {'color'}, num2cell(sz, 2));
                    
                    yticks([squeeze(gridspacemat(t,1,:))]);
                    yticklabels({num2str(lfpstack(ian).ntrodes)});
                    title(sprintf('%s %s Hz',Fp.LFPtypes{t}, Fp.LFPrangesHz{t}), ...
                        'FontSize',8,'FontWeight','bold', 'FontName','Arial');
                    axis tight
                    yl = ylim;
                   
                    if strcmp(Fp.LFPtypes{t}, 'eeg')
                        % plot noise events in window with std label
                        noiseindsinwin = find(noise{ian}{day}{epoch}{1}.starttime>winstarttime & ...
                            noise{ian}{day}{epoch}{1}.endtime<winendtime);
                        noiseinwinTimes = [noise{ian}{day}{epoch}{1}.starttime(noiseindsinwin) ...
                            noise{ian}{day}{epoch}{1}.endtime(noiseindsinwin)];
                        nxs = noiseinwinTimes(:,1)-lfpstack(ian).ripStartTime(irip);
                        nxe = noiseinwinTimes(:,2)-lfpstack(ian).ripStartTime(irip);
                        plot([nxs'; nxe'],zeros(1,length(nxs)), '.r')
%                         plot(nxe,[zeros(1,length(nxe))-1], 'dk')
%                         patch([nxs'; nxe'; nxe'; nxs'], repmat([yl(1)-1 yl(1)-1 yl(1) yl(1)], 1, length(nxe)),'k', ...
%                             'FaceAlpha',1, 'edgecolor','none');
                        for i = 1:length(nxe)
                            text(nxe(i), yl(1), ...
                                sprintf('%.0f',noise{ian}{day}{epoch}{1}.maxthresh(noiseindsinwin(i))), ...
                                'FontName', 'Arial', 'FontSize',Pp.patchtxtsz);
                        end
                        
                    elseif strcmp(Fp.LFPtypes{t}, 'ripple')     
                        % plot refrip events in window with std label
                        refripindsinwin = find(refrip{ian}{day}{epoch}{1}.starttime>winstarttime & ...
                            refrip{ian}{day}{epoch}{1}.endtime<winendtime);
                        refripinwinTimes = [refrip{ian}{day}{epoch}{1}.starttime(refripindsinwin) ...
                            refrip{ian}{day}{epoch}{1}.endtime(refripindsinwin)];
                        rrxs = refripinwinTimes(:,1)-lfpstack(ian).ripStartTime(irip);
                        rrxe = refripinwinTimes(:,2)-lfpstack(ian).ripStartTime(irip);
                        plot([rrxs'; rrxe'],zeros(1,length(rrxs)), '.g')
%                         plot(rrxe,[zeros(1,length(nxe))-1], 'r')
%                         patch([rrxs'; rrxe'; rrxe'; rrxs'], repmat([yl(1)-1 yl(1)-1 yl(1) yl(1)]', 1, length(rrxe)),'b', ...
%                             'FaceAlpha',1, 'edgecolor','none');
                        for i = 1:length(rrxe)
                            text(rrxe(i), yl(1), ...
                                sprintf('%.0f',refrip{ian}{day}{epoch}{1}.maxthresh(refripindsinwin(i))), ...
                                'FontName', 'Arial', 'FontSize',Pp.patchtxtsz);
                        end
                    end
                    % plot rips
                    patch([xs'; xe'; xe'; xs'], repmat([0 0 yl(2) yl(2)]', 1, length(xe)),'y', ...
                        'FaceAlpha',.2, 'edgecolor','none');
                    for i = 1:length(xe)
                        text(xe(i), yl(2)+.5, ...
                            sprintf('%.0f',rips{ian}{day}{epoch}{1}.maxthresh(ripindsinwin(i))), ...
                            'FontName', 'Arial', 'FontSize',Pp.patchtxtsz);
                    end
                end
                %% single unit
                subaxis(2,Pp.nrow, [9 10], 'Spacing', Pp.Spacing, ...
                    'SpacingVert', Pp.SpacingVert, 'SpacingHoriz', Pp.SpacingHoriz, 'Padding', ...
                    Pp.Padding, 'MarginLeft', Pp.MarginLeft, 'MarginRight', Pp.MarginRight,...
                    'MarginTop', Pp.MarginTop, 'MarginBottom', Pp.MarginBottom);
                
                irip_suspikes = squeeze(spikestack(ian).suspikes{day}{epoch}(irelrip,:,:))';
                iripsu_ntcl = spikestack(ian).sukeys{day}{epoch}(:,[3 4]);
                %             irip_muspikes(irip_muspikes>0) = 1;
                %             f1 = imagesc(time, iripmuntrodes, irip_muspikes);
                %             colormap(1-gray);
                [sx, sy] = find(irip_suspikes');
                ca1y = find(ismember(iripsu_ntcl(:,1), ca1tets));
                sz = zeros(length(sy), 3);
                sz(ismember(sy, ca1y), 3) = 1;
                f1 = scatter(sx/1000-1.001, sy,10,sz,'+'); % could ad 'z' data for coloring
                %             f1.MarkerFaceAlpha = 0.1;
                f1.MarkerEdgeAlpha = .4;
                %             f1.MarkerEdgeColor = [0 0 0];
                yticks(1:length(iripsu_ntcl(:,1)));
                %             yticklabels(, 'FontSize', 6);
                set(gca,'YTickLabel',{num2str(iripsu_ntcl)},'fontsize',6)
                %             set(gca,'YTickLabelMode','auto')
                ylabel('ntrode - cluster','FontSize',8,'FontWeight','bold', 'FontName','Arial')
                xticks([-1 -.5 0 .5 1]);
                xlabel('time','FontSize',8,'FontWeight','bold', 'FontName','Arial');
                title('su', 'FontSize',8,'FontWeight','bold', 'FontName','Arial')
                axis tight
                xlim([-1 1])
                yl = ylim;
                patch([xs'; xe'; xe'; xs'], repmat([0 0 yl(2) yl(2)]', 1, length(xe)),'y', ...
                    'FaceAlpha',.15, 'edgecolor','none');
                
                %% 2d position
                %             axes(ha(t+3));
                subaxis(2,Pp.nrow,[11 12], 'Spacing', Pp.Spacing, ...
                    'SpacingVert', Pp.SpacingVert, 'SpacingHoriz', Pp.SpacingHoriz, 'Padding', ...
                    Pp.Padding, 'MarginLeft', Pp.MarginLeft, 'MarginRight', Pp.MarginRight,...
                    'MarginTop', Pp.MarginTop, 'MarginBottom', Pp.MarginBottom);
                
                posidx = knnsearch(pos{ian}{day}{epoch}.data(:,1), ...
                    lfpstack(ian).ripStartTime(irip));
                posirip = pos{ian}{day}{epoch}.data(posidx,:);
                xpos = pos{ian}{day}{epoch}.data(:,6);
                ypos = pos{ian}{day}{epoch}.data(:,7);
                p = scatter(xpos, ypos, 1, 'filled');
                p.MarkerFaceColor = [0 0 0];
                p.MarkerFaceAlpha = .1;
                hold on
                
                s = posidx - 300;
                if s<1; s=1; end
                e = posidx + 300;
                if e > length(xpos(:,1)); e = length(xpos(:,1)); end
                
                ixpos = xpos(s:e,1);
                iypos = ypos(s:e,1);
                jt = cool(length(iypos));
                be = scatter(ixpos, iypos, 10, jt);
                %             af = scatter(pos{ian}{day}{epoch}.data(posidx:e,6), ...
                %                 pos{ian}{day}{epoch}.data(posidx:e,7), 5, 'filled');
                %             be.MarkerFaceColor = [0 .5 1];
                %             be.MarkerFaceAlpha = 0;
                %             af.MarkerFaceColor = [.5 0 1];
                %             af.MarkerFaceAlpha = .5;
                plot(posirip(6), posirip(7),'.r', 'MarkerSize', 30)
                title('position', 'FontSize',8,'FontWeight','bold', 'FontName','Arial')
                xlabel('cm');
                axis tight
                xl = xlim;
                yl = ylim;
                if xl(2) > 250 || xl(1) < 0
                    xlim([50 250])
                end
                if yl(1) > 250 || yl(1) < 0
                    ylim([20 150])
                end
                linposidx = knnsearch(linpos{ian}{day}{epoch}.statematrix.time(:,1), ...
                    lfpstack(ian).ripStartTime(irip));
                wellexitenter = linpos{ian}{day}{epoch}.statematrix.wellExitEnter(linposidx,:);
                %             linpos{ian}{day}{epoch}.statematrix.segmentIndex(linposidx,1)
                %             sc = scatter(linpos{ian}{day}{epoch}.wellSegmentInfo.wellCoord(wellexitenter,1) ...
                %                 ,linpos{ian}{day}{epoch}.wellSegmentInfo.wellCoord(wellexitenter,2), '.');
                wsx = linpos{ian}{day}{epoch}.wellSegmentInfo.wellCoord(wellexitenter(1),1);
                wsy = linpos{ian}{day}{epoch}.wellSegmentInfo.wellCoord(wellexitenter(1),2);
                wex = linpos{ian}{day}{epoch}.wellSegmentInfo.wellCoord(wellexitenter(2),1);
                wey = linpos{ian}{day}{epoch}.wellSegmentInfo.wellCoord(wellexitenter(2),2);
                text(wsx, wsy,'\wedge','FontSize',20,'FontWeight','bold', ...
                    'Color', [0 .6 0], 'horizontalAlignment', 'center')
                text(wex, wey,'\vee','Fontsize',20,'FontWeight','bold',...
                    'Color', [.6 .3 .6], 'horizontalAlignment', 'center')
                %             linpos{ian}{day}{epoch}.statematrix.segmentIndex(linposidx,:);
                hold off
                
                %% Linear position
                % plot animal's lindist from center well and colorize by track segment
                subaxis(2,Pp.nrow,[13 16], 'Spacing', Pp.Spacing, ...
                    'SpacingVert', Pp.SpacingVert, 'SpacingHoriz', Pp.SpacingHoriz, 'Padding', ...
                    Pp.Padding, 'MarginLeft', Pp.MarginLeft, 'MarginRight', Pp.MarginRight,...
                    'MarginTop', Pp.MarginTop, 'MarginBottom', Pp.MarginBottom);
                
                postime = linpos{ian}{day}{epoch}.statematrix.time;
                traj = linpos{ian}{day}{epoch}.statematrix.traj;
                segmentnum = linpos{ian}{day}{epoch}.statematrix.segmentIndex;
                lindist = linpos{ian}{day}{epoch}.statematrix.linearDistanceToWells(:,1);
                plot(postime,lindist,'linewidth',1,'color',[.8 .8 .8]); %background pos
                hold on
                
                clr1 = [.6 .5 .2; .8 .5 .8; .6 0 .6; 0 0.4980 0.3255; 0 0.8 0.5216];
                %             line(postime, lindist, segmentnum);
                for iseg = 1:5
                    inds = [];
                    inds = double(segmentnum == iseg);
                    inds(find(inds==0))=nan;
                    plot(postime.*inds,lindist.*inds,'-','linewidth',2,'color',clr1(iseg,:));
                    hold on;
                end
                
                %% performance patches over linpos
                trialIO = behave{ian}.BehaveState.statechanges{day}{epoch}.statechangeseq;
                outBcol = 9; lastTimecol= 2; currTimecol=3; corrcol=7; timeportoutcol=10;
                outBCorr = trialIO((trialIO(:,outBcol)==1 & trialIO(:,corrcol)==1), ...
                    [lastTimecol currTimecol corrcol timeportoutcol]);
                outBMist = trialIO((trialIO(:,outBcol)==1 & trialIO(:,corrcol)==0), ...
                    [lastTimecol currTimecol corrcol timeportoutcol]);
                inBCorr = trialIO((trialIO(:,outBcol)==0 & trialIO(:,corrcol)==1), ...
                    [lastTimecol currTimecol corrcol timeportoutcol]);
                inBMist = trialIO((trialIO(:,outBcol)==0 & trialIO(:,corrcol)==0), ...
                    [lastTimecol currTimecol corrcol timeportoutcol]);
                clr2 = ([.5 .5 1; .5 .7 .7; .3 .8 .3; .5 .5 1; .5 .7 .7; .3 .8 .3;]);
                yl = ylim;
                xl = xlim;
                for itrial = 1:length(outBCorr(:,1)) %patch correct outbound trials
                    patch([outBCorr(itrial,1) outBCorr(itrial,2) outBCorr(itrial,2) ...
                        outBCorr(itrial,1)],  [yl(2) yl(2) yl(1) yl(1)] , 'y', 'edgecolor', ...
                        'none', 'FaceAlpha',.2);
                end
                for itrial = 1:length(outBMist(:,1)) %patch incorrect outbound trials
                    patch([outBMist(itrial,1) outBMist(itrial,2) outBMist(itrial,2) ...
                        outBMist(itrial,1)],  [yl(2) yl(2) yl(1) yl(1)] , 'k', 'edgecolor', ...
                        'none', 'FaceAlpha',.1);
                end
                for itrial = 1:length(outBCorr(:,1)) %is a portOUT
                    plot(outBCorr(itrial,4), yl(2)-1, 'd', 'linewidth', 1, 'MarkerEdgeColor', ...
                        'k','MarkerFaceColor', 'y'); hold on;
                end
                % inbound patches
                for itrial = 1:length(inBCorr(:,1)) %patch correct outbound trials
                    patch([inBCorr(itrial,1) inBCorr(itrial,2) inBCorr(itrial,2) ...
                        inBCorr(itrial,1)],  [yl(2) yl(2) yl(1) yl(1)] , 'y', 'edgecolor', ...
                        'none', 'FaceAlpha',.2);
                end
                for itrial = 1:length(inBMist(:,1)) %patch incorrect outbound trials
                    if inBMist(itrial,1) < 5 % ignore any patches within first 5 seconds
                        continue
                    end
                    patch([inBMist(itrial,1) inBMist(itrial,2) inBMist(itrial,2) ...
                        inBMist(itrial,1)],  [yl(2) yl(2) yl(1) yl(1)] , 'k', 'edgecolor', ...
                        'none', 'FaceAlpha',.1);
                end
                for itrial = 1:length(inBCorr(:,1)) %is a portOUT
                    plot(inBCorr(itrial,4), yl(1)-2, 'd', 'linewidth', 1, 'MarkerEdgeColor', ...
                        'k','MarkerFaceColor', 'y'); hold on;
                end
                
                % plot all rips in ep over linpos
                de_ripidx = find(ismember([lfpstack(ian).day lfpstack(ian).epoch], ...
                    [day epoch], 'rows'));
                de_ripidxlinpos = knnsearch(postime, lfpstack(ian).ripStartTime(de_ripidx));
                f = scatter(lfpstack(ian).ripStartTime(de_ripidx), lindist(de_ripidxlinpos), ...
                    10, 'filled');
                f.MarkerFaceColor = 'k';
                % highlight current rip
                hold on
                linposripidx = knnsearch(postime, lfpstack(ian).ripStartTime(irip));
                plot(lfpstack(ian).ripStartTime(irip), lindist(linposripidx), ...
                    '.r', 'MarkerSize', 30)
                title('linpos + performance', 'FontSize',8,'FontWeight','bold', 'FontName', ...
                    'Arial')
                xlabel('time s');
                axis tight
                
                % plot noise and refrip events over linpos
                f = scatter(noise{ian}{day}{epoch}{1}.starttime, zeros(length(noise{ian}{day}{epoch}{1}.starttime),1)+yl(1)+1, ...
                    1, 'filled');
                f.MarkerFaceColor = [1 0 0];
                f = scatter(refrip{ian}{day}{epoch}{1}.starttime, zeros(length(refrip{ian}{day}{epoch}{1}.starttime),1)+yl(1)+5, ...
                    1, 'filled');
                f.MarkerFaceColor = [0 1 0];
               
                % plot prefirst and extendedimmobility period patches
                pre_s = firstwell{ian}{day}{epoch}.prefirst_list(:,1);
                pre_e = firstwell{ian}{day}{epoch}.prefirst_list(:,2);
                for i = 1:length(pre_s)
                    patch([pre_s'; pre_e'; pre_e'; pre_s'], repmat([0 0 yl(2) yl(2)]', 1, length(pre_s)),'b', ...
                        'FaceAlpha',.2, 'edgecolor','none');
                end
                
                imm_s = immobile{ian}{day}{epoch}.immobile_list(:,1);
                imm_e = immobile{ian}{day}{epoch}.immobile_list(:,2);
                for i = 1:length(imm_s)
                    patch([imm_s'; imm_e'; imm_e'; imm_s'], repmat([0 0 yl(2) yl(2)]', 1, length(imm_s)),'r', ...
                        'FaceAlpha',.2, 'edgecolor','none');
                end
                
                %% super
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('%s %d %d %s %s rip%d start%dms', animal, day, epoch, ...
                    Fp.paths.filenamesave(5:end), Fp.epochEnvironment, irip, ...
                    round(1000*lfpstack(ian).ripStartTime(irip)));
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center','FontSize', 12);
                %% ---- pause, save figs ----
                if pausefigs
                    pause
                end
                if savefigs
                    save_figure(Fp.paths.figdirectory, 'riptrig_wExcl', sprtit)
                    close all
                end
            end
        end
    end
end