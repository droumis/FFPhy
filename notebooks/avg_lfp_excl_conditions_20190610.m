
% need to plot heatrasters and averages of lfp using the exclusion
% timefilters/eventfilters, and have each rip tagged with condition labels
% so i can compare across conditions, such as correct/incorrect, etc.
% first i just need to create the average LFP plots per epoch..

% allnt riptriglfp heatrasters, 1 lfptype: viz_reference_20190528
% mean lfp 1lfptype waterfall: phasereset_20190524
% 1 ntrode all lfptypes: phasereset_20190521

load_lfp = 0;
stack_lfp = 0;
load_events = 0;
get_filter_vecs = 0;


plotfigs = 1;
plot_allntrodes = 0;
plot_all_lfptypes = 0;
plot_avg_traces = 1;
savefigs = 1;
pausefigs = 0;

animals = {'D10'};
add_params = {'wtrack'};
Fp.animals = animals;
Fp.add_params = add_params;
Fp.filtfunction = 'dfa_riptriglfp';
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);

%% load lfp
if load_lfp
    F = load_filter_output(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', ['_', Fp.epochEnvironment]);
end
%% stack LFP
if stack_lfp
    lfpstack = stack_riptriglfp(F);
end

% now that i have an lfp stack.. i want to create event vecs that tag each
% event, so that i can use that for inclusion/grouping

% Intervals:

% extended immobile
% prefirst well

% within 1 second of a large reference rip event (rr > 15std)
% within 1 second of a large noise event (noise > 15 std)
% within .5 sec of a ref rip bigger than the ca1 rip (rr > ca1r)

%% load events, infostructs, timefilters
if load_events
    for ian = 1:length(Fp.animals)
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntinfoAll = cellfetch(ntinfo, '', 'alltags', 1);
        rips{ian} = loaddatastruct(aninfo{2}, animal, 'ca1rippleskons');
        noise{ian} = loaddatastruct(aninfo{2}, animal, 'ca1noisekons');
        refrip{ian} = loaddatastruct(aninfo{2}, animal, 'refrippleskons');
        immobile{ian} = getextendedimmobility(aninfo{2},aninfo{3}, ...
            F(1).epochs{1});
        firstwell{ian} = getpriortofirstwell(aninfo{2},aninfo{3}, ...
            F(1).epochs{1});
    end
end

%% create vec with entry for each lfp rip event. save alongside lfpstack
if get_filter_vecs
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


%% plot
if plotfigs
    if plot_allntrodes
        Pp = load_plotting_params('riptriglfp_perLFPtype_allntrodes');
        for ian = 1:numel(Fp.animals)
            animal = Fp.animals{ian};
            anidx = find(strcmp({lfpstack.animal}, animal));
            ntrodes = lfpstack(anidx).ntrodes;
            for t = 1:length(lfpstack(anidx).lfptypes)
                %         time = lfpstack(ian).time;
                %% ---- init fig----
                if savefigs && ~pausefigs
                    close all
                    ifig = figure('Visible','off','units','normalized','position', ...
                        Pp.position);
                else
                    ifig = figure('units','normalized','position',Pp.position);
                end
                set(gcf,'color','white')
                %%
                ha = tight_subplot(2, ceil(max(ntrodes)/2), [.05 .005],[.05 .1],[.05 .05]);
                exclude_rips = any(lfpstack(anidx).filtervecs,2);
                scratchstack = lfpstack(anidx).data{t};
                scratchstack(~exclude_rips,:,:) = [];
                
                
                for nti = 1:length(ntrodes)
                    ntrode = ntrodes(nti);
                    %                 subaxis(2, length(ntrodes)/2, nti);
                    axes(ha(ntrode));
                    d = scratchstack(:,:,nti);
                    %trim and nan zscore
                    mididx = ceil(size(d,2)/2); % right now assumes center is rip start
                    srate = 1500;
                    d = double(d(:,mididx-(Pp.pwin(1)*srate):mididx+(Pp.pwin(2)*srate)));
                    ptime = Fp.time(mididx-(Pp.pwin(1)*srate):mididx+(Pp.pwin(2)*srate));
                    m = nanmean(d,2);
                    s = nanstd(d, [], 2);
                    z = (d-m)./s;
                    % plot
                    imagesc(ptime, 1:size(z,1), z)
                    colormap(jet)
                    % day/epoch lines
                    %                 dayind = find(diff(lfpstack(anidx).dayeps(:,1))>0);
                    %                 epbounds = cumsum(lfpstack(anidx).numrips_perep);
                    %                 daybounds = epbounds(dayind);
                    %                 line([-Pp.pwin(1) Pp.pwin(2)], [epbounds'; epbounds'], 'color', [.9 .9 .9])
                    %                 line([-Pp.pwin(1) Pp.pwin(2)], [daybounds'; daybounds'], 'color', [0 0 0])
                    line([0 0], [1 size(z,1)], 'color', [0 0 0], 'linestyle', '--')
                    title(sprintf('%d', ntrode))
                    if nti ~= length(ntrodes)/2+1
                        set(gca, 'ytick', []);
                        set(gca, 'xtick', []);
                    end
                end
                %% super
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('%s %s %s %s', animal, lfpstack(anidx).lfptypes{t}, ...
                    Fp.paths.filenamesave(5:end), Fp.epochEnvironment);
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
    if plot_all_lfptypes
        Pp = load_plotting_params('riptriglfp_allLFPtype_perntrode');
        for ian = 1:numel(Fp.animals)
            animal = Fp.animals{ian};
            anidx = find(strcmp({lfpstack.animal}, animal));
            ntrodes = lfpstack(anidx).ntrodes;
            dayep = [lfpstack(anidx).day lfpstack(anidx).epoch];
            for nt = 1:length(ntrodes)
                ntrode = ntrodes(nt);
                %% ---- init fig----
                if savefigs && ~pausefigs
                    close all
                    ifig = figure('Visible','off','units','normalized','position', ...
                        Pp.position);
                else
                    ifig = figure('units','normalized','position',Pp.position);
                end
                set(gcf,'color','white')
                for t = 1:length(lfpstack(anidx).lfptypes)
                    subplot(1,length(lfpstack(anidx).lfptypes), t)
                    %                 ha = tight_subplot(2, ceil(max(ntrodes)/2), [.05 .005],[.05 .1],[.05 .05]);
                    exclude_rips = any(lfpstack(anidx).filtervecs,2);
                    scratchstack = lfpstack(anidx).data{t};
                    scratchstack(~exclude_rips,:,:) = [];
                    d = scratchstack(:,:,nti);
%                     a = squeeze(d)';
%                     z = 10*log10(abs(hilbert(a))');
%                     imagesc(z)
                    
                    % trim and nan zscore
                    mididx = ceil(size(d,2)/2); % right now assumes center is rip start
                    srate = 1500;
                    d = double(d(:,mididx-(Pp.pwin(1)*srate):mididx+(Pp.pwin(2)*srate)));
                    ptime = Fp.time(mididx-(Pp.pwin(1)*srate):mididx+(Pp.pwin(2)*srate));
                    m = nanmean(d,2);
                    s = nanstd(d, [], 2);
                    z = (d-m)./s;
                    imagesc(ptime, 1:size(z,1), z)
                    colormap(parula)
                    % day/epoch lines
%                     lfpstack(anidx).dayeps(:,1)
                    exdayep = dayep(~exclude_rips,:);
                    daybounds = find(diff(exdayep(:,1)));
                    epbounds = find(abs(diff(exdayep(:,2))));
                    line([-Pp.pwin(1) Pp.pwin(2)], [epbounds'; epbounds'], 'color', [.9 .9 .9])
                    line([-Pp.pwin(1) Pp.pwin(2)], [daybounds'; daybounds'], 'color', [0 0 0])
                    line([0 0], [1 size(z,1)], 'color', [0 0 0], 'linestyle', '--')
                    title(sprintf('%s', lfpstack(anidx).lfptypes{t}))
                end
                %% super
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('heatraster riptriglfp pernt %s nt%d %s', animal, ntrode, ...
                    Fp.epochEnvironment);
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center','FontSize', 12);
                %                     spylabel = text(.02, .5, sprintf('ripnum'), 'Parent', sprtitleax, 'Units', ...
                %                         'normalized', 'Rotation', 90);
                %                     set(spylabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                %                         'horizontalAlignment', 'center', 'FontSize', 12);
                %
                %                     spxlabel = text(.5, .02, sprintf('time'), 'Parent', sprtitleax, 'Units', ...
                %                         'normalized');
                %                     set(spxlabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                %                         'horizontalAlignment', 'center', 'FontSize', 12);
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
    if plot_avg_traces
        Pp = load_plotting_params('riptriglfp_allLFPtype_perntrode');
        for ian = 1:numel(Fp.animals)
            animal = Fp.animals{ian};
            anidx = find(strcmp({lfpstack.animal}, animal));
            ntrodes = lfpstack(anidx).ntrodes;
            dayep = [lfpstack(anidx).day lfpstack(anidx).epoch];
            for nt = 1:length(ntrodes)
                ntrode = ntrodes(nt);
                %% ---- init fig----
                if savefigs && ~pausefigs
                    close all
                    ifig = figure('Visible','off','units','normalized','position', ...
                        Pp.position);
                else
                    ifig = figure('units','normalized','position',Pp.position);
                end
                set(gcf,'color','white')
                for t = 1:length(lfpstack(anidx).lfptypes)
                    subplot(3,length(lfpstack(anidx).lfptypes), [t t+length(lfpstack(anidx).lfptypes)])
                    %                 ha = tight_subplot(2, ceil(max(ntrodes)/2), [.05 .005],[.05 .1],[.05 .05]);
                    exclude_rips = any(lfpstack(anidx).filtervecs,2);
                    scratchstack = lfpstack(anidx).data{t};
                    scratchstack(~exclude_rips,:,:) = [];
                    d = scratchstack(:,:,nt);
%                     a = squeeze(d)';
%                     z = 10*log10(abs(hilbert(a))');
%                     imagesc(z)
                    
                    % trim and nan zscore
                    mididx = ceil(size(d,2)/2); % right now assumes center is rip start
                    srate = 1500;
                    d = double(d(:,mididx-(Pp.pwin(1)*srate):mididx+(Pp.pwin(2)*srate)));
                    ptime = Fp.time(mididx-(Pp.pwin(1)*srate):mididx+(Pp.pwin(2)*srate));
                    m = nanmean(d,2);
                    s = nanstd(d, [], 2);
                    z = (d-m)./s;
                    imagesc(ptime, 1:size(z,1), z)
                    colormap(parula)
                    % day/epoch lines
%                     lfpstack(anidx).dayeps(:,1)
                    exdayep = dayep(~exclude_rips,:);
                    daybounds = find(diff(exdayep(:,1)));
                    epbounds = find(abs(diff(exdayep(:,2))));
                    line([-Pp.pwin(1) Pp.pwin(2)], [epbounds'; epbounds'], 'color', [.9 .9 .9])
                    line([-Pp.pwin(1) Pp.pwin(2)], [daybounds'; daybounds'], 'color', [0 0 0])
                    line([0 0], [1 size(z,1)], 'color', [0 0 0], 'linestyle', '--')
                    title(sprintf('%s', lfpstack(anidx).lfptypes{t}))
                    
                    
                    subplot(3,length(lfpstack(anidx).lfptypes), t+(length(lfpstack(anidx).lfptypes)*2))
                    plot(ptime, nanmean(d))
                end
                
                %% super
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('heatraster riptriglfp pernt %s nt%d %s', animal, ntrode, ...
                    Fp.epochEnvironment);
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center','FontSize', 12);
                %                     spylabel = text(.02, .5, sprintf('ripnum'), 'Parent', sprtitleax, 'Units', ...
                %                         'normalized', 'Rotation', 90);
                %                     set(spylabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                %                         'horizontalAlignment', 'center', 'FontSize', 12);
                %
                %                     spxlabel = text(.5, .02, sprintf('time'), 'Parent', sprtitleax, 'Units', ...
                %                         'normalized');
                %                     set(spxlabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                %                         'horizontalAlignment', 'center', 'FontSize', 12);
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

