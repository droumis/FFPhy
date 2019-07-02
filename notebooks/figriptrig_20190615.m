%{

Figure 1:
% Sen Chen replication.. but not just novelty, also during learning.

% "Theta Amplitude in MEC During awake ca1 SWR's"
- show example of riptrig CA1, MEC LFP and SU spiking
- show average riptrig CA1, MEC LFP and SU spiking
- show rip trig theta power

Figure 2:
% "Theta Phase in MEC During awake ca1 SWR's"
- show rip trig theta phase
- show CA1, MEC ITPC over freq, time

Figure 3:
% "Theta Increases according to novelty and learning rate"
- % show

Figure 4:
% "Theta reset at time of SWRs"
- show

Figure 5:
% "Theta reset at time of SWRs across MEC grid spacing modules"
-
%}

% but first plot traces that go along with the all ntrodes heatrasters
% start with D10, D12, then when D13, JZ1, JZ4 is done, use those too

run_ff = 1;
savedata = run_ff;

load_lfp = 0;
stack_lfp = 0;
load_events = 0;
make_filter_vecs = 0;

plotfigs = 0;
plot_heatrast_traces_perep_allnt = 0;

savefigs = 1;
pausefigs = 0;

use_filters = {'firstwell', 'noise'};
animals = {'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'}; %, 'D13'}; % next run ff jz1, jz3, jz4
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

if 0
    for ian = 1:length(Fp.animals)
        aninfo = animaldef(Fp.animals{ian});
%         refnt = evaluatefilter(ntinfo{1}, 'isequal($area, ''ref'')');
%         mecnt = evaluatefilter(ntinfo{1}, 'isequal($area, ''mec'') || isequal($subarea, ''mec'')');
%         mecnt = unique(mecnt(:,3));
%         mecntidx = find(ismember(lfpstack(ian).ntrodes,mecnt));
%         mecrefnt = evaluatefilter(ntinfo{1}, 'isequal($area, ''ref'') && isequal($subarea, ''mec'')');
%         mecrefnt = unique(mecrefnt(:,3));
%         mecrefntidx = find(lfpstack(ian).ntrodes == mecrefnt);
%         ca1nt = evaluatefilter(ntinfo{1}, 'isequal($area, ''ca1'') || isequal($subarea, ''ca1'')');
%         ca1nt = unique(ca1nt(:,3));
%         ca1ntidx = find(ismember(lfpstack(ian).ntrodes,ca1nt));
%         ca1refnt = evaluatefilter(ntinfo{1}, 'isequal($area, ''ref'') && isequal($subarea, ''ca1'')');
%         ca1refnt = unique(ca1refnt(:,3));
%         ca1refntidx = find(lfpstack(ian).ntrodes == ca1refnt);
        newstack = struct;
%         newstack(ian).data{1} = lfpstack(ian).data{1};
        for t = 1
            mecrefd1 = lfpstack(ian).data{t}(:,:,1:4) - lfpstack(ian).data{t}(:,:,5);
            mecrefd2 = lfpstack(ian).data{t}(:,:,6:15) - lfpstack(ian).data{t}(:,:,5);
            ca1refd = lfpstack(ian).data{t}(:,:,16:30) - lfpstack(ian).data{t}(:,:,5);
            newstack(ian).data{t} = cat(3, mecrefd1, lfpstack(ian).data{t}(:,:,5), mecrefd2, ca1refd);
        end
    end
end
%%

if plotfigs
    Pp = load_plotting_params({'defaults', 'riptriglfp_perLFPtype_allntrodes'});
    for ian = 1:numel(Fp.animals)
        animal = Fp.animals{ian};
        anidx = find(strcmp({lfpstack.animal}, animal));
        ntrodes = lfpstack(anidx).ntrodes;
        dayep = [lfpstack(anidx).day lfpstack(anidx).epoch];

        use_filts = find(any(cell2mat(cellfun(@(x) strcmp(x, lfpstack(anidx).filterfields), ...
            use_filters, 'un', 0)), 2));
        exclude_rips = any(lfpstack(anidx).filtervecs(:,use_filts),2);

        invalidtets = evaluatefilter(ntinfo{ian}, 'isequal($valid, ''no'')');
        invalidtets = unique(invalidtets(:,3));
        for t = 1:length(lfpstack(anidx).lfptypes)
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
                
                subaxis(2,ceil(max(ntrodes)/2), nti, 'SpacingVert', Pp.SpVt, ...
                'SpacingHoriz', Pp.SpHz, 'MarginLeft', Pp.MgLt, 'MarginRight', ...
                Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);

                exdayep = dayep(~exclude_rips,:);
                de = unique(exdayep, 'rows');
                daybounds = find(diff(exdayep(:,1)));
                epbounds = find(abs(diff(exdayep(:,2))));

                excld_stack = squeeze(double(lfpstack(anidx).data{t}(~exclude_rips,:,nti)));
                
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
                title(sprintf('%d', ntrode), 'FontSize',Pp.FontS, 'FontWeight',Pp.FontW, ...
                    'FontName', Pp.FontNm)
                caxis([-1 1])
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
            sprtit = sprintf('%s %s %s %s %s 20190615', animal, lfpstack(anidx).lfptypes{t}, ...
                Fp.paths.filenamesave(5:end), Fp.epochEnvironment, filetail);
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





