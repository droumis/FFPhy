
% need to fix the referencing for D12, JZ1, JZ4.. 

% I need to confirm that there aren't SU spikes on the new ref channels. 
% Once i confirm a good new ref, I need to remake the referenced lfp
% then rerun riptrig lfp
% then plot heatrasters

% :: D12 :: 
% current ref is 11
% the problem is that i can see a lfp response in the unreferenced lfp.. 
% ref 5 seems flatter. but does it have spikes?
% how to check the validity of a new ref?

%{ 
- check my notes about that ntrode over days
- check the preprocess script for the animal to make sure it matches what i
intended as the ref in my notes
- check the histology figures to make sure the new ref is in CC
- load the cellinfo, check if that ntrode has any non-mu
- look at the raw traces from the rip plots to see how it looks. 
%}


% i don't think nt5 is a good ref for D12
% it's in POR
% the only tet near CC is 11.. but it is really close to deep MEC
% are there any quiet CA1 tets? maybe 18 or 27

% - ***my notes say that the w track ports were changed for D12 from D10.. see
% sheet15 in the exp log. is this accounted for in the
% linpos/preprocessing??? Also D10 day 9 the ports are different.. so need
% to check on that for D10.. 

% nt 18 looks like a good option as a reference from my notes.. 
% how does it look SU-wise?

if 0
    andef = animaldef('D12');
    cellinfo = loaddatastruct(andef{2}, andef{3}, 'cellinfo');
    su = evaluatefilter(cellinfo, 'isempty($tags)');
    %%
    checknt = 11;
    hassu = ~isempty(su(su(:,3)==checknt,:));
    ntsWithSu = unique(su(:,3));
    
end
% cool.. no single units on nt 18.. how does the ripple filtered lfp look?
% since the existing rip filtered lfp is referenced.. i'll need to check
% there's a slight bump in the riptrig lfp of 18.. although that's somewhat
% to be expected in the CC... nt 5 looks flat and also has no SU, but i'm
% quite sure it's not in the CC. 

% now that I've re-referenced D12 LFP using nt 18 as the reference.. plot
% the heatrasters

srate = 1500;

Fp.animals = {'D12'}; 
Fp.filtfunction = 'dfa_riptriglfp';
epochEnvironment = 'wtrack';
Fp.add_params = {epochEnvironment};
me = animaldef('demetris');
use_filters = {'firstwell', 'noise'};

run_ff = 0;
savedata = run_ff;
load_data = 0;
stack_data = 0;
load_events = 0;
make_filter_vecs = 0;

plotfigs = 1;
pausefigs = 0;
savefigs = 1;
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
filetail = 'ref18Unref'; % mini description to prevent overwriting something else
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
        ['_', Fp.epochEnvironment, '_', filetail])
end
%% load data
if load_data
    F = load_filter_output(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', ['_', Fp.epochEnvironment, '_', filetail]);
end
%% stack
if stack_data
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

%% create vec with entry for each lfp rip event. save alongside lfpstack
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

%% plot
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
            %         time = datastack(ian).time;
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
%             ha = tight_subplot(2, ceil(max(ntrodes)/2), [.05 .005],[.05 .1],[.05 .05]);
            
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
                numeps = length(de(:,1));
                days = unique(de(:,1));
                numdays = length(days);
                daybounds = find(diff(exdayep(:,1)));
                epbounds = find(abs(diff(exdayep(:,2))));
                
                exclude_rips = any(lfpstack(anidx).filtervecs(:,use_filts),2);
                ntd = squeeze(double(lfpstack(anidx).data{t}(:,:,nti)));
                %d = double(d(:,mididx-(Pp.pwin(1)*srate):mididx+(Pp.pwin(2)*srate)));
                excld = ntd;
                excld(exclude_rips,:) = [];

                d = lfpstack(anidx).data{t}(:,:,nti);
                %trim and nan zscore
                mididx = ceil(size(excld,2)/2); % right now assumes center is rip start
                ptime = Fp.time(mididx-(Pp.pwin(1)*srate):mididx+(Pp.pwin(2)*srate));
                m = nanmean(excld,2);
                s = nanstd(excld, [], 2);
                z = (excld-m)./s;
                imagesc(ptime, 1:size(z,1), z)
                colormap(parula)
                % day/epoch lines
%                 dayind = find(diff(lfpstack(anidx).dayeps(:,1))>0);
%                 epbounds = cumsum(lfpstack(anidx).numrips_perep);
%                 daybounds = epbounds(dayind);
                line([-Pp.pwin(1) Pp.pwin(2)], [epbounds'; epbounds'], 'color', [.9 .9 .9])
                line([-Pp.pwin(1) Pp.pwin(2)], [daybounds'; daybounds'], 'color', [0 0 0])
                line([0 0], [1 size(z,1)], 'color', [0 0 0], 'linestyle', '--')
                title(sprintf('%d', ntrode), 'FontSize',Pp.FontS, ...
                    'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
                caxis([-1 1])
                if nti ~= 1
                    set(gca, 'ytick', []);
                    set(gca, 'xtick', []);
                end
            end
            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('%s %s %s %s %s', animal, lfpstack(anidx).lfptypes{t}, ...
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

