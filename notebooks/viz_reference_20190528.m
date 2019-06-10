


% Diagnose and fix the issue with referencing on many of the days of data
% the average LFP responses across animals and across days within animals
% is quite variable, my guess is reference bleed that is introducing big
% artifacts across the channels
% multiunit spikes are largely immune to this (as long as the mu is not on
% the reference tetrode) so I could use that as a sanity check as long as I
% confirm that multiunit on the reference is not an issue..
% alternatively, maybe I could just go in an remove the mu spikes on the
% ref from the rest.
% ACTUALLY i think that i gave Mountainsort unreferenced MDA so this
% would't be an issue, at the cost at potentially introducing more noise
% into the clustering pipeline.. which might change the cluster metrics and
% even act too conservatively

% TODO: visualize the reference on the average rip-trig lfp traces, both with
% unreferenced (eeggnd) and referenced lfp (eeg)
% - does it look like there's a dominant theta rhythim, or any other rhythm
% being introduced into the LFP?
% - which tetrode are quiet except for large artifactual fluctuations on the
% unreferenced lfp and could potentially be used as the reference tet?
% - for the noisy tetrodes, are there other channels from this tetrode that
% are not noisy? to resovle this I would need to either load the rec file back up into
% trodes and observe there.. or extract all the channels with the rec
% extractor then convert into filterframework, and then visualize here...
%
% NOOOO the lfp problem is hard.. so use multiunit instead for now..
% for D10.. there's definitely multiunit spikes on the reference channel..

% wait no i didn't use a ref for clustering right?? confirm this now
% ok it looks like i did use the default setting for 'userefs' of exportmda
% it was the lfp that i extracted specifically with 'userefs' 0 so that I
% could just reference myself later (eeggnd to eeg)
% so it's very likely i guess that spikes designated as 'multiunit' on the
% reference tetrode made it into the multiunit cluster of the other
% tetrodes.. so how should i deal with this? There's probably a bit of
% wiggle time between the exact spike times of reference tet spikes on the
% others, since it they would be at different baselines... but at most this
% would be the duration of a spike.. so i guess i could just subtract out
% multiunit spikes within the same ms on the reference to those on the
% other tetrodes.. although it's definately possible and likely that this
% will erase a lot of 'real' spikes that happend within 1ms of the
% reference spike... so i have to be careful to not introduce structure
% into the multiunit by subtracting out the multiunit spikes..

% so replot D10 riptrigmu but with tetrode 14 multiunit spikes (per ms)
% removed..
% what was the script i used to make those other plots?
% riptrigrasters_alltetsoverdays_20190514.m

% basically need to treat the multiunit firing rate as lfp and do phase
% reset analysis, cross-freq-coupling, etc

% fuck ok this needs to be fixed since previously i was
% stacking per day, across epochs.. and now i'm not..
% should i just start doing what mike cohen does? using
% a single cube of data for all the riptrigmu (rip x
% time x ntrode) ?? i think this is a good idea.. it's
% basically the same idea that xarray employs.. and it
% makes things a lot easier in terms of data management
% to have things in a single multidim and multindexed
% array.. so i should create a function that takes the
% output of the filter framework analysis functions
% (dfa...) and reformat into an array block.. that
% would also make life a lot easier when converting
% back and forth with python.. to what extent do i need
% to standardize? i should have a single struct with
% the .data field.. then any other fields as a tag or as an index
% vector the length a matching dim of the corresponding
% data dimension.. i.e. ripplenum should be
% something like size (1089 x 1 x 1) whereas ntrode should be
% size (1 x 1 x 30) and time should be (1 x 1001 x 1)
% i should have another 'fields' field that name the
% data dims {'ripplenum', 'time', 'ntrode}
% i could also have a day,ep vector or an epochnum vec
% i guess i could also have across all epochtypes and
% then have an epochtype vector to slice into the
% result data.. i like this a lot as an
% analysis/plotting standard format to get the analysis
% results into.. of would be kind of weird
% to do for 2d or gridded data like place fields,
% remember that i need to focus on getting riptriggered multiunit
% firing rate rasters..
% thought about this again... i can't just use
% multiunit firing because i want to estimate phase
% reset.. which is dependent on frequency, which is
% dependent on smoothing/bandwidth parameters

% but i still really want to use spikes...
% maybe just go forward with multiunit spikes right now..
% and what? no i need to fix rip triggreed lfp.

% PLOT UNREFERENCE AND REF LFP

Fp.animals = {'D12', 'JZ1'}; 
Fp.filtfunction = 'dfa_riptriglfp';
epochEnvironment = 'sleep';
Fp.add_params = {epochEnvironment, 'ripples', '<1cm/s'};
me = animaldef('demetris');

run_ff = 0;
savedata = 0;
load_data = 0;
stack_data = 1;
plotfigs = 1;
pausefigs = 0;
savefigs = 1;
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
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
        ['_', Fp.epochEnvironment, '_', 'refUnref'])
end
%% load data
if load_data
    F = load_filter_output(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', ['_', Fp.epochEnvironment, '_', 'refUnref']);
end
%% combine per ntrode
if stack_data
    datastack = stack_riptriglfp(F);
end
%% plot
if plotfigs
    Pp = load_plotting_params('riptriglfp_waterfall');
    for ian = 1:numel(Fp.animals)
        animal = Fp.animals{ian};
        anidx = find(strcmp({datastack.animal}, animal));
        ntrodes = datastack(anidx).ntrodes;
        for t = 1:length(datastack(anidx).lfptypes)
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
            ha = tight_subplot(2, ceil(max(ntrodes)/2), [.05 .005],[.05 .1],[.05 .05]);
            for nti = 1:length(ntrodes)
                ntrode = ntrodes(nti);
%                 subaxis(2, length(ntrodes)/2, nti);
                axes(ha(ntrode));
                d = datastack(anidx).data{t}(:,:,nti);
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
                % day/epoch lines
                dayind = find(diff(datastack(anidx).dayeps(:,1))>0);
                epbounds = cumsum(datastack(anidx).numrips_perep);
                daybounds = epbounds(dayind);
                line([-Pp.pwin(1) Pp.pwin(2)], [epbounds'; epbounds'], 'color', [.9 .9 .9])
                line([-Pp.pwin(1) Pp.pwin(2)], [daybounds'; daybounds'], 'color', [0 0 0])
                line([0 0], [1 size(z,1)], 'color', [0 0 0], 'linestyle', '--')
                title(sprintf('%d', ntrode))
                if nti ~= length(ntrodes)/2+1
                    set(gca, 'ytick', []);
                    set(gca, 'xtick', []);
                end
            end
            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('%s %s %s %s', animal, datastack(anidx).lfptypes{t}, ...
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


%%
% for t = 1:length(datastack(ian).data)+1
%                 subaxis(1,length(datastack(ian).data)+1, t);
%                 if t <3
%                     d = datastack(ian).data{t}(:,:,ntrode);
%                      %trim and nan zscore
%                     mididx = ceil(size(d,2)/2); % right now assumes center is rip start
%                     srate = 1500;
%                     d = double(d(:,mididx-(Pp.pwin(1)*srate):mididx+(Pp.pwin(2)*srate)));
%                     ptime = time(mididx-(Pp.pwin(1)*srate):mididx+(Pp.pwin(2)*srate));
%                     m = nanmean(d,2);
%                     s = nanstd(d, [], 2);
%                     z = (d-m)./s;
%                     % plot
%                     imagesc(ptime, 1:size(z,1), z)
%                     % day/epoch lines
%                     dayind = find(diff(datastack(ian).dayeps(:,1))>0);
%                     epbounds = cumsum(datastack(ian).numrips_perep);
%                     daybounds = epbounds(dayind);
%                     line([-Pp.pwin(1) Pp.pwin(2)], [epbounds'; epbounds'], 'color', [.9 .9 .9])
%                     line([-Pp.pwin(1) Pp.pwin(2)], [daybounds'; daybounds'], 'color', [0 0 0])
%                     line([0 0], [1 size(z,1)], 'color', [0 0 0], 'linestyle', '--')
%                     title(datastack(ian).lfptypes{t})
%                 else
%                     refntrode = 14;
%                     d = datastack(ian).data{2}(:,:,ntrode);
%                     r = datastack(ian).data{2}(:,:,refntrode);
%                     d = d+r;
%                      %trim and nan zscore
%                     mididx = ceil(size(d,2)/2); % right now assumes center is rip start
%                     srate = 1500;
%                     d = double(d(:,mididx-(Pp.pwin(1)*srate):mididx+(Pp.pwin(2)*srate)));
%                     ptime = time(mididx-(Pp.pwin(1)*srate):mididx+(Pp.pwin(2)*srate));
%                     m = nanmean(d,2);
%                     s = nanstd(d, [], 2);
%                     z = (d-m)./s;
%                     % plot
%                     imagesc(ptime, 1:size(z,1), z)
%                     % day/epoch lines
%                     dayind = find(diff(datastack(ian).dayeps(:,1))>0);
%                     epbounds = cumsum(datastack(ian).numrips_perep);
%                     daybounds = epbounds(dayind);
%                     line([-Pp.pwin(1) Pp.pwin(2)], [epbounds'; epbounds'], 'color', [.9 .9 .9])
%                     line([-Pp.pwin(1) Pp.pwin(2)], [daybounds'; daybounds'], 'color', [0 0 0])
%                     line([0 0], [1 size(z,1)], 'color', [0 0 0], 'linestyle', '--')
%                     title('ref+unref')
%
%                 end