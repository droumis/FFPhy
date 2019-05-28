


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

animals = {'D10'};
filtfunction = 'dfa_riptriglfp';
env = 'wtrack';
me = animaldef('demetris');

loaddata = 1;
combine_per_ntrode = 0;
plotfigs = 0;
pausefigs = 0;
savefigs = 0;
%% load data
paths = make_paths(filtfunction, env);
if loaddata
    data = load_filter_output(paths.filtOutputDirectory, paths.filenamesave, ...
        animals);
end
%% combine per ntrode
if combine_per_ntrode
%     psth = combine_riptrigspiking_perntrode(data);
    datastack = stack_riptriglfp(data);
end

%% plot
if plotfigs
    Pp = load_plotting_params('all_nts_days_riptrigspiking');
    for ian = 1:numel(animals)
        animal = animals{ian};
        fprintf('%s\n',animal);
        %% ---- init fig----
        if savefigs && ~pausefigs
            close all
            ifig = figure('Visible','off','units','normalized','position', ...
                Pp.position);
        else
            ifig = figure(2);%'units','normalized','position',Pp.position);
        end
        set(gcf,'color','white')
        allnts = [data{ian}.ntrode];
        ntrodes = allnts;
        idx = cell2mat({data{ian}(1).mudata.index}');
        days = unique(idx(:,1), 'stable');
        numntrodes = numel(ntrodes);
        numdays = numel(days);
        ha = tight_subplot(numntrodes, numdays,[.01 .005],[.05 .05],[.05 .05]);
        for int = 1:numntrodes
            nt = ntrodes(int);
            ntind = find(allnts==nt);
            for idy = 1:numdays
                axes(ha(idy+((int-1)*numdays)));
                if ~isempty(data{ian}(ntind).mudata)
                    
                    day_mu_psth = cell2mat(data{ian}(ntind).mudata(days(idy)));
                    [sx, sy] = find(day_mu_psth');
                    f1 = scatter(sx,sy,2, 'filled');
                    
                    f1.MarkerFaceAlpha = 0.1;
                    
                    f1.MarkerFaceColor = [0 0 0];
                    axis tight
                    set(gca,'XTick',[]);
                    set(gca,'YTick',[]);
                    if idy == 1
                        ylabel(sprintf('%d',nt))
                    end
                    if int == numntrodes
                        xlabel(sprintf('%d',days(idy)))
                        set(gca,'XTick', [1 size(day_mu_psth,2)-1]);
                        if idy == 1
                            set(gca,'XTickLabel', round(data{ian}(ntind).mudata(1).time([1, end]),1));
                        end
                    end
                    % ripple line
                    centerline = round(size(day_mu_psth,2)/2);
                    lh = line([centerline centerline], [1 size(day_mu_psth,1)]);
                    lh.LineWidth = 1;
                    lh.Color = [.8 0 0 .8];
                    % epoch lines
                    le = line([1 size(day_mu_psth,2)], ...
                        repmat(data{ian}(ntind).eplengths{days(idy)}(1),2,1), 'Color', [0 0 0]);
                end
            end
        end
        
        %% ---- super title and colorbar----
        sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
        %         if justSU
        sprtit = sprintf('%s %s riptrigmu_unrefed', env, animal); % just change the stitle manually for diff conditions
        %         else
        %             sprtit = sprintf('%s %s riptrigmusu', env, animal);
        %         end
        iStitle = text(.5, .99, {sprtit}, 'Parent', sprtitleax, 'Units', ...
            'normalized');
        set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
            'horizontalAlignment', 'center','FontSize', 16);
        
        spylabel = text(.01, .5, sprintf('ntrode'), 'Parent', sprtitleax, 'Units', ...
            'normalized', 'Rotation', 90);
        set(spylabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
            'horizontalAlignment', 'center', 'FontSize', 16);
        
        spxlabel = text(.5, .01, sprintf('day'), 'Parent', sprtitleax, 'Units', ...
            'normalized');
        set(spxlabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
            'horizontalAlignment', 'center', 'FontSize', 16);
        %% ---- pause, save figs ----
        if pausefigs
            pause
        end
        if savefigs
            save_figure(paths.figdirectory, paths.filenamesave, sprtit)
            close all
        end
    end
end













