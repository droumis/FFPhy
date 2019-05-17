
% the data to run this is made by dfa_perripspikingcorr
% as called in notebooks/perSWR_corr_20190507.m

% this notebook gathers the paircorr data per nt pair, across days, finds
% performance state for each rip time, and computes corr 


% load riptrig data, cellinfo
me = animaldef('demetris');
animals = {'D12'};
env = 'wtrack';
loaddata = 1;

get_perf_perrip = 0;
paircorr_perform_corr = 0;
save_results = 0;

plotfigs = 0;
plotperpair = 0;
savefigs = 0;
pausefigs = 0;
plotperday = 0;

paths = make_paths('paircorrVperformance', 'wtrack');

if loaddata
    for an =1:length(animals);
        % load behavestate
        andef = animaldef(animals{an});
        bstate{an} = load(sprintf('%s/%sBehaveState.mat',andef{2},animals{an}));
        % load riptrig_paircorr
        f = 'dfa_perripspikingcorr';
        tmp = load(sprintf('%s/%s/%s_%s_%s.mat',me{2}, f, f, env, animals{an}));
        paircorr{an} = tmp.F.output{1};
        % get [day ep tet clA clB] indices for the data output struct array
        data_indices{an} = cell2mat({paircorr{an}.index}');
        % filter for only wtrack days (no rotated wtrack)
        days = unique(data_indices{an}(:,1));
        task = loaddatastruct(andef{2},andef{3},'task');
        for d = 1:length(days)
            isrotated = cellfun(@(x) strcmp(x.environment,'wtrackrotated'), ...
                task{days(d)},'un',1);
            if any(isrotated)
                days(d) = -1;
            end
        end
        days = days(days>0);
        usedatainds = ismember(data_indices{an}(:,1),days);
        paircorr{an} = paircorr{an}(usedatainds);
        data_indices{an} = data_indices{an}(usedatainds,:);

        %load eventcons
        ripcons{an} = loaddatastruct(andef{2},andef{3},'ca1rippleskons');
    end
end

% i can definately move this part to the analysis function, since it's per
% epoch
if get_perf_perrip
    % for each an, day - from ripcorr get the rip times from any of the pairs on this day, all epochs
    % from bstate, get an day's statespace allbound time, mode, and diffmode (need to add this col)
    % for each nt pair in ripcorr, get rip times
    % interp1 with (allbound.time, allbound.mode(and seperately -
    % diffmode), rip.times) to get interped performance at rip times
    % concat the rip pair corr val with the interped perf val and run
    % corrcoef
    pcorrbs = paircorr;
    for an = 1:length(animals)
        an_pairs = unique(data_indices{an}(:,[3 4]),'rows');
        days = unique(data_indices{an}(:,1));
        andef = animaldef(animals{an});
        for d = 1:length(days)
            day = days(d);
            epochs = unique(data_indices{an}(data_indices{an}(:,1)==day,2));
            for e = 1:length(epochs)
                epoch = epochs(e);
                iep_datainds = find(ismember(data_indices{an}(:,[1 2]), [d epoch], 'rows'));
                times = ripcons{an}{d}{epoch}{1}.eegtimesvec_ref ;
                
                % the excludeperiods, based on ca1 ripples, are consistent
                % for all pairs for a given day epoch, so just use first
                excludeperiods = paircorr{an}(iep_datainds(1)).excludeperiods;
                % need to reconstruct the eventtimes from the exclude periods
                % i could/should have saved the eventimes into the result,
                % but at least this way is a good validation of getting
                % rip windows from excludetimes, as was done in the
                % analysis function
                includetimes = ~isExcluded(times, excludeperiods);
                includetimes = includetimes(:);
                
                %% process events
                % reconstitute the cons events
                % use this instead of the eventcons directly because it's already been
                % filtered to certain ripples according to the time filter options
                eventstart = times([0 diff(includetimes')]==1);
                eventend = times(diff(includetimes')==-1);
                
                % throw out last event if interrupted by end of epoch
                if (eventend(end)-eventstart(end) < 0)
                    eventstart = eventstart(1:(end-1));
                end
                % throw out first event if occurred before beginning of epoch
                if (eventend(1)-eventstart(1) < 0)
                    eventend = eventend(2:end);
                end
                % throw out any event that begins less than window(2) (i.e. 0.5 seconds) from end of epoch
                win = paircorr{an}(iep_datainds(1)).win;
                while (eventend(end) > times(end)-win(2))
                    eventend(end) = [];
                    eventstart(end) = [];
                end
                
                % choose event time to be the START of the event period
                ripstarttimes = eventstart';
                
                allbound = bstate{an}.BehaveState.statespace.allbound;
                allbound_inep_inds = ismember(allbound(:,[5 6]), [d epoch], 'rows');
                bstimes = bstate{an}.BehaveState.statespace.allepsMat(allbound_inep_inds,3);
                if isempty(bstimes)
                    fprintf('no behavestate data data found for %s day %d epoch %d...skipping\n', animals{an},day,epoch);
                    continue
                end
                bsmode = allbound(allbound_inep_inds,1);
                bsmodediff = abs(diff(allbound(allbound_inep_inds,1)));
                % why are a bunch of the below results NaN? i think they
                % are rips outside the first and last trial time
                rip_bs = interp1(bstimes, bsmode, ripstarttimes);
                rip_bsdiff = interp1(bstimes(2:end), bsmodediff, ripstarttimes);
                if length(rip_bs) ~= length(rip_bsdiff)
                    pause
                end
                for p = 1:length(an_pairs(:,1))
                    ntA = an_pairs(p,1);
                    ntB = an_pairs(p,2);
                    dep_ind = find(ismember(data_indices{an}, ...
                        [day epoch ntA ntB], 'rows'));
                    if ~isempty(dep_ind)
                        pcorrbs{an}(dep_ind).performance = rip_bs;
                        pcorrbs{an}(dep_ind).performance_diff = rip_bsdiff;
                        pcorrbs{an}(dep_ind).ripstarttimes = ripstarttimes;
                    else
                        fprintf('no paircorr data found for %s day %d epoch %d ntA %d ntB %d...skipping\n', animals{an},day,epoch,ntA,ntB);
                        continue
                    end
                end
            end
        end
    end
end

if paircorr_perform_corr
    for an = 1:length(animals)
        an_pairs = unique(data_indices{an}(:,[3 4]),'rows');
        
        for p = 1:length(an_pairs(:,1))
            ntA=an_pairs(p,1);
            ntB=an_pairs(p,2);
            ipair_inds = find(ismember(data_indices{an}(:,[3 4]),[ntA ntB], 'rows'));
%             ipair_inds = ipair_inds
            ipair_data = pcorrbs{an}(ipair_inds);
            %             inds_mat = cell2mat({ipair_data.index}');
            % vertcat the results across day eps for this pair
            try
                perf_mat = [ipair_data.performance]';
                perfdiff_mat = [ipair_data.performance_diff]';
            catch % if size M x 1 instead of 1 x M
                perf_mat = cell2mat({ipair_data.performance}');
                perfdiff_mat = cell2mat({ipair_data.performance_diff}');
            end
            try
                xc = [ipair_data.ntAB_excesscorr]';
            catch
                xc = cell2mat({ipair_data.ntAB_excesscorr}');
            end
            
            % get corr of ntAB_excesscorr VS bs,bsdiff across
            % day eps (which is why is can't be in the above loop)
            perf_mask = ~isnan(perf_mat);
            [R,P] = corrcoef(perf_mat(perf_mask), ...
                xc(perf_mask));
            perfdiff_mask = ~isnan(perfdiff_mat);
            [dR, dP] = corrcoef(perfdiff_mat(perfdiff_mask), ...
                xc(perfdiff_mask));

            % output results
            Fout(an).animal = animaldef(animals{an});
            Fout(an).output(p).pair_data = ipair_data;
            Fout(an).output(p).perf = perf_mat;
%             Fout(an).output(p).perf_mask = perf_mask;
%             Fout(an).output(p).perfdiff_mask = perfdiff_mask;
            Fout(an).output(p).perfdiff = perfdiff_mat;
            Fout(an).output(p).xc = xc;
            Fout(an).output(p).ntA = ntA;
            Fout(an).output(p).ntB = ntB;
            Fout(an).output(p).xc_perf_r = R(2,1);
            Fout(an).output(p).xc_perf_p = P(2,1);
            Fout(an).output(p).xc_diff_r = dR(2,1);
            Fout(an).output(p).xc_diff_p = dP(2,1);
        end
        if save_results
            save_filter_output(Fout, paths.filtOutputDirectory, paths.filenamesave)
        end
    end
end
%% continue with the plotting stuff below in a new script so i can load the 
% result and also the rasters? what if i saved the rasters along with each
% rip time so i wouldn't have to load more data files.. 
% i do def want to load the behavestate again in the new plotting script so
% that i can plot all the behave state along with the riptime-behavestate
% that is in the result.. or do i??
% i think all i really want for each pair is the day epoch riptimes, ...
% riptrig spike rasters, riptrig pair excess corr score, riptime perf/diff, ...
% and the P/R value between XC and perf/diff
% from that data i can make the psth per day/event, boxplot xc per day, scatter
% xc VS perf/diff with a labeled score.
% the important parts are to visually confirm that the riptrig spiking on
% the change day looks more coherent for the pair. 
% so now save the ntAB_results with the spike rasters included.

alpha = .001;

if plotfigs
    Pp = load_plotting_params('riptrigspiking');
    if plotperpair
        for an = 1:length(animals)
            pairs = [[ntAB_results{an}.ntA]' [ntAB_results{an}.ntB]'];
            pz = cell2mat({ntAB_results{an}.xc_diff_p}');
            pzsig = find(pz<alpha);
            ppairs = pairs(pzsig,:);
            for p = 1:length(ppairs)
                if savefigs && ~pausefigs;
                    close all
                    ifig = figure('Visible','off','units','normalized','position', ...
                        Pp.position);
                else
                    clf
                    ifig = figure('units','normalized','position',Pp.position);
                end
                set(gcf,'color','white')
                ppdataind = find(ismember(pairs, ppairs(p,:),'rows'));
                perfdiff = ntAB_results{an}(ppdataind).performance_change;
                perf = ntAB_results{an}(ppdataind).performance;
                xc = ntAB_results{an}(ppdataind).ntAB_excess_corr;
                R = ntAB_results{an}(ppdataind).xc_perf_r;
                P = ntAB_results{an}(ppdataind).xc_perf_p;
                dR = ntAB_results{an}(ppdataind).xc_diff_r;
                dP = ntAB_results{an}(ppdataind).xc_diff_p;
                
                clf
                subplot(1,2,1)
                scatter(xc,perf, 4, 1:length(xc))
                xlabel('performance')
                ylabel('excesscorr')
                title(sprintf('%03d %03d',R, P));
                axis tight
                lsline
                subplot(1,2,2)
                scatter(xc,perfdiff, 4, 1:length(xc))
                xlabel('perfdiff')
                ylabel('excesscorr')
                title(sprintf('%03d %03d',dR, dP));
                axis tight
                lsline
                %% ---- super title and colorbar----
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('riptrig MU-SU %s %s - nt%d', env, animals{an}, tet_ids(1,3));
                iStitle = text(.5, .95, {sprtit}, 'Parent', sprtitleax, 'Units', ...
                    'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center');
                Sylabel1 = text(.06, .7, 'swr #', 'Parent', sprtitleax, 'Units', ...
                    'normalized');
                set(Sylabel1,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center', 'Rotation',90);
                Sylabel2 = text(.06, .3, 'spike probability', 'Parent', sprtitleax, 'Units', ...
                    'normalized', 'Rotation',90);
                set(Sylabel2,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center');
                
                %% ---- pause, save figs ----
                if pausefigs
                    pause
                end
                if savefigs
                    save_figure(paths.figdirectory, paths.filenamesave, sprtit)
                end
            end
        end
    end
    
    
    if plotperday
        for an = 1:length(animals)
            indaybs = [];
            indaybsdiff = [];
            bsallb = bstate{an}.BehaveState.statespace.allbound;
            days = unique(bsallb(:,5));
            for d = 1:length(days);
                indaybs{d} = bsallb(bsallb(:,5)==days(d),[1 5]);
                indaybsdiff{d} = [diff(indaybs{d}(:,1)) indaybs{d}(1:end-1,2)];
            end
            figure(an)
            alldiff = cell2mat(indaybsdiff');
            allbs = cell2mat(indaybs');
            
            subplot(2,1,1)
            boxplot(allbs(:,1), allbs(:,2))
            ylabel('est prob correct')
            title(sprintf('%s',animals{an}))
            
            subplot(2,1,2)
            boxplot(alldiff(:,1), alldiff(:,2))
            ylabel('diff')
        end
    end
    
end