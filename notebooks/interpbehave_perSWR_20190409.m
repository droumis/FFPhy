
% load riptrig data, cellinfo
me = animaldef('demetris');
animals = {'D10', 'D13', 'JZ1', 'JZ3'};
env = 'wtrack';
loaddata = 0;

get_perf_perrip = 1;

plotfigs = 0;
savefigs = 1;
pausefigs = 0;
plotperday = 0;


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

        %load eventcons
        ripcons{an} = loaddatastruct(andef{2},andef{3},'ca1rippleskons');
    end
end

if get_perf_perrip
    % for each an, day - from ripcorr get the rip times from any of the pairs on this day, all epochs
    % from bstate, get an day's statespace allbound time, mode, and diffmode (need to add this col)
    % for each nt pair in ripcorr, get rip times 
    % interp1 with (allbound.time, allbound.mode(and seperately -
    % diffmode), rip.times) to get interped performance at rip times
    % concat the rip pair corr val with the interped perf val and run
    % corrcoef
    paircorr_wPerform = paircorr;
    for an = 1:length(animals)
        an_pairs = unique(data_indices{an}(:,[3 4]),'rows');
        days = unique(data_indices{an}(:,1));
        for d = 1:length(days)
            day = days(d);
            epochs = unique(data_indices{an}(data_indices{an}(:,1)==days(d),2));
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
                ripstarttimes = times([0 diff(includetimes')]==1);
                allbound = bstate{an}.BehaveState.statespace.allbound;
                allbound_inep_inds = ismember(allbound(:,[5 6]), [d epoch], 'rows');
                bstimes = bstate{an}.BehaveState.statespace.allepsMat(allbound_inep_inds,3);
                bsmode = allbound(allbound_inep_inds,1);
                bsmodediff = diff(allbound(allbound_inep_inds,1));
                % why are a bunch of the below results NaN? i think they
                % are rips outside the first and last trial time
                rip_bs = interp1(bstimes, bsmode, ripstarttimes);
                rip_bsdiff = interp1(bstimes(2:end), bsmodediff, ripstarttimes);
                for p = 1:length(an_pairs(:,1))
                    ntA = an_pairs(p,1);
                    ntB = an_pairs(p,2);
                    dayeppair_dataind = find(ismember(data_indices{an}, ...
                        [day epoch ntA ntB], 'rows'));
%                     xspcorr = paircorr{an}(dayeppair_dataind).excesscorr;
                    if ~isempty(dayeppair_dataind)
                    paircorr_wPerform{an}(dayeppair_dataind).performance = rip_bs;
                    paircorr_wPerform{an}(dayeppair_dataind).performance_diff = rip_bsdiff;
                    paircorr_wPerform{an}(dayeppair_dataind).ripstarttimes = ripstarttimes;
                    else
                        fprintf('no paircorr data found for %s day %d epoch %d ntA %d ntB %d...skipping\n', animals{an},day,epoch,ntA,ntB);
                        continue
                    end
                end
            end
        end
    end
end

if plotfigs
    
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

