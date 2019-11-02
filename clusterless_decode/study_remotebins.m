datadir_superlin = '/opt/data13/kkay/Superlin_data/';
savedir = '/opt/data13/kkay/__WG';

runcollect = 1;
if runcollect
    plot_probe = 1;
        plot_epoch = 1;
            plotwin_length = 60;        
                 tracekern_ms = [10];           
    binsize = 0.3;   % in sec
    min_remote_units = 2;  % minimum remote cells active 
    ratethresh_field = 2 ;  % in Hz, place field (epoch-active) threshold
                            % 0.5, 1, 2, 5 Hz
    epochstring = 'runW';
    ripple_removal = 1;
    % load active indicator matrix
    load('/opt/data13/kkay/__WG/ACTIVEMATRIX.mat','ACT')
    load('/opt/data13/kkay/__WG/ACTIVEMATRIX_REST.mat','ACT_R')
    load('/opt/data13/kkay/__WG/ACTIVEMATRIX_BOTH.mat','ACT_B')
        threshind = find(ratethresh_field == ACT.ratethresh_field);
    animals_torun = {'Bond'} ; {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
    manual_epochs = [4 6];
        remote_epinds = [1 2]; % for rest box case: [1] if just latest RUN epoch in day, [2] if second latest RUN epoch of different type
                               %  [1 2] if require units to be clustered
                               %  for both W-track RUN epochs
end

if runcollect
  
    for aa = 1:length(animals_torun)
        
        animalname = animals_torun{aa};
        an = find(strcmp(animalname,ACT.animals_order));
        animalinfo = animaldef(animalname);
    
        disp(animalname)

        % Load data
        task = loaddatastruct(animalinfo{2},animalinfo{3},'task');
        pos = loaddatastruct(animalinfo{2},animalinfo{3},'pos');    
        spikes = loaddatastruct(animalinfo{2},animalinfo{3},'spikes');  
        
        % Retrieve epochs here
        if isempty(manual_epochs)
            epochs = evaluatefilter(task,epochmaker(epochstring));
        else
            epochs = manual_epochs;
        end
        epochs_w =  evaluatefilter(task,epochmaker('runW'));       
        numeps = size(epochs,1);
        
        clear R
        R.animalname = animalname;
        R.epochs = epochs;
        R.binsize = binsize;
        R.min_remote_units = min_remote_units;
        R.epochstring = epochstring;
        R.ratethresh_field = ratethresh_field;
        R.remotebins = [];
        
        for ee = 1:numeps
            
            d = epochs(ee,1);
            ep = epochs(ee,2);
            
            epstart = pos{d}{ep}.data(1,1);
            epend = pos{d}{ep}.data(end,1);
                postimevec = pos{d}{ep}.data(:,1);
            
            binvec = epstart:binsize:epend;             % bin edges
            binvec_c = binvec(1:(end-1)) + binsize/2 ;  % bin centers
       
            numbins = length(binvec) - 1 ;
            
            % Identify remote epoch for this current epoch  (latest in the day)
            day_eps = epochs_w(epochs_w(:,1)==d,2)';
            if ~isfield(task{d}{ep},'type')
                continue
            else
                eptype = task{d}{ep}.type;
            end
            
            ep_r = [];
            if strcmp(eptype,'run')
                runsleep = 1;
                currenv = task{d}{ep}.environment;
                remoteeps = [];
                for xx = day_eps
                    if ~strcmp(currenv,task{d}{xx}.environment)
                        remoteeps = [remoteeps xx];
                    end
                end       
                ep_r = max(remoteeps);  %take the latest epoch (most efficient running behavior)
            elseif strcmp(eptype,'sleep')
                % if sleep, find the two latest different run epochs
                latestrun_ep = nan;
                for xx = sort(day_eps,'descend')
                    if isfield(task{d}{xx},'type') && strcmp(task{d}{xx}.type,'run')
                        latestrun_env = task{d}{xx}.environment;
                        ep_r = [ep_r  xx];
                        break
                    end
                end      
                for yy = sort(day_eps,'descend')
                    if isfield(task{d}{yy},'type') && strcmp(task{d}{yy}.type,'run')
                        if ~strcmp(task{d}{yy}.environment,latestrun_env)
                            latestrun_ep = task{d}{yy}.environment;
                            ep_r = [ep_r  yy];
                            break
                        end
                    end
                end
                % specify which RUN epochs you want to include with "remote_epinds"
                ep_r = ep_r(remote_epinds);    
            else
                keyboard
            end
            
            if isempty(ep_r)
                disp([num2str([d ep]) ' has no remote ep'])
                continue
            end
            
           numrems = length(ep_r);
            
            % Load SWR data
            timefilterscript
            ripout = kk_getconstimes(animalinfo{2},animalinfo{3},[d ep ],'ripplescons',1,...
                'consensus_numtets',3,'minthresh',2,...
                'exclusion_dur',0,'minvelocity',0,'maxvelocity',4);
            consvec_rip2 = ripout{d}{ep}.cons;
            consvectimes_rip2 = ripout{d}{ep}.time;
            ripperiods = vec2list(consvec_rip2,consvectimes_rip2);
           
            % Identify units clustered in these epochs
            adtca = [];   % "a" (5th column)  -- ctive flag value
                          % 1: local only , 2: remote1 only, 3: remote2 only, 4: remote1 & remote2 only, 5: something active, 6: nothing active
                for ccc = 1:size(ACT_B.adtc,1)
                    adtc = ACT_B.adtc(ccc,:);
                    if all( [an d] == [ACT_B.adtc(ccc,1:2)]  )  % only look at units clustered on same day as this epoch
                        actmat = ACT_B.actmatrix{ccc};
                        if all ( ismember([ep ep_r], actmat(:,1)) )  % check that the unit is clustered in both local + remote epochs
                            % local
                            rownum_loc = find(ep == actmat(:,1));
                            loc = actmat(rownum_loc,threshind);  % 0: inactive, 1: active
                            % remote(s)
                            clear rem
                            for mm = 1:numrems
                                rownum_rem = find(ep_r(mm) == actmat(:,1));
                                rem(mm) = actmat(rownum_rem,threshind);  % 0: inactive, 1: active
                            end
                            mat = [loc rem];
                                single = (sum(mat) == 1);
                            if loc && single
                                active = 1;
                            elseif rem(1) && single
                                active = 2;
                            elseif numrems == 2 && rem(2) && ~rem(1)
                                active = 3;
                            elseif numrems == 2 && all(rem == 1) && ~loc
                                active = 4;
                            elseif sum(mat) > 0
                                active = 5;
                            elseif sum(mat) == 0
                                active = 6;
                            else
                                keyboard
                            end
                            adtca = [adtca ; adtc active];
                        end
                    end
                end

            
            % Sort unit list by active type
            [~,I] = sort(adtca(:,5));
            adtca = adtca(I,:);
            
            numunits = size(adtca,1);
           
            % Go through each cell and determine spike indicator matrix (fired at least one spike) for each bin
            spikemat = zeros(numunits,numbins);
            for c = 1:numunits
                tet = adtca(c,3);
                cellnum = adtca(c,4);
                if ~isempty(spikes{d}{ep}{tet}{cellnum}.data)
                    spiketimes = spikes{d}{ep}{tet}{cellnum}.data(:,1);
                else
                    spiketimes = [];
                end
                % Assign each spike into a bin
                if ~isempty(spiketimes)
                    [spikeobs,~] = histc(spiketimes, binvec);
                        spikeobs = spikeobs(1:(end-1));
                    spikemat(c,:) = (spikeobs > 1);
                end
            end
            
            % Identify candidate remote bins
            rows_loc = (adtca(:,5) == 1);
            rows_rem = (adtca(:,5) == 2) | (adtca(:,5) == 3) | (adtca(:,5) == 4);        
            numloc = sum(rows_loc);
            numrem = sum(rows_rem);
            remsum = sum(spikemat(rows_rem,:),1);
            locsum = sum(spikemat(rows_loc,:),1);
            minbins = find(remsum >= min_remote_units); %   % index of bins that have at least min # of remote units firing
                rem_bin = remsum(minbins);
                loc_bin = locsum(minbins);
            if ~isempty(minbins)
                disp([num2str([d ep]) ' : ' num2str(rem_bin)])
            else
                continue
            end

            % Remove remote bins that contain or overlap with SWRs
            if ripple_removal
                for bb = length(minbins):-1:1
                    bind = minbins(bb);
                    bin_start = binvec(bind);
                    bin_end = binvec(bind+1);
                    ripoverlap =  logical(  isExcluded(bin_start,ripperiods) || isExcluded(bin_end,ripperiods) || ...
                        any(isExcluded(ripperiods(:),[bin_start bin_end]))   ) ;
                    if ripoverlap
                        minbins(bb) = [];
                        rem_bin(bb) = [];
                        loc_bin(bb) = [];
                        disp('rip removal')
                    end
                end
            end
            
            % Calculate binomial probability for each bin
            if 0
                numrembins = length(minbins);
                pvalue = nan(1,numrembins);
                for rr = 1:numrembins
                    num_firingrem = rem_bin(rr);
                    num_firingloc = loc_bin(rr);
                    % treat these two scenarios differently
                    p_binomial = numrem / (numrem + numloc);
                    pvalue(rr) = binocdf(num_firingrem,numloc + numrem,p_binomial);
                end
            end
            
            if plot_probe
            
                plotwinsize = 5;  % seconds
                
                % Ripple and WG power traces
                out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'rippletrace', d, ep);
                riptrace = zscoretrace(out{d}{ep}{1},[]);
                riptrace_timevec = out{d}{ep}{1}.eegtimesvec_ref;
                out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'fastgammatrace', d, ep);
                wgtrace = zscoretrace(out{d}{ep}{2},tracekern_ms);
                wgtrace_timevec = out{d}{ep}{2}.eegtimesvec_ref;
                out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'lowgammatrace', d, ep);
                lgtrace_ca1 = zscoretrace(out{d}{ep}{1},tracekern_ms);
                out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'highgammatrace', d, ep);
                hgtrace_ca1 = zscoretrace(out{d}{ep}{1},tracekern_ms);
                lgtrace_ca3 = zscoretrace(out{d}{ep}{3},tracekern_ms);
                lgtrace_timevec = out{d}{ep}{1}.eegtimesvec_ref;
                if 0
                    vte = loaddatastruct(animalinfo{2}, animalinfo{3}, 'vte',d);
                    vte_start = vte{d}{ep}.starttime;
                    vte_end = vte{d}{ep}.endtime;
                end

                %if ~isempty(tracekern_ms)
                %    sig_bins = 1500 * tracekern_ms/1000;
                %    gauskern = gaussian(sig_bins,8*sig_bins);
                %    wgtrace = smoothvect(wgtrace,gauskern);
                %    lgtrace_ca1 = smoothvect(lgtrace_ca1,gauskern);
                %    lgtrace_ca3 = smoothvect(lgtrace_ca3,gauskern);
                %end
            
                % scroll plot

                numwins = floor(binvec_c(end) / plotwin_length);
                for ww = 1:numwins
                
                    figure('units','normalized','outerposition',[0 0 1 .5])
                    
                    bintime = binvec_c(ww);

                        plotstart = (ww-1)*plotwin_length + binvec_c(1) ;
                        plotend = plotstart + plotwin_length;

                    
                    % First plot ripple patches and wg traces
                    ax(1) = subplot(5, 2, [1 3 5 7]);
                    numunits;
                    ylevel = numunits+20;
                    % Unit spike raster
                    
                    % first plot active type backgrounds
                    for cc = 1:numunits
                        tet = adtca(cc,3);
                        cellnum = adtca(cc,4);
                        active = adtca(cc,5);
                        % 1: local only , 2: remote1 only, 3: remote2 only, 4: remote1 & remote2 only, 5: something active, 6: nothing active
                        if active == 1
                            clr = [.95 .95 .95];
                        elseif active == 2
                            clr = [.8 .8 1];
                        elseif active == 3
                            clr = [1 .8 .8];
                        elseif active == 4
                            clr = [.8 .6 .8];
                        elseif active == 5
                            clr = [1 1 1];
                        elseif active == 6
                            clr = [.6 .6 .6];
                        end
                        patch([plotstart plotstart plotend plotend]- epstart,[cc-0.5 cc+0.5 cc+0.5 cc-0.5],clr,'edgecolor','none');
                        hold on
                    end
                    
                    % next plot ripples and traces
                    % plot 2 SD SWRs that occur in the widnow %%%%%%%%%
                    plotwinvec = logical(list2vec([plotstart plotend],ripout{d}{ep}.time))';
                    consvec_rip2_toplot = consvec_rip2 & plotwinvec;
                    rip_toplot = vec2list(consvec_rip2_toplot,ripout{d}{ep}.time);
                    numrip_toplot = size(rip_toplot,1);
                    % plot all 2 SD ripples in the raster window
                    for rp = 1:numrip_toplot
                        patch([rip_toplot(rp,1) rip_toplot(rp,2) rip_toplot(rp,2) rip_toplot(rp,1)] - epstart,...
                            [0 0 ylevel ylevel],...
                            [1 .89 .89],'edgecolor','none'); hold on
                    end
                    
                    for cc = 1:numunits
                        tet = adtca(cc,3);
                        cellnum = adtca(cc,4);
                        % plot spikes
                        if ~isempty(spikes{d}{ep}{tet}{cellnum}.data)
                            spiketimes = spikes{d}{ep}{tet}{cellnum}.data(:,1);
                            spiketimes = spiketimes(logical(isExcluded(spiketimes,[plotstart plotend]))) ;  % spikes within window
                            numunitspk = length(spiketimes);
                            if numunitspk > 0
                                for ll = 1:numunitspk
                                    plot([spiketimes(ll) spiketimes(ll)]- epstart,[cc-0.48 cc+0.48],...
                                        'linestyle','-','Color',[0 0 0],'LineWidth',1,'Parent',ax(1)); hold on
                                end
                            end
                        end
                    end
                 
                    % plot WG and SWR power traces at bottom%%%%%%%%%%%%%%%%
                    ind1 = lookup(plotstart,lgtrace_timevec) ;
                    ind2 = lookup(plotend,lgtrace_timevec) ;
                    ind3 = lookup(plotstart,lgtrace_timevec) ;
                    ind4 = lookup(plotend,lgtrace_timevec) ;
                    plot((lgtrace_timevec(ind1:ind2))- epstart, 4 * -riptrace(ind1:ind2)/2 + ylevel + 10,'-','Color',[.8 .75 .75],'linewidth',2,'Parent',ax(1) ); hold on
                    plot((lgtrace_timevec(ind3:ind4))- epstart, 5 * -wgtrace(ind3:ind4)/2 + ylevel + 10,'-','Color',[0 .3 .9],'linewidth',2,'Parent',ax(1)); hold on
                    plot((lgtrace_timevec(ind1:ind2))- epstart, 5 * -lgtrace_ca1(ind1:ind2)/2 + ylevel + 10,'-','Color',[1 .1 .1],'linewidth',2,'Parent',ax(1)); hold on
                    plot((lgtrace_timevec(ind1:ind2))- epstart, 5 * -hgtrace_ca1(ind1:ind2)/2 + ylevel + 10,'-','Color',[0 .7 0],'linewidth',2,'Parent',ax(1)); hold on
                    %plot((lgtrace_timevec(ind1:ind2) - epstart), 4 * -lgtrace_ca3(ind1:ind2)/2 + ylevel + 10,'-','Color',[1 .4 .4],'linewidth',2 ,'Parent',ax(2)); hold on
                    
                    
                    xlim([plotstart plotend]- epstart)
                    ylim([0  (ylevel + 15)])
                    set(gca,'ydir','reverse')
                    set(gca,'fontsize',12)
                    set(gca,'tickdir','out');
                    xlabel('Time (ms)') ;   ylabel('unit #');
                    
                    
                    
                    % Plot speed
                    ax(2) = subplot(5, 2, 9);
                    a = lookup(plotstart,postimevec);
                    b = lookup(plotend,postimevec);
                    plot(pos{d}{ep}.data(a:b,1) - epstart,pos{d}{ep}.data(a:b,5),'-k','linewidth',2,'Parent',ax(2))
                    ylim([0 20])
                    
                    linkaxes(ax,'x')
                    %break

                    pause
                    close all
                    
                end
                
                
                
                
            end
            
        end
       
        
        
    end
    
end