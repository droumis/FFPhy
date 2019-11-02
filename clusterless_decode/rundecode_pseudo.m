% Cluster decoding on W-track -- psuedo 1D "armdists" (Wu--Foster-2014 Fig 1)


%datadir = '/opt/data13/kkay/Superclustdecode_data';
datadir = '/opt/data13/kkay/Superclustdecode_data/100ms_decodes'

calculate = 0;
plot_continuous = 1;
    if plot_continuous 
       animal_toplot = 'Frank';
            animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
       dayep = [7 2];
       starttime = [0]  ;  % time (within epoch, not clock time) to begin plot
       endtime = [900] ;
       plot_remoteW = 0;
       tracekern_ms = [];       
       plot_gammadiff = 0;
       plot_onlyminspike = 0;
       plot_hg = 0;
       plot_delta = 0;
       plot_vte = 0;
       plot_riptrace = 1;
       plot_gammamap = 0;
           mindur_gamma = 1;
       running_corr = 0;
       omitraster = 0;
       omitposterior = 0;
       plot_nonlocalbins = 1;
    end
plot_events = 0;
    if plot_events
       animal_toplot = 'Bond';
       eventname = 'ripples';
       dayep = [3 2];
       loaddata = 1;
    end
if calculate
    %%% select data %%
    animals_tocalc = {'Bond'}; %{'Frank','Government','Dave','Corriander','Egypt'}; {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'}; %{'Government','Corriander','Dave','Egypt','Chapati'}; %{'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
    cellchoice =  1;  % look inside function to see what is specified
    dayeps = [9 2; 9 4; 9 6; 10 2; 10 4; 10 6];  %[6 2; 6 4; 6 6];  % leave empty if want to do all epochs
    manual_period = [];
    winsize = .1 %.015;  % in seconds, decoding bin size
    remoteW = 1;   % set to 1 if want to decode in other W as well
    exclude_inters = 1;
    plot_infunction = 0;
        selected_events = [];
    %%% select what to decode %%
    decodemode = 1; % 1: entire epoch, 2: SWRs
    extratime = 500;  % if not doing decodemode 1, the ms around event start to plot
    %%% select decoding parameters %%
    modelnum = 3;   % 1: all spikes, 2: trajencode, 3: exclude SWR only   4: 10 cm/s and above
end

if calculate
    for aa = 1:length(animals_tocalc)
        animalname = animals_tocalc{aa};
            animalinfo = animaldef(animalname);
        kk_clusterdecode(animalinfo{2},animalinfo{3},dayeps,animalname,...
            'savedir',datadir,'decodemode',decodemode,'cellchoice',cellchoice,'modelnum',modelnum,...
            'plot_powertraces',0,'calctraj',1,...
            'extratime',extratime,'remoteW',remoteW,'plot_infunction',plot_infunction,...
            'winsize',winsize,'manual_period',manual_period,'exclude_inters',exclude_inters);
    end
end

if plot_continuous
    
    animpref = animal_toplot(1:3);
        animalinfo = animaldef(animal_toplot);
        daydir = getdaydir(animal_toplot);   
    d = dayep(1);
    ep = dayep(2);


    % load data    
    if ~exist('P','var') || ~strcmp(P.animalname,animal_toplot) || ~all(P.dayep == dayep)
        
        % posteriors
        cd(datadir)
        load(sprintf('%s_%d_%d.mat',animpref,d,ep),'P');
        
        % position
        pos = loaddatastruct(animalinfo{2},animalinfo{3},'pos',d);
        postimevec = pos{d}{ep}.data(:,1);
        if size(pos{d}{ep}.data,2) == 5
            xtrace = pos{d}{ep}.data(:,2);
            ytrace = pos{d}{ep}.data(:,3);
            veltrace = pos{d}{ep}.data(:,5);
            dirtrace = pos{d}{ep}.data(:,4);
        else
            xtrace = pos{d}{ep}.data(:,6);
            ytrace = pos{d}{ep}.data(:,7);
            veltrace = pos{d}{ep}.data(:,9);
            dirtrace = pos{d}{ep}.data(:,4);       % use old head dir      
        end
        epstart = postimevec(1);
        epend = postimevec(end);
        
        % clustered spikes
        spikes = loaddatastruct(animalinfo{2},animalinfo{3},'spikes',d);
        cellinfo = loaddatastruct(animalinfo{2},animalinfo{3},'cellinfo');
        
        % ripples, 2 SD
        ripout = kk_getconstimes(animalinfo{2},animalinfo{3},[d ep ],'ripplescons',1,...
            'consensus_numtets',3,'minthresh',2,...
            'exclusion_dur',0,'minvelocity',0,'maxvelocity',4);
        consvec_rip2 = ripout{d}{ep}.cons;
        consvectimes_rip2 = ripout{d}{ep}.time;
        periodtimes_rip2 = vec2list(consvec_rip2,consvectimes_rip2);
        
        % Ripple and WG power traces
        out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'rippletrace', d, ep);
        [riptrace,meanz,stdz] = zscoretrace(out{d}{ep}{1},[]);
        riptrace_timevec = out{d}{ep}{1}.eegtimesvec_ref;
        out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'fastgammatrace', d, ep);
        wgtrace = zscoretrace(out{d}{ep}{2},tracekern_ms);
        wgtrace_timevec = out{d}{ep}{2}.eegtimesvec_ref;
        out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'lowgammatrace', d, ep);
        lgtrace_ca1 = zscoretrace(out{d}{ep}{1},tracekern_ms);
        lgtrace_ca3 = zscoretrace(out{d}{ep}{3},tracekern_ms);
        lgtrace_timevec = out{d}{ep}{1}.eegtimesvec_ref;  
        if plot_hg
            out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'highgammatrace', d, ep);
            hgtrace = zscoretrace(out{d}{ep}{1},tracekern_ms);
        end
        if plot_delta
            out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'deltatrace', d, ep);
            dtrace_timevec = out{d}{ep}{2}.eegtimesvec_ref;  
            dtrace = zscoretrace(out{d}{ep}{2}.powertrace,tracekern_ms);
        end
        if plot_vte
             vte = loaddatastruct(animalinfo{2}, animalinfo{3}, 'vte',d);
                vte_start = vte{d}{ep}.starttime;
                vte_end = vte{d}{ep}.endtime;
        end
      
    end
    
    if iscell(P.binvec_c)
        timevec = P.binvec_c{1} - P.binvec_c{1}(1);
    else
        timevec = P.binvec_c - P.binvec_c(1);
    end
        binsize = timevec(2) - timevec(1);
        timevec = timevec + binsize/2;
        if starttime == endtime
            a = lookup(P.replay{1}(starttime,1),timevec);
            b = lookup(P.replay{1}(starttime,2),timevec);
        else
                a = lookup(starttime,timevec);
                b = lookup(endtime,timevec);
        end
    armdists = P.armdists;
    linposbins = P.linposbins{1};
    centerarmmax = P.centerarmmax;
    rightarmmax = P.rightarmmax;
    detc = P.detc{1};
    
   figure('units','normalized','outerposition',[0 0 .8 .8])
    
    % First plot ripple patches and wg traces

    ymax = 20;
    
    % obtain all clustered units this epoch
    an = find(strcmp(animal_toplot,animals_order));
    %   [adtc_list] = clusteredunits(spikes,cellinfo,an,d,ep,[],1);
    
    maxcell = size(detc,1);
    ax(1) = subplot(1,3,1);
    ylevel = maxcell+20;

    % Unit spike raster
    if ~omitraster
    ax(1) = subplot(12,1,[1:5]);
    celllist = {};
    % iterate
    for cc = 1:maxcell
        tet = detc(cc,3);
        cellnum = detc(cc,4);
        region = cellinfo{d}{ep}{tet}{cellnum}.area;
            [~,~,~,bgclr] = regionfunc(region);
        type = cellinfo{d}{ep}{tet}{cellnum}.type;
        celllist = [celllist ; num2str(cc) type region];
        disp(type)
        if ~strcmp(type,'principal')
            bgclr = [1 1 1];
        end
        
        % plot region (CA1, CA3) color patch for region
        patch([starttime starttime endtime endtime],[cc-0.5 cc+0.5 cc+0.5 cc-0.5],bgclr,'edgecolor','none');
            hold on
       
    end
    end
    
    % Plot candidate bins (if in P structure already)
    if plot_nonlocalbins
        nonlocalclr = [.85 1 .9];
    else
        nonlocalclr = [1 1 1];
    end
    if plot_nonlocalbins
        ax(1) = subplot(12,1,[1:5]);
        candtimes = P.binvec_c{1}(P.candbins) - epstart;
        candtimes(  (candtimes < starttime) & (candtimes > endtime)  ) = [];  % ignore those outside plotting window
        for ll = candtimes
            time_1 = ll - binsize/2;
            time_2 = ll + binsize/2;
            patch([time_1 time_1 time_2 time_2],...
                  [0 (ylevel + ymax + 2) (ylevel + ymax + 2) 0],nonlocalclr,'edgecolor','none');
            hold on                 
        end
    end
    
     % Plot 2 SD SWRs that occur in the widnow %%%%%%%%%
    plotwinvec = logical(list2vec([starttime endtime] + epstart,ripout{d}{ep}.time))';
    consvec_rip2_toplot = consvec_rip2 & plotwinvec;
    rip_toplot = vec2list(consvec_rip2_toplot,ripout{d}{ep}.time) - epstart;
    numrip_toplot = size(rip_toplot,1);
    % plot all 2 SD ripples in the raster window
    for rp = 1:numrip_toplot
        patch([rip_toplot(rp,1) rip_toplot(rp,2) rip_toplot(rp,2) rip_toplot(rp,1)] ,...
            [0 0 (ylevel + ymax) (ylevel + ymax)],...
            [1 .84 .84],'edgecolor','none','Parent',ax(1) ); hold on
    end
    
    % plot spike raster
    if ~omitraster
        for cc = 1:maxcell
            tet = detc(cc,3);
            cellnum = detc(cc,4);
            region = cellinfo{d}{ep}{tet}{cellnum}.area;
            [~,~,~,bgclr] = regionfunc(region);
            type = cellinfo{d}{ep}{tet}{cellnum}.type;
            disp(type)
            if ~strcmp(type,'principal')
                bgclr = [1 1 1];
            end
            % plot spikes
            if ~isempty(spikes{d}{ep}{tet}{cellnum}.data)
                spiketimes = spikes{d}{ep}{tet}{cellnum}.data(:,1);
                spiketimes = spiketimes(logical(isExcluded(spiketimes,[starttime endtime]+epstart))) ;  % spikes within window
                numunitspk = length(spiketimes);
                if numunitspk > 0
                    for ll = 1:numunitspk
                        plot([spiketimes(ll) spiketimes(ll)] - epstart,[cc-0.48 cc+0.48],...
                            'linestyle','-','Color',[0 0 0],'LineWidth',1,'Parent',ax(1) ); hold on
                    end
                end
            end
        end
    else
        set(gca,'ytick',0)
    end
   
    % plot WG and SWR power traces at bottom%%%%%%%%%%%%%%%%
    ind1 = lookup(starttime,riptrace_timevec - epstart) ;
    ind2 = lookup(endtime,riptrace_timevec - epstart) ;
    ind3 = lookup(starttime,wgtrace_timevec - epstart) ;
    ind4 = lookup(endtime,wgtrace_timevec - epstart) ;
    if plot_delta
        ind5 = lookup(starttime,dtrace_timevec - epstart) ;
        ind6 = lookup(endtime,dtrace_timevec - epstart) ;
    end
    if ~isempty(tracekern_ms)
        GFACTOR = 8;
    else
        GFACTOR = 6;
    end
    
    if plot_riptrace
        RIP = 10;  % ripple factor
        plot((riptrace_timevec(ind1:ind2) - epstart), RIP * -riptrace(ind1:ind2)/2 + ylevel + 10,'-','Color',[.8 .75 .75],'linewidth',2,'Parent',ax(1) ); hold on
        % mean line
        plot([riptrace_timevec(ind1) riptrace_timevec(ind2)] - epstart, [(ylevel+10+meanz*RIP)  (ylevel+10+meanz*RIP)]  ,'--','Color',[.1 .1 .1],'linewidth',2,'Parent',ax(1) ); hold on
        % +/- 1 SD
        plot([riptrace_timevec(ind1) riptrace_timevec(ind2)] - epstart, [(ylevel+10+(meanz+stdz)*RIP)  (ylevel+10+(meanz+stdz)*RIP)]  ,'--','Color',[.4 .4 .4],'linewidth',2,'Parent',ax(1) ); hold on
        plot([riptrace_timevec(ind1) riptrace_timevec(ind2)] - epstart, [(ylevel+10+(meanz-stdz)*RIP)  (ylevel+10+(meanz-stdz)*RIP)]  ,'--','Color',[.4 .4 .4],'linewidth',2,'Parent',ax(1) ); hold on
    end
    
    plot((wgtrace_timevec(ind3:ind4) - epstart), GFACTOR * -wgtrace(ind3:ind4)/2 + ylevel + 10,'-','Color',[0 .3 .9],'linewidth',2,'Parent',ax(1) ); hold on
    plot((lgtrace_timevec(ind3:ind4) - epstart), GFACTOR * -lgtrace_ca1(ind3:ind4)/2 + ylevel + 10,'Color',[1 .1 .1],'linewidth',2,'Parent',ax(1) ); hold on
    %plot((lgtrace_timevec(ind3:ind4) - epstart), GFACTOR * -lgtrace_ca3(ind3:ind4)/2 + ylevel + 10,'Color',[1 .4 .4],'linewidth',2,'Parent',ax(1) ); hold on
    
    
    plot([lgtrace_timevec(ind3) lgtrace_timevec(ind4)]-epstart,[ ylevel + 10   ylevel + 10],'--','Color',[0 0 0],'linewidth',1,'Parent',ax(1));
    % +/- 1 SD
    plot([lgtrace_timevec(ind3) lgtrace_timevec(ind4)]-epstart,[ ylevel + 10 + GFACTOR   ylevel + 10 + GFACTOR],'--','Color',[.6 .6 .6],'linewidth',1,'Parent',ax(1));
    plot([lgtrace_timevec(ind3) lgtrace_timevec(ind4)]-epstart,[ ylevel + 10 - GFACTOR   ylevel + 10 - GFACTOR],'--','Color',[.6 .6 .6],'linewidth',1,'Parent',ax(1));
    
    if plot_hg
        plot((lgtrace_timevec(ind3:ind4) - epstart), GFACTOR * -hgtrace(ind3:ind4)/2 + ylevel + 10,'Color',[.5 1 .6],'linewidth',2,'Parent',ax(1) ); hold on
    end
    if plot_delta
        plot((dtrace_timevec(ind5:ind6) - epstart), 16 * -dtrace(ind5:ind6)/2 + ylevel + 10,'Color',[.2 .2 .2],'linewidth',2,'Parent',ax(1) ); hold on
    end
    if 1
        gammadiff = lgtrace_ca1(ind3:ind4) -  wgtrace(ind3:ind4);
            antitimes_1 = (lgtrace_ca1(ind3:ind4) > 0) & (wgtrace(ind3:ind4) < 0);
            antitimes_2 = (lgtrace_ca1(ind3:ind4) < 0) & (wgtrace(ind3:ind4) > 0);
            antitimes = antitimes_1 | antitimes_2;
                gammadiff2 = gammadiff;
                gammadiff2(~antitimes) = nan;
                
        % plot gamma power difference
        if plot_gammadiff
            plot([lgtrace_timevec(ind3)  lgtrace_timevec(ind4)] - epstart,...
                [ylevel ylevel] - 10,'--','Color',[.4 .4 .4],'linewidth',1,'Parent',ax(1) ); hold on            
            plot(lgtrace_timevec(ind3:ind4) - epstart, - 2 * gammadiff2 + ylevel - 10,'Color',[0 0 0],'linewidth',3,'Parent',ax(1) ); hold on
        end
        % plot running correlation
        
    end
    
    if starttime ~= endtime
        xlim([starttime endtime])
    else
    end
    ylim([0  (ylevel + ymax )])
    set(gca,'ydir','reverse')
    set(gca,'fontsize',12)
    set(gca,'tickdir','out');
    xlabel('Time (ms)') ;
    
    if ~omitraster
        ylabel('unit #');
    end
    
    if running_corr
        timevec_trace = lgtrace_timevec(ind3:ind4) - epstart;
        win =   0.5 ;
        numwins = floor(  (timevec_trace(end) - timevec_trace(1)) / win ) ;
        
        corrtrace = nan(1,numwins);
        timetrace = nan(1,numwins);
        for w = 1:numwins
            j = floor(   (w-1) * win * 1500 + 1   );
            k = floor(    w*win*1500              );
            m = round(mean([j k]));
                timetrace(w) = timevec_trace(m);
            [RHO,PVAL] = corr(lgtrace_ca1(j:k)',wgtrace(j:k)');
            if PVAL < 0.05
                corrtrace(w) = RHO;
            end
        end
        if 1
            %figure;
            plot(timetrace,-10 * corrtrace + ylevel + 40,'ok','linewidth',2); hold on
            plot(timetrace,-10 * corrtrace + ylevel + 40,'-k','linewidth',2); hold on
            plot([timetrace(1) timetrace(end)],[-10*1+ylevel+40   -10*1+ylevel+40],'k--','linewidth',1)
            plot([timetrace(1) timetrace(end)],[10*1+ylevel+40     10*1+ylevel+40],'k--','linewidth',1)
            disp(num2str(sum(isnan(corrtrace))))
        end
        %keyboard
    end
    
    % 2D map looking at where gamma diff occurs on track
    if plot_gammamap
        figure
        if 0
            ax(2) = subplot(12,1,[6:8]);
            plot(lgtrace_timevec(ind3:ind4) - epstart, gammadiff,...
                'Color',[0 0 0],'linewidth',3,'Parent',ax(2) ); hold on
            linkaxes(ax,'x');
        end
        if 1
            clear gperiods
            gperiods{1} = vec2list(antitimes_1,lgtrace_timevec(ind3:ind4));
            gperiods{2} = vec2list(antitimes_2,lgtrace_timevec(ind3:ind4));

            for mm = 1:2
                for p = size(gperiods{mm},1):-1:1
                    if (gperiods{mm}(p,2) - gperiods{mm}(p,1)) < mindur_gamma
                        gperiods{mm}(p,:) = [];
                    end
                end
            end
            
            % sort periods chronologically
            bothperiods = [gperiods{1}   1 * ones(size(gperiods{1},1),1)   ;   gperiods{2}  2 * ones(size(gperiods{2}),1)];
            bothperiods = sortrows(bothperiods,1);
            
            % plot 2D map
            plot(xtrace,ytrace,'.','Color',[.85 .85 .85]); hold on
            for bp = 1:size(bothperiods,1)
                mmm = bothperiods(bp,3);
                if mmm == 1
                    clr = [1 .1 .1];
                elseif mmm == 2
                    clr = [0 .3 .9];
                end
                g = lookup(bothperiods(bp,1),postimevec);
                f = lookup(bothperiods(bp,2),postimevec);
                plot(xtrace(g:f),ytrace(g:f),'.','Color',clr); hold on
                
            end
            
        end
        
    end
    
    
    % Plot posteriors
    for n = 1:2
        
        if plot_remoteW
            if n == 1
                ax(2) = subplot(1,3,2);
            elseif n == 2
                ax(3) = subplot(1,3,3);
            end
        else
            if n == 1
                ax(2) = subplot(12,1,[6:8]);
            else
                continue
            end
        end
        
        
        imageblock = P.posteriors{n}(:,a:b);
        
        if plot_onlyminspike
            blackoutinds = ~ismember(a:b,P.activebins{1});
            imageblock(:,blackoutinds) = 0;
        end
        
        if ~omitposterior
            imagesc(timevec(a:b),linposbins,imageblock,[0 .05]);
        end
        hold on
        colormap hot
        
        %plot position on local plot
        if n == 1
           pos_a = lookup(starttime,postimevec - postimevec(1));
           pos_b = lookup(endtime,postimevec - postimevec(1));
           if ~omitposterior
               armdistclr = [.7 .7 .7];
           else
              armdistclr = [0 0 0]; 
           end
           plot(postimevec(pos_a:pos_b) - postimevec(1),armdists(pos_a:pos_b),'.','Color',armdistclr,'Parent',ax(n+1))
           ylabel('Position','fontsize',16)
           set(gca,'fontsize',16)
        end
        
        % plot lines demarcating positions of arms
        % between Center and Right
        plot([-timevec(a) timevec(b)],[centerarmmax(n) centerarmmax(n)],'--','Color',[.8 .8 .8],'linewidth',.5,'Parent',ax(n+1))
        % between Right and Left
        plot([-timevec(a) timevec(b)],[centerarmmax(n)+rightarmmax(n) ...
            centerarmmax(n)+rightarmmax(n)],'--','Color',[.8 .8 .8],'linewidth',0.5,'Parent',ax(n+1))
        

        
    end
    
    
    % Plot speed + vte
    ax(3) = subplot(12,1,9:10);
    
    % plot vte
    if 0
    if plot_vte
        for vv = 1:length(vte_start)
            patch([vte_start(vv) vte_start(vv) vte_end(vv) vte_end(vv)] - postimevec(1),...
                [0 50 50 0],...
                [.95 1 .95],'edgecolor','none','Parent',ax(3)); hold on
        end
    end
    end
        
    postimevec2 = postimevec - postimevec(1);
    pos_a = lookup(starttime,postimevec2);
    pos_b = lookup(endtime,postimevec2);
    plot(postimevec2(pos_a:pos_b),veltrace(pos_a:pos_b),'-k','linewidth',3);
    ylabel('Speed','fontsize',16)
    set(gca,'fontsize',16)
    ylim([0 50])
    
    % Plot Head direction
    ax(4) = subplot(24,1,21:22);
    if plot_vte
    for vv = 1:length(vte_start)
        patch([vte_start(vv) vte_start(vv) vte_end(vv) vte_end(vv)] - postimevec(1),...
            [-pi pi pi -pi],...
            [.9 1 .9],'edgecolor','none','Parent',ax(4)); hold on
    end
    end
    postimevec2 = postimevec - postimevec(1);
    pos_a = lookup(starttime,postimevec2);
    pos_b = lookup(endtime,postimevec2);
    plot(postimevec2(pos_a:pos_b),dirtrace(pos_a:pos_b),'-','Color',[.5 .5 .5],'linewidth',3,'Parent',ax(4)); 
    ylabel('Headdir','fontsize',16)
    set(gca,'fontsize',16)
    set(gca,'ytick',[-pi 0 pi])
     set(gca,'yticklabel',{'pi','0','pi'})
    ylim([-pi pi])
    
    % Plot Difference in decode positions from actual
    if ~omitposterior
    if plot_nonlocalbins
        ax(5) = subplot(24,1,23:24);
            timevec = P.binvec_c{1} - P.binvec_c{1}(1) + binsize/2;
        inda = lookup(starttime,timevec);
        indb = lookup(endtime,timevec);
        plot(timevec(inda:indb),P.diffdist(inda:indb),'-','Color',[.5 .5 .5],'linewidth',3,'Parent',ax(5));
        ylim([0 175])
    end
    end
    
    % link axes
    linkaxes(ax,'x');
    
end




if plot_events

    animpref = animal_toplot(1:3);
        animalinfo = animaldef(animal_toplot);
        
    d = dayep(1);
    ep = dayep(2);
    if loaddata
        cd(datadir)
        load(sprintf('%s_%d_%d',animpref,d,ep),'P')
        pos = loaddatastruct(animalinfo{2},animalinfo{3},'pos',d);
            % obtain times from pos
            ep_start = pos{d}{ep}.data(1,1);
            ep_end = pos{d}{ep}.data(end,1);
            if 1
                out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'rippletrace', d, ep);
                riptrace = zscorer(out{d}{ep}{1}.powertrace);
                riptrace_timevec = out{d}{ep}{1}.eegtimesvec_ref;
                out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'fastgammatrace', d, ep);
                wgtrace = zscorer(out{d}{ep}{2}.powertrace);
                wgtrace_timevec = out{d}{ep}{2}.eegtimesvec_ref;
            end
    end         
    
    numevents = P.numevents;
    eventtimes = P.eventtimes;
    spike_mats = P.spike_mats;
    maxtet = max(P.selected_tets);
    mstimevec_ep = round(ep_start):.001:round(ep_end);
    
    for ww = 1:numevents

    %figure; hold on
    
    ev_starttime = eventtimes(ww,1);
    ev_endtime = eventtimes(ww,2);
    
    mstimevec_win = P.xvecms/1000 + ev_starttime;  
    
    % ms time of event (e.g. ripple) start
    ev_startms = lookup(ev_starttime,mstimevec_win);
    ev_endms = lookup(ev_endtime,mstimevec_win);
    
    % Spike raster plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3,3,[1 4 7]);
    
    % plot 2 SD SWRs that occur in the widnow %%%%%%%%%
        % identify which 2 SD ripples to plot
        % first find 2 SD ripples again (brief calculation)
    ripout = kk_getconstimes(animalinfo{2},animalinfo{3},[d ep],'ripplescons',1,...
        'consensus_numtets',3,'minthresh',2,...
        'exclusion_dur',0,'minvelocity',0,'maxvelocity',4);
    consvec_rip2 = ripout{d}{ep}.cons;
    consvectimes_rip2 = ripout{d}{ep}.time;
    periodtimes_rip2 = vec2list(consvec_rip2,consvectimes_rip2);
    plotwinvec = logical(list2vec([ev_starttime-0.5   ev_starttime+0.5],ripout{d}{ep}.time))';
    consvec_rip2_toplot = consvec_rip2 & plotwinvec;
    rip_toplot = 1000 * (vec2list(consvec_rip2_toplot,ripout{d}{ep}.time) - ev_starttime);
    numrip_toplot = size(rip_toplot,1);
    % plot all 2 SD ripples in the raster window
    for rp = 1:numrip_toplot
        patch([rip_toplot(rp,1) rip_toplot(rp,2) rip_toplot(rp,2) rip_toplot(rp,1)],...
            [1 1 (maxtet + 4) (maxtet + 4)],...
            [1 .9 .9],'edgecolor','none'); hold on
    end

    % plot EVENT as a thick red line (i.e. ripple or wave gamma) %%%%%%%%%%%
    hold on
    plot([ev_startms ev_endms] - ev_startms,[.5 .5],'Color','r','linewidth',4)
    if 0
        set(gca,'xtick',[-300:100:300]);
        set(gca,'xticklabel',{'-300','','','0','','','+300'});
        xlim([-300 300]);
    else
        set(gca,'xtick',[-500:100:500]);
        set(gca,'xticklabel',{'-500','','','','','0','','','','','+500'});
        xlim([-500 500]);
    end
    
    % plot SPIKES raster %%
    for tet = 1:maxtet
        spikebins = find(spike_mats(tet,:,ww) > 0);
        numlines = length(spikebins);
        for ll = 1:numlines
            plot([spikebins(ll) spikebins(ll)] - ev_startms,[tet-0.5 tet+0.5],...
                'linestyle','-','Color',[0 0 0],'LineWidth',1)
        end
    end
    set(gca,'ytick',1:maxtet)
    %set(gca,'YTick',1:5:30,'YTickLabel',1:5:30,'FontSize',14);
    xlabel('Time (s)') ;   ylabel('Tetrode');
    title(sprintf('%s %d %d %s # : %d',animpref,d,ep,eventname(1:(end-1)),ww),...
        'fontweight','bold','fontsize',14);
    set(gca,'tickdir','out');
    
    box off
    
    % plot WG and SWR power traces at bottom%%%%%%%%%%%%%%%%
    a = lookup(ev_starttime - 0.5,riptrace_timevec) ;
    b = lookup(ev_starttime + 0.5,riptrace_timevec) ;
    plot(1000 * (riptrace_timevec(a:b) - ev_starttime), -riptrace(a:b)/2 + maxtet + 4,'r-','linewidth',2); hold on
    a = lookup(ev_starttime - 0.5,wgtrace_timevec) ;
    b = lookup(ev_starttime + 0.5,wgtrace_timevec) ;
    plot(1000 * (wgtrace_timevec(a:b) - ev_starttime), 2 * -wgtrace(a:b)/2 + maxtet + 4,'-','Color',[0 .5 1],'linewidth',2);
    
    set(gca,'ydir','reverse')
    set(gca,'fontsize',12)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Decoded posterior image plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for n = 1:2
        if n == 1
            subplot(3,3,[2 5 8]);
        elseif n == 2
            subplot(3,3,[3 6 9]);
        end
        linposbins = P.linposbins{n} ;
        imagesc(P.xvecms,linposbins,P.posteriors{n}(:,:,end));
        
        %title('postx, clusterless','FontSize',12);
        ylabel('Linearized position','FontSize',12);
        %set(gca,'FontSize',14,'YTick',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3],'YTickLabel',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]);
        %colormap(flipud(hot(256)));
        colormap(hot);
        caxis([0 0.1]);
        ylim([-.01 max(linposbins)])
        set(gca,'tickdir','out');
        
        set(gca,'xtick',[-500:100:500]);
        set(gca,'xticklabel',{'-500','','','','','0','','','','','+500'});
        xlim([-500 500])
        set(gca,'fontsize',12)
        xlabel('Time (ms)')
        
        hold on
        % in local environment, plot current position (dark grey line)
        if n == 1
         %   winpos = armdists_cat{n}( ev_startind:ev_endind );  % armdists_cut (linear) positions over the course of the ripple
         %   winpos_timevec = 1000 * ( postimevec{1}(ev_startind:ev_endind) - ev_starttime );  % 1-ms indices over the course of the ripple
         %   plot(winpos_timevec,winpos,'-','linewidth',6,'Color',[.8 .8 .8])
        end
        
        % plot line indicating start and end of ripple
        %plot([ripms_start ripms_end] - ripms_start,[50 50],'Color','r','linewidth',4)
        plot([0 0],[0 max(linposbins)],'Color','r','linewidth',1)
        plot([ev_endms - ev_startms ev_endms - ev_startms],...
            [0 max(linposbins)],'Color','r','linewidth',1)
        
        % plot lines demarcating positions of arms
        % between Center and Right
        plot([-500 500],[centerarmmax(nn) centerarmmax(nn)],'--','Color','w','linewidth',2)
        % between Right and Left
        plot([-500 500],[centerarmmax(nn)+rightarmmax(nn) ...
            centerarmmax(nn)+rightarmmax(nn)],'--','Color','w','linewidth',2)
    end
    
    end
end


    
    




