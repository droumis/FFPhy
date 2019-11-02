% Cluster decoding on W-track -- psuedo 1D "armdists" (Wu--Foster-2014 Fig 1)


datadir = '/opt/data13/kkay/Superclustdecode_dir_data';

plot_events = 1;
    if plot_events
       animal_toplot = 'Bond';
            animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
       eventname = 'ripples';
       dayep = [4 6];
        minthresh_rip = 5;
       loaddata = 1;
            % other
       plot_remoteW = 0;
       tracekern_ms = 10;       
       omitraster = 0;
       omitposterior = 0;
       plot_nonlocalbins = 1;    
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
            postimevec = pos{d}{ep}.data(:,1);
            if 0
                out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'rippletrace', d, ep);
                riptrace = zscorer(out{d}{ep}{1}.powertrace);
                riptrace_timevec = out{d}{ep}{1}.eegtimesvec_ref;
                out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'fastgammatrace', d, ep);
                wgtrace = zscorer(out{d}{ep}{2}.powertrace);
                wgtrace_timevec = out{d}{ep}{2}.eegtimesvec_ref;
            elseif 1
                out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'rippletrace', d, ep);
                [riptrace,~,~] = zscoretrace(out{d}{ep}{1},[]);
                riptrace_timevec = out{d}{ep}{1}.eegtimesvec_ref;
                out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'fastgammatrace', d, ep);
                [wgtrace,~,~] = zscoretrace(out{d}{ep}{2},tracekern_ms);
                wgtrace_timevec = out{d}{ep}{2}.eegtimesvec_ref;           
                out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'lowgammatrace', d, ep);
                [lgtrace,~,~] = zscoretrace(out{d}{ep}{1},tracekern_ms);
                lgtrace_timevec = out{d}{ep}{1}.eegtimesvec_ref;          
                out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'highgammatrace', d, ep);
                [hgtrace,~,~] = zscoretrace(out{d}{ep}{1},tracekern_ms);
                hgtrace_timevec = out{d}{ep}{1}.eegtimesvec_ref;                
            end
    end
    
    % Clustered spike data
    spikes = loaddatastruct(animalinfo{2},animalinfo{3},'spikes',d);
    pos = loaddatastruct(animalinfo{2},animalinfo{3},'pos',d);
        epstart = pos{dayep(1)}{dayep(2)}.data(1,1);
    cellinfo = loaddatastruct(animalinfo{2},animalinfo{3},'cellinfo');
        
    % Identify 2 SD SWRs (here, to use to exclude from encoding model)
    ripout = kk_getconstimes(animalinfo{2},animalinfo{3},[d ep],'ripplescons',1,...
                             'consensus_numtets',3,'minthresh',minthresh_rip,...
                             'exclusion_dur',0,'minvelocity',0,'maxvelocity',4);
    ripconsvec = ripout{d}{ep}.cons;
    consvectimes = ripout{d}{ep}.time;
    riptimes = vec2list(ripconsvec,consvectimes);
    
    
    % Define times to evaluate replay
    eventtimes = riptimes;
    
    numevents = size(eventtimes,1);
        disp(sprintf('Plotting %d events',numevents))
    
    %spike_mats = P.spike_mats;

    mstimevec_ep = round(ep_start):.001:round(ep_end);
    
    for ww = 1:numevents
        
        figure('units','normalized','outerposition',[0 0 .5 .5]);  hold on
        
        ev_starttime = eventtimes(ww,1);
        ev_endtime = eventtimes(ww,2);
            eventtime = round(ev_starttime - epstart);  % time of event, from beginning of epoch
        
        win_start = ev_starttime - 0.5;
        win_end = ev_starttime + 0.5;        
        
        mstimevec_win = P.binvec_c + ev_starttime;
        
        % ms time of event (e.g. ripple) start
        ev_startms = lookup(ev_starttime,mstimevec_win);
        ev_endms = lookup(ev_endtime,mstimevec_win);
        
        % Spike raster plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for nn = 1:2
  
            numunits = size(P.detc{nn},1);
            detc = P.detc{nn};            
            
            if nn == 1
                subplot(3,4,[1 5 9]);
            elseif nn == 2
                subplot(3,4,[3 7 11]);
            end
        
        plotwinvec = logical(list2vec([ev_starttime-0.5   ev_starttime+0.5],consvectimes))';
        consvec_rip2_toplot = ripconsvec & plotwinvec;
        rip_toplot = 1000 * (vec2list(consvec_rip2_toplot,consvectimes) - ev_starttime);
        numrip_toplot = size(rip_toplot,1);
        % plot all 2 SD ripples in the raster window
        for rp = 1:numrip_toplot
            patch([rip_toplot(rp,1) rip_toplot(rp,2) rip_toplot(rp,2) rip_toplot(rp,1)],...
                [0 0 (numunits + 4) (numunits + 4)],...
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
        
        %     % plot SPIKES raster %%
        %     for tet = 1:numunits
        %         spikebins = find(spike_mats(tet,:,ww) > 0);
        %         numlines = length(spikebins);
        %         for ll = 1:numlines
        %             plot([spikebins(ll) spikebins(ll)] - ev_startms,[tet-0.5 tet+0.5],...
        %                 'linestyle','-','Color',[0 0 0],'LineWidth',1)
        %         end
        %     end
        % Unit spike raster
        if ~omitraster
            %ax(1) = subplot(12,1,[1:5]);
            celllist = {};
            % iterate
            for cc = 1:numunits
                tet = detc(cc,3);
                cellnum = detc(cc,4);
                region = cellinfo{d}{ep}{tet}{cellnum}.area;
                [~,~,~,bgclr] = regionfunc(region);
                type = cellinfo{d}{ep}{tet}{cellnum}.type;
                celllist = [celllist ; num2str(cc) type region];
                %disp(type)
                if ~strcmp(type,'principal')
                    bgclr = [1 1 1];
                end
                
                % plot region (CA1, CA3) color patch for region
                %patch([starttime starttime endtime endtime],[cc-0.5 cc+0.5 cc+0.5 cc-0.5],bgclr,'edgecolor','none');
                %    hold on
                
            end
        end
        % plot spike raster
        if ~omitraster
            for cc = 1:numunits
                tet = detc(cc,3);
                cellnum = detc(cc,4);
                region = cellinfo{d}{ep}{tet}{cellnum}.area;
                [~,~,~,bgclr] = regionfunc(region);
                type = cellinfo{d}{ep}{tet}{cellnum}.type;
                %disp(type)
                if ~strcmp(type,'principal')
                    bgclr = [1 1 1];
                end
                % plot spikes
                if ~isempty(spikes{d}{ep}{tet}{cellnum}.data)
                    spiketimes = spikes{d}{ep}{tet}{cellnum}.data(:,1);
                    spiketimes = spiketimes(logical(isExcluded(spiketimes,[win_start win_end]))) ;  % spikes within window
                    numunitspk = length(spiketimes);
                    if numunitspk > 0
                        for ll = 1:numunitspk
                            plot(1000 * ([spiketimes(ll) spiketimes(ll)] - ev_starttime),[cc-0.48 cc+0.48],...
                                'linestyle','-','Color',[0 0 0],'LineWidth',1 ); hold on
                        end
                    end
                end
            end
        else
            set(gca,'ytick',0)
        end
        
        
        set(gca,'ytick',1:numunits)
        %set(gca,'YTick',1:5:30,'YTickLabel',1:5:30,'FontSize',14);
        xlabel('Time (s)') ;   ylabel('Unit #');
            
        title(sprintf('%s %d %d %s # : %d (time: %d s)',animpref,d,ep,eventname(1:(end-1)),ww,eventtime),...
            'fontweight','bold','fontsize',14);
        set(gca,'tickdir','out');
        
        box off

        % plot WG and SWR power traces at bottom%%%%%%%%%%%%%%%%
        a = lookup(win_start,riptrace_timevec) ;
        b = lookup(win_end,riptrace_timevec) ;
        plot(1000 * (riptrace_timevec(a:b) - ev_starttime), -3 * riptrace(a:b)/2 + numunits + 4,'-','Color',[.7 .7 .7],'linewidth',2); hold on
        a = lookup(win_start,wgtrace_timevec) ;
        b = lookup(win_end,wgtrace_timevec) ;
        plot(1000 * (wgtrace_timevec(a:b) - ev_starttime), 3 * -wgtrace(a:b)/2 + numunits + 4,'-','Color',[0 .5 1],'linewidth',2);
        a = lookup(win_start,lgtrace_timevec) ;
        b = lookup(win_end,lgtrace_timevec) ;
        plot(1000 * (lgtrace_timevec(a:b) - ev_starttime), 3 * -lgtrace(a:b)/2 + numunits + 4,'-','Color',[1 .5 .5],'linewidth',2);
        a = lookup(win_start,hgtrace_timevec) ;
        b = lookup(win_end,hgtrace_timevec) ;
        plot(1000 * (hgtrace_timevec(a:b) - ev_starttime), 3 * -hgtrace(a:b)/2 + numunits + 4,'-','Color',[.2 .6 .2],'linewidth',2);
        
        set(gca,'ydir','reverse')
        set(gca,'fontsize',12)
        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ind_a = lookup(win_start,P.binvec_c) ;
        ind_b = lookup(win_end,P.binvec_c) ;
        
        % Decoded posterior image plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for n = 1:length(P.linposbins)
            
            if n == 1
                subplot(3,4,[2 6 10]);
            elseif n == 2
                subplot(3,4,[4 8 12]);
            end
            linposbins = P.linposbins{n} ;
            imagesc(1000*(P.binvec_c(ind_a:ind_b) - P.binvec_c(ind_a))-500,linposbins,P.posteriors{n}(:,ind_a:ind_b));
            
            %title('postx, clusterless','FontSize',12);
            ylabel('Linearized position','FontSize',12);
            %set(gca,'FontSize',14,'YTick',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3],'YTickLabel',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]);
            %colormap(flipud(hot(256)));
            colormap(gray); colormap(flipud(colormap))
            caxis([0 0.1]);
            ylim([-.01 max(linposbins)]);
            set(gca,'tickdir','out');
            
            set(gca,'xtick',[-500:100:500]);
            set(gca,'xticklabel',{'-500','','','','','0','','','','','+500'});
            xlim([-500 500])
            set(gca,'fontsize',12)
            xlabel('Time (ms)')
            
            hold on
            % in local environment, plot current position (dark grey line)
            if n == 1
                j = lookup(win_start,postimevec);
                k = lookup(win_end,postimevec);
                   winpos = P.armdists( j:k );  % armdists_cut (linear) positions over the course of the ripple
                   plot(1000*(postimevec(j:k)-postimevec(j))-500,winpos,'-','linewidth',6,'Color',[.8 .8 .8]);
            end
            
            % plot line indicating start and end of ripple
            %plot([ripms_start ripms_end] - ripms_start,[50 50],'Color','r','linewidth',4)
            plot([0 0],[0 max(linposbins)],'Color','r','linewidth',1)
            plot([ev_endms - ev_startms ev_endms - ev_startms],...
                [0 max(linposbins)],'Color','r','linewidth',1)
            
            % plot lines demarcating positions of arms
            % between Center and Right
            plot([-500 500],[P.centerarmmax(n) P.centerarmmax(n)],'--','Color','k','linewidth',2)
            % between Right and Left
            plot([-500 500],[P.centerarmmax(n)+P.rightarmmax(n) ...
                        P.centerarmmax(n)+P.rightarmmax(n)],'--','Color','k','linewidth',2)
                   
            % plot Dividing line between inbound (toward center junct) vs.
            % outbound (away from center junct)
            divide_dist = max(P.linposbins{n})/2;
            plot([-500 500],[divide_dist divide_dist],'-','Color','k','linewidth',4)
                % second set of armdist demarcations
            %between Center and Right
            plot([-500 500],[P.centerarmmax(n) P.centerarmmax(n)] + divide_dist,'--','Color','k','linewidth',2)
            % between Right and Left
            plot([-500 500],[P.centerarmmax(n)+P.rightarmmax(n) ...
                        P.centerarmmax(n)+P.rightarmmax(n)] + divide_dist,'--','Color','k','linewidth',2)            
            
            
                    
        end
            pause
            close all
            
    end
    

end


    
    




