datadir = '/opt/data13/kkay/Superposteriors_data';

calculate = 0;
plot_continuous = 1;
    if plot_continuous 
       animal_toplot = 'Bond';
            animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
       dayep = [4 4];
       starttime = 64   ;  % time (within epoch, not clock time) to begin plot
       endtime = 78;
       omit_tetraster = 1;
       plot_hg = 0;
       tracekern_ms = [1000];
    end
plot_events = 0;
    if plot_events
       animal_toplot = 'Bond';
       eventname = 'ripples';
       dayep = [3 2];
    end
if calculate
    %%% select data %%
    animals_tocalc = {'Bond'};  %{'Bond'}; {'Frank','Bond','Government','Egypt','Dave'};
    tetfilter =  '( isequal($area, ''CA1'') || isequal($area, ''CA3'') || isequal($area, ''CA2'') || isequal($area, ''DG'')   )';
    dayeps = [4 4 ]; %[3 2; 3 4; 3 6 ; 5 2 ; 5 4 ; 5 6 ];  %[6 2; 6 4; 6 6];  % leave empty if want to do all epochs
    remoteW = 0;   % set to 1 if want to decode in other W as well
    plot_infunction = 0;
        selected_events = [];
    %%% select what to decode %%
    decodemode = 1; % 1: entire epoch, 2: SWRs
    extratime = 500;  % if not doing decodemode 1, the ms around event start to plot
    %%% select decoding parameters %%
    old_version = 1;
    modelnum = 3;   % 1: all spikes, 2: trajencode, 3: exclude SWR only
    spikethresh = 100;  % leave [] if want to use the amplitudes in the spike data
    dt = .001; 0.0334; %postimevec(2) - postimevec(1);
    sigma_transmat = 5 ; % sigma for the gaussian used to smooth the transition matrix
    mdel = 2;  % uV, spacing in amplitude mark space
    xdel = 1;  % cm, spacing in positional mark space
        smker =  6 * mdel;  % gaussian kernel sigma in amplitude mark space
        sxker = 3 * xdel;   % gaussian kernel sigma in positional mark space
end

if calculate
    for aa = 1:length(animals_tocalc)
        animalname = animals_tocalc{aa};
            animalinfo = animaldef(animalname);
        kk_clusterlessdecode2(animalinfo{2},animalinfo{3},dayeps,animalname,...
            'savedir',datadir,'decodemode',decodemode,'tetfilter',tetfilter,'modelnum',modelnum,'spikethresh',spikethresh,'old_version',old_version,...
            'plot_powertraces',1,'calctraj',1,'xdel',xdel,'mdel',mdel,...
            'sxker',sxker,'smker',smker,'extratime',extratime,'sigma_transmat',sigma_transmat,...
            'dt',dt,'remoteW',remoteW,'plot_infunction',plot_infunction);
    end
end

if plot_continuous
    
    animpref = animal_toplot(1:3);
        animalinfo = animaldef(animal_toplot);
        daydir = getdaydir(animal_toplot);   
        an = find(strcmp(animal_toplot,animals_order));
 
    d = dayep(1);
    ep = dayep(2);
    a = round(starttime*1000);
    b = round(endtime*1000);
        num_mssteps = b - a + 1 ;

    % load data    
    if ~exist('D','var') || ~strcmp(D{1}.animalname,animal_toplot) && ...
                                         ~all(D{1}.dayep == dayep)
        
        % posteriors
        D = {};
        cd(datadir)
        for nn = 1 %:2
            load(sprintf('%s_%d_%d_fullepoch_%d_thresh100.mat',animpref,d,ep,nn),'P');
            D{nn} = P;
        end
        
        % position
        pos = loaddatastruct(animalinfo{2},animalinfo{3},'pos',d);
        postimevec = pos{d}{ep}.data(:,1);
        veltrace = pos{d}{ep}.data(:,5);
        epstart = postimevec(1);
        epend = postimevec(end);
        
        % clustered spikes
        spikes = loaddatastruct(animalinfo{2},animalinfo{3},'spikes',d);
        cellinfo = loaddatastruct(animalinfo{2},animalinfo{3},'cellinfo');
                [adtc_list] = clusteredunits(spikes,cellinfo,an,d,ep,[],1);
                maxcell = size(adtc_list,1);   
        
        
        % ripples, 2 SD
        ripout = kk_getconstimes(animalinfo{2},animalinfo{3},[d ep ],'ripplescons',1,...
            'consensus_numtets',3,'minthresh',2,...
            'exclusion_dur',0,'minvelocity',0,'maxvelocity',4);
        consvec_rip2 = ripout{d}{ep}.cons;
        consvectimes_rip2 = ripout{d}{ep}.time;
        periodtimes_rip2 = vec2list(consvec_rip2,consvectimes_rip2);
        

        

    end
    
        % Ripple and WG power traces
        out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'rippletrace', d, ep);
        riptrace = zscorer(out{d}{ep}{1}.powertrace);
        riptrace_timevec = out{d}{ep}{1}.eegtimesvec_ref;
        out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'fastgammatrace', d, ep);
        wgtrace = zscorer(out{d}{ep}{2}.powertrace);
        wgtrace_timevec = out{d}{ep}{2}.eegtimesvec_ref;
        out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'lowgammatrace', d, ep);
        lgtrace_ca1 = zscorer(out{d}{ep}{1}.powertrace);  % ca1
        lgtrace_ca3 = zscorer(out{d}{ep}{2}.powertrace);  % ca3 dg
        lgtrace_timevec = out{d}{ep}{2}.eegtimesvec_ref;        
        if plot_hg
            out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'highgammatrace', d, ep);
            hgtrace = zscorer(out{d}{ep}{1}.powertrace);
        end    
    
    if ~isempty(tracekern_ms)
        sig_bins = 1500 * tracekern_ms/1000;
        gauskern = gaussian(sig_bins,8*sig_bins);
        wgtrace = smoothvect(wgtrace,gauskern);
        lgtrace_ca1 = smoothvect(lgtrace_ca1,gauskern);
        lgtrace_ca3 = smoothvect(lgtrace_ca3,gauskern);
        if plot_hg
            hgtrace = smoothvect(hgtrace,gauskern);
        end
        %riptrace = smoothvect(riptrace,gauskern);
    end
    
    
    
    selected_tets = D{1}.selected_tets;
        maxtet = max(selected_tets);
    spikethresh = D{1}.spikethresh;
    %xvecms = D{1}.xvecms;
    xvecms = round(epstart)*1000:1:round(epend)*1000 - epstart;
    armdists = D{1}.armdists;
    linposbins = D{1}.linposbins{1};
    centerarmmax = D{1}.centerarmmax;
    rightarmmax = D{1}.rightarmmax;
    
    figure;
    
    % Unit regional patches
    ax(1) = subplot(6,1,[1 2 3]);
    celllist = {};
    % iterate
    for cc = 1:maxcell
        tet = adtc_list(cc,3);
        cellnum = adtc_list(cc,4);
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
    
    % First plot ripple patches and wg traces
        for xx = 2  %:2
            ymax = 20;
            if xx == 1
                ax(1) = subplot(6,1,[1 2 3]); 
                ylevel = maxtet;
            elseif xx == 2
                % obtain all clustered units this epoch
                an = find(strcmp(animal_toplot,animals_order));
                [adtc_list] = clusteredunits(spikes,cellinfo,an,d,ep,[],1);
                maxcell = size(adtc_list,1);
                ax(1) = subplot(6,1,[1 2 3]); 
                ylevel = maxcell;
            end
            % plot 2 SD SWRs that occur in the widnow %%%%%%%%%
            plotwinvec = logical(list2vec([starttime endtime] + epstart,ripout{d}{ep}.time))';
            consvec_rip2_toplot = consvec_rip2 & plotwinvec;
            rip_toplot = vec2list(consvec_rip2_toplot,ripout{d}{ep}.time) - epstart;
            numrip_toplot = size(rip_toplot,1);
            % plot all 2 SD ripples in the raster window
            for rp = 1:numrip_toplot
                patch([rip_toplot(rp,1) rip_toplot(rp,2) rip_toplot(rp,2) rip_toplot(rp,1)] ,...
                    [0 0 (ylevel + ymax) (ylevel + ymax)],...
                    [1 .9 .9],'edgecolor','none','Parent',ax(1) ); hold on
            end
            % plot WG and SWR power traces at bottom%%%%%%%%%%%%%%%%
            ind1 = lookup(starttime,riptrace_timevec - epstart) ;
            ind2 = lookup(endtime,riptrace_timevec - epstart) ;
            ind3 = lookup(starttime,wgtrace_timevec - epstart) ;
            ind4 = lookup(endtime,wgtrace_timevec - epstart) ;
            plot((riptrace_timevec(ind1:ind2) - epstart), 4 * -riptrace(ind1:ind2)/2 + ylevel + 10,'-','Color',[.8 .75 .75],'linewidth',2,'Parent',ax(1) ); hold on
            plot((wgtrace_timevec(ind3:ind4) - epstart), 8 * -wgtrace(ind3:ind4)/2 + ylevel + 10,'-','Color',[0 .3 .9],'linewidth',2,'Parent',ax(1) ); hold on
            plot((lgtrace_timevec(ind1:ind2) - epstart), 8 * -lgtrace_ca1(ind1:ind2)/2 + ylevel + 10,'-','Color',[1 .1 .1],'linewidth',2,'Parent',ax(1) ); hold on
            %plot((lgtrace_timevec(ind1:ind2) - epstart), 4 * -lgtrace_ca3(ind1:ind2)/2 + ylevel + 10,'-','Color',[1 .4 .4],'linewidth',2,'Parent',ax(1) ); hold on
            
            
            xlim([starttime endtime])
            ylim([0  (ylevel + ymax)])
            set(gca,'ydir','reverse')
            set(gca,'fontsize',12)
            set(gca,'tickdir','out');
            xlabel('Time (ms)') ;   ylabel('unit #');
        end    
    
    % Tetrode spike raster
    if 0
    ax(1) = subplot(2,3,1);  
    if ~omit_tetraster
    spike_mats = zeros(maxtet,num_mssteps);
    for tet = selected_tets
        filedata = loadparamsfile(daydir,d,tet);  % spike data from day directory   
        inds_ep =  (  filedata.params(:,1)/10000  >=  epstart  )  & ...  %  Spikes in epoch
                      ( filedata.params(:,1)/10000 <= epend );
        inds_thresh = any(filedata.params(:,2:5) > spikethresh, 2) ;  %  Spikes with min amplitude
            spiketimes = filedata.params(inds_thresh & inds_ep,1)/10000;
        spiketimes = spiketimes(logical(isExcluded(spiketimes,[starttime endtime]+epstart))) ;  % spikes within window
        numspk = length(spiketimes);
        for ll = 1:numspk
            plot([spiketimes(ll) spiketimes(ll)] - epstart,[tet-0.48 tet+0.48],...
                'linestyle','-','Color',[0 0 0],'LineWidth',1,'Parent',ax(1) ); hold on
        end
    end
    xlim([starttime endtime])
    %ylim([0  (maxtet + 7)])    
    %set(gca,'ytick',1:maxtet,'fontsize',8)
    xlabel('Time (s)') ;   ylabel('Tetrode');
    title(sprintf('%s %d %d (%d - %d s)',animpref(1:3),d,ep,starttime,endtime),...
        'fontweight','bold','fontsize',14);
    set(gca,'tickdir','out');
    set(gca,'ydir','reverse')
    set(gca,'fontsize',12)
   
    box off
    end    
    end
    
    
    % Unit spike raster
    ax(1) = subplot(6,1,[1 2 3]); 
        % iterate
        for cc = 1:maxcell
            tet = adtc_list(cc,3);
            cellnum = adtc_list(cc,4);
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
        
        

  
    
    % Plot posteriors
    for n = 1   %:2
        
        if n == 1
            %ax(3) = subplot(2,3,[2 3 5 6]);
            ax(2) = subplot(6,1,[4 5]);
        elseif n == 2
            ax(4) = subplot(2,3,[3 6]);
        end
        
        xvec = (round(starttime*1000):1:round(endtime*1000))/1000;
       
        imagesc(xvec,linposbins,D{n}.posteriors{n}(:,a:b),[0 .1]);
        hold on
        colormap hot
        
        %plot position on local plot
        if n == 1
           pos_a = lookup(starttime,postimevec - postimevec(1));
           pos_b = lookup(endtime,postimevec - postimevec(1));
           plot(postimevec(pos_a:pos_b) - postimevec(1),armdists(pos_a:pos_b),'.','Color',[.7 .7 .7],'Parent',ax(2))
        end
        
        % plot lines demarcating positions of arms
        % between Center and Right
        plot([-xvecms(a) xvecms(b)],[centerarmmax(n) centerarmmax(n)],'--','Color','w','linewidth',2,'Parent',ax(2))
        % between Right and Left
        plot([-xvecms(a) xvecms(b)],[centerarmmax(n)+rightarmmax(n) ...
            centerarmmax(n)+rightarmmax(n)],'--','Color','w','linewidth',2,'Parent',ax(2))
        
        
        
        
    end

    
    % plot speed
    ax(3) = subplot(6,1,6);
    
    postimevec2 = postimevec - postimevec(1);
    pos_a = lookup(starttime,postimevec2);
    pos_b = lookup(endtime,postimevec2);
    plot(postimevec2(pos_a:pos_b),veltrace(pos_a:pos_b),'-k','linewidth',3);
    
        
    
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
        load(sprintf('%s_%d_%d_%s',animpref,d,ep,eventname),'P')
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


    
    




