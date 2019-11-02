datadir = '/opt/data50/kkay/__Decode';
datadir = '/opt/data50/kkay/__Decode/Randwalk';


plotter = 1;
analyze = 0;
plot_placemap = 0;
place_map_only = 0;
    EXAMPLE = 100;   % 100, 1005, 1001
    MAXVAL = 0.05;
    PHASECHOICE = 1;
    clusterless_examples;  % <<< select examples here
    if place_map_only
        plotter = 0;
        analyze = 0;
        plot_placemap = 1;
    end
    MOVINGPERIODS  = 3 ;  % 1: vel4, 2: nonimmobile05, 3: nonimmobile1
    omit_mainplot = 0;
    

    if analyze
        ELIMINATE_LR_EXCLUSIVE_MINPHASE = 0;
        EXCURPERIODS = 1;  % 1: Excursions + <CPbuff // 2: Head out + <CP buff
        % Spatial index
        LR_prop_thresh = 0.1;   % minimum posterior density that is either in L or R (versus C)        
        % Alternation detection parameters
        ALT_THRESHOLD       = 0.5;
        NUMSAMP            = 3000 + 1;  % # of resampling methods
            SAMPMETHOD = 1;   % 1: permutation, 2: bootstrap
        ALTDISTBINS         = 0:0.05:1;
        MAXALTS             = 3;    
        P_THRESH            = 0.1;
    end
    if plotter 
        % Study figures %%%%%%%%%%%%
        plot_Xcorr_LR = 0;
        plot_LRprop_hist = 0;
        plot_LRprop_phasehist = 0;
        plot_altspeed_hist = 1;
        plot_pvalue_scatter = 1;
        plot_altspeed_counts = 1;        
        % LFP plot %%%%%%%%%%%%%%%
        LFPMAX = 500;
        plot_LFP = 0;
        plot_LOWFREQ_LFP = 0;      
        plot_DELTA_LFP = 0;  
    end
        % place plot %%%%%%%%%%%%%%%%%%%%%%%%%   
        ANGVEL_SCALE = 2;
        NUMSAMP_HEADPOS = 2 ;    
        
    % Data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotter || analyze
        
        cd(datadir)
        
        % loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander',...
            'Iceland','Justice','Kapital','Laplace','Eight','Ten','Conley','Miles'};
        animpref = animal_toplot(1:3);
        animalinfo = animaldef(animal_toplot);
        daydir = getdaydir(animal_toplot);
        an = find(strcmp(animal_toplot,animals_order));
        
        d = dayep(1);
        ep = dayep(2);

        % load Decode file
        if ~exist('P','var') || isempty(P) || ~strcmp(P.animalname,animal_toplot) || ...
               ~all(P.dayep == dayep) || P.EXTYPE ~= EXTYPE_toplot
           
           % Decode file %%%%%%%%%%%%%%
           filestr = sprintf('%s_%d_%d_Pos_%d.mat',animpref,d,ep,EXTYPE_toplot);
           disp(sprintf('loading %s',filestr))
           load(filestr,'P');
            
           % Positional %%%%%%%%%%%%%%%%%%%
           pos = loaddatastruct(animalinfo{2},animalinfo{3},'pos',P.dayep(1));
            [veltrace,hdir,xpos,ypos] = posparser(pos{d}{ep}.data);
           linpos = loaddatastruct(animalinfo{2},animalinfo{3},'linpos',P.dayep(1));
           lindist = linpos{d}{ep}.statematrix.lindist;
                CPbuff = choicepointer(linpos{d}{ep});
                CPprevec = lindist < CPbuff;
                CPvec = (lindist < CPbuff );   % 10 cm zone  (in green)  %
                Lvec =  ismember(linpos{d}{ep}.statematrix.segmentIndex,[4 5]);
                Rvec =  ismember(linpos{d}{ep}.statematrix.segmentIndex,[2 3]);
                
           postimevec = pos{d}{ep}.data(:,1);    
            epstart = postimevec(1);
            epend = postimevec(end);
            headang = loaddatastruct(animalinfo{2}, animalinfo{3}, 'headang', d);
                angvel = abs(headang{d}{ep}.angvel)  * 180 / pi;    % convert to degrees
            % tets, clust spikes
            spikes = loaddatastruct(animalinfo{2},animalinfo{3},'spikes',d);
            cellinfo = loaddatastruct(animalinfo{2},animalinfo{3},'cellinfo');
            tetinfo = loaddatastruct(animalinfo{2},animalinfo{3},'tetinfo');
            [adtc_list] = clusteredunits(spikes,cellinfo,an,d,ep,[],1);
            maxcell = size(adtc_list,1);
            selected_tets = unique(P.selected_tets);
                numtets = length(selected_tets);           
           
           % LFP %%%%%%%%%%%%%%%%%%%%%%%%%%%
            thetabins = loaddatastruct(animalinfo{2},animalinfo{3},'thetabins',d);
            thetatet = thetabins{d}.thetatet;
            if strcmp(animpref,'Dav') && d == 3
                disp('Dav thetatet = 1, d 3 ep 2');
                thetatet = 1;
            end
            theta = loadeegstruct(animalinfo{2},animalinfo{3},'theta',d,ep,thetatet);
                  thetatrace = theta{d}{ep}{thetatet}.data(:,1);
                  thetaphase = double(theta{d}{ep}{thetatet}.data(:,2))/10000;  % hilbert phase
                  thetatimevec = geteegtimes(theta{d}{ep}{thetatet});
                  if strcmp(animpref,'Cor')
                      LFPSIGN = 1;
                  else
                      LFPSIGN = -1;
                  end
            if plotter
                if plot_LFP
                    eeg = loadeegstruct(animalinfo{2},animalinfo{3},'eeg',d,ep,thetatet);
                    eegtrace = eeg{d}{ep}{thetatet}.data;
                    eegtimevec = geteegtimes(eeg{d}{ep}{thetatet});
                end
                if plot_LOWFREQ_LFP
                    lowfreq = loadeegstruct(animalinfo{2},animalinfo{3},'lowfreq',d,ep,thetatet);
                    lowfreqtrace = lowfreq{d}{ep}{thetatet}.data;
                    lowfreqtimevec = geteegtimes(lowfreq{d}{ep}{thetatet});
                end
                if plot_DELTA_LFP
                    delta = loadeegstruct(animalinfo{2},animalinfo{3},'delta',d,ep,thetatet);
                    deltatrace = delta{d}{ep}{thetatet}.data;
                    deltatimevec = geteegtimes(delta{d}{ep}{thetatet});
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           timefilterscript
           if MOVINGPERIODS == 1
               [dummy,~] = evaluatetimefilter(animalinfo{2},animalinfo{3},{vel4},[d ep]);   
           elseif MOVINGPERIODS == 2
               [dummy,~] = evaluatetimefilter(animalinfo{2},animalinfo{3},{nonimmobile05},[d ep]);    
           elseif MOVINGPERIODS == 3
               [dummy,~] = evaluatetimefilter(animalinfo{2},animalinfo{3},{nonimmobile1},[d ep]);   
           end
            movingperiods = dummy{d}{ep};   %
            movingvec = logical(list2vec(movingperiods,postimevec)); 
            
            % Process position data
            [posvec_all] = posvecmaker(lindist,linpos{d}{ep});
            [posvec_all,outersep] = posvecmaker_compress(posvec_all);
            [posvec,outersep] = posvecmaker_compress(P.posvec); %,'inds200','inds400',);
            
        end
        
        % Identify what to plot (excursions or manual) %%%%%%%%%%%%%%%%%%%%
        
        if ~isempty(plot_start)
            numplots = 1;
%             foundflag = 0;
%             for exx = 1:numexcur
%                 decbins = [ P.Tvecs{exx}(1:(end-1))   P.Tvecs{exx}(2:end) ] ;
%                 if any(logical(isExcluded( plot_start:plot_end , [decbins(1,1) decbins(end,2)] )))
%                     foundflag = 1;
%                     tvec = P.Tvecs{exx};
%                     break
%                 end
%             end
%             if ~foundflag
%                 disp('target period not found')
%                 keyboard
%             end
%             disp(sprintf('%s: day %d ep %d (time: %d - %d) (EXCUR: %d)',P.animalname(1:3),d,ep,numexcur,plot_start,plot_end,EXTYPE_toplot))            
        else
            numplots = length(P.Posts);
            disp(sprintf('%s: day %d ep %d (%d excurs) (EXCUR: %d)',P.animalname(1:3),d,ep,numplots,EXTYPE_toplot))
        end
               
    end
    

% Plotter of posterior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotter

    for PL = 1:numplots
        
        K = figure('units','normalized','outerposition',[.1 .06 .6 .8]);
        AX(1) = subplot(7,1,[2 3]);
        
        tvec = P.Tvecs{PL};
        decbins = [ tvec(1:(end-1))   tvec(2:end) ];
        plot_start = tvec(1) ;
        plot_end = tvec(end);
            a = lookup(plot_start,tvec);
            b = lookup(plot_end,tvec);
            aa = lookup(plot_start,postimevec - postimevec(1));
            bb = lookup(plot_end,postimevec - postimevec(1));
        
         %%% mua %%%   
         muatvec = plot_start:0.001:plot_end;
         muavec = zeros(1,length(muatvec));
         for tet = selected_tets
             filedata = loadparamsfile(daydir,d,tet);  % spike data from day directory
             inds_thresh = any(filedata.params(:,2:5) > P.spikethresh, 2) ;  %  Spikes with min amplitude
             spiketimes = filedata.params(inds_thresh,1)/10000 - epstart;
             onehist = histc(spiketimes,muatvec);
             onehist(end) = 0;
             muavec = muavec + onehist(:)';
         end
         muavec = muavec * 1000;   % in Hz
         msmuakern = 20;  % in ms
         muakern = gaussian(msmuakern,msmuakern*10);
         muavec = smoothvect(muavec,muakern);
         if 0
             figure;
             plot(tvec,muavec,'k-')
         end
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        posvec2 = posvec(aa:bb);
        
       % Process posterior into image
       post        = P.Posts{PL}(:,a:b);
       postfull    = P.Posts{PL}(:,:);
       
       % convert to RGB
       post_rgb = nan(size(post,1),size(post,2),3);  % x, y, (r,g,b)
       yvec = 1:size(post,1);
       Cinds = yvec <= CPbuff;                          % Center arm
       Linds = (yvec > CPbuff) & (yvec <= outersep);    % Left arm
       Rinds = ~Cinds & ~Linds;                         % Right arm
       % C arm
       %maxval = max(post(:));
       maxval = MAXVAL;
       post_rgb(Cinds,:,1) = (maxval - post(Cinds,:))/maxval;
       post_rgb(Cinds,:,2) = (maxval - post(Cinds,:))/maxval;
       post_rgb(Cinds,:,3) = (maxval - post(Cinds,:))/maxval;
       % L arm  (red)
       post_rgb(Linds,:,1) = 1;
       post_rgb(Linds,:,2) = (maxval - post(Linds,:))/maxval;
       post_rgb(Linds,:,3) = (maxval - post(Linds,:))/maxval;
       % R arm  (blue)
       post_rgb(Rinds,:,1) = (maxval - post(Rinds,:))/maxval;
       post_rgb(Rinds,:,2) = (maxval - post(Rinds,:))/maxval;
       post_rgb(Rinds,:,3) = 1;
       
       % saturation
       post_rgb(post_rgb < 0) = 0;
       posvec2(posvec2 <= 0) = nan; 
       

    % plot image of posterior %%%%%%%%%%%%%%%%%%%%%%%%%
    maxpos = size(post_rgb,1);
    image(tvec(a:b),1:maxpos,post_rgb); hold on
    colormap jet
    %colormap(flipud(gray))
    set(gca,'ydir','normal')
    set(gca,'xticklabel','')

    % Plot animal's position %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % lindist (all)  (green stripe)
    jumpinds = find(   diff([nan posvec_all(:)']) > 30  );  % where goes from center arm to R arm
    posvec_all_2 = posvec_all;
    posvec_all_2(jumpinds) = nan;
    plot(postimevec(aa:bb) - postimevec(1),posvec_all_2(aa:bb),'k-','linewidth',5);
    plot(postimevec(aa:bb) - postimevec(1),posvec_all_2(aa:bb),'-','Color',[40 200 40]/255,'linewidth',2); hold on
    % plot that single jumped sample!
    for jj = jumpinds(:)'
        linkinds = jj:(jj+1);
        plot(postimevec(linkinds) - postimevec(1),posvec_all(linkinds),'k-','linewidth',5);
        plot(postimevec(linkinds) - postimevec(1),posvec_all(linkinds),'-','Color',[40 200 40]/255,'linewidth',2); hold on
    end

    if 0
        % posvec (encoding positions, taken from P file) (stripe)
        % identify jumps and put a single nan for sake of plotting
        jumpinds = find(diff([nan posvec(:)']) > 30 );  % where goes from center arm to R arm
        posvec2 = posvec;
        posvec2(jumpinds) = nan;
        plot(postimevec(aa:bb) - postimevec(1),posvec2,'k-','linewidth',5); hold on
        plot(postimevec(aa:bb) - postimevec(1),posvec2,'-','Color',[40 200 40]/255,'linewidth',2); hold on
        % plot that single jumped sample!
        for jj = jumpinds(:)'
            linkinds = ((jj-1):(jj)) + aa;
            plot(postimevec(linkinds) - postimevec(1),posvec(linkinds - aa + 1),'k-','linewidth',5);
            plot(postimevec(linkinds) - postimevec(1),posvec(linkinds - aa + 1),'-','Color',[40 200 40]/255,'linewidth',2); hold on
        end
    end

    % Plot arm boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot([plot_start plot_end],[CPbuff CPbuff],'k--','linewidth',2)
    plot([plot_start plot_end],[outersep outersep],'k--','linewidth',2)
    % Plot 1 s grid
    for sc = 1:length(plot_start:plot_end)
        plot([plot_start+sc-1 plot_start+sc-1],[0 max(posvec_all)+1],'color',[.5 .5 .5]); hold on
    end
    set(gca,'fontsize',14)

    % Plot posterior traces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AX(2) = subplot(7,1,6);
    %plot(tvec(a:b),sum(postfull(Cinds,a:b),1),'k-','linewidth',2); hold on
    if 1
        plot(tvec(a:b),sum(postfull(Linds,a:b),1),'r-','linewidth',2); hold on
        plot(tvec(a:b),sum(postfull(Rinds,a:b),1),'b-','linewidth',2); hold on
    else
        smoothL = smoothvect(sum(postfull(Linds,a:b),1),muakern);
        smoothR = smoothvect(sum(postfull(Rinds,a:b),1),muakern);
        plot(tvec(a:b),smoothL,'r-','linewidth',2); hold on
        plot(tvec(a:b),smoothR,'b-','linewidth',2); hold on
    end
    ylim([0 1.1])  ;      set(gca,'xticklabel','')
    set(gca,'fontsize',14)

    % Plot LFP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AX(3) = subplot(7,1,1);
    if plot_DELTA_LFP
        lfp_a = lookup(plot_start,deltatimevec-epstart);
        lfp_b = lookup(plot_end,deltatimevec-epstart);
        plot(deltatimevec(lfp_a:lfp_b)-epstart,LFPSIGN *deltatrace(lfp_a:lfp_b),'-','color',[.75 .75 .9],'linewidth',2); hold on
    end
    if plot_LOWFREQ_LFP
        lfp_a = lookup(plot_start,lowfreqtimevec-epstart);
        lfp_b = lookup(plot_end,lowfreqtimevec-epstart);
        plot(lowfreqtimevec(lfp_a:lfp_b)-epstart,LFPSIGN *lowfreqtrace(lfp_a:lfp_b),'-','color',[.8 .8 .9],'linewidth',2); hold on
    end
    if 1 %plot_THETA_LFP
        lfp_a = lookup(plot_start,thetatimevec-epstart);
        lfp_b = lookup(plot_end,thetatimevec-epstart);
        plot(thetatimevec(lfp_a:lfp_b)-epstart,LFPSIGN * thetatrace(lfp_a:lfp_b),'k-','linewidth',2); hold on
        ylim([-LFPMAX LFPMAX])
    end
    if plot_LFP
        lfp_a = lookup(plot_start,eegtimevec-epstart);
        lfp_b = lookup(plot_end,eegtimevec-epstart);
        plot(eegtimevec(lfp_a:lfp_b)-epstart,LFPSIGN *eegtrace(lfp_a:lfp_b),'k-','linewidth',2); hold on
    end

    set(gca,'fontsize',14)
    set(gca,'xticklabel','')

    % Unit MUA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    AX(4) = subplot(7,1,4);
    area(muatvec,muavec,'facecolor','k');
    set(gca,'xticklabel','')
    set(gca,'fontsize',14)
    xlim([plot_start plot_end])

    % Plot Speed (linear + angular) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AX(5) = subplot(7,1,5);
    postimevec2 = postimevec - postimevec(1);
    pos_a = lookup(plot_start,postimevec2) - 1;
    pos_b = lookup(plot_end,postimevec2);
    %  speed
    area(postimevec2(pos_a:pos_b),veltrace(pos_a:pos_b),'facecolor',[.8 .8 .8],'linewidth',2);   hold on
    %  ang speed
    area(postimevec2(pos_a:pos_b),ANGVEL_SCALE * angvel(pos_a:pos_b),'facecolor',[250 120 0]/255,'edgecolor','none')
    set(gca,'fontsize',14)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    linkaxes(AX,'x');

    % plot title
    AX(3) = subplot(7,1,1);
    title(sprintf('%s day %d ep %d',animpref(1:3),d,ep),'fontweight','bold','fontsize',16)

    keyboard
    %pause
    %close all
    end
end
        
    
    if analyze
        
        excurperiods = [];

        if EXCURPERIODS == 1
            % Excursion + before CP
            excursions = loaddatastruct(animalinfo{2},animalinfo{3},'excursions',P.dayep(1));
            excurlist = excursions{d}{ep}.excurlist;
                inds = ismember(excurlist(:,3),[1 3 -11]); % outbound trajs (1,3) + outbound trackbacks
                excurlist = excurlist(inds,1:3);
            excurvec = list2vec(excurlist(:,[1 2]),postimevec);
            intersectvec = CPprevec & excurvec & movingvec;
            excurperiods = vec2list(intersectvec,postimevec);  
        elseif EXCURPERIODS == 2
%            %  vel4 + head out + before CP
%            headdir = linpos{d}{ep}.statematrix.segmentHeadDirection(:,1);   % head direction relative to center well -- values > 0 are outbound
%                 postimevec_nonnan = postimevec(~isnan(headdir));
%                 headdir_nonnan = headdir(~isnan(headdir));
%                 headdir2 = interp1(postimevec_nonnan,headdir_nonnan,postimevec,'linear');
%            outdirvec = logical(headdir2 >= 0);
%             excurperiods = vec2list(movingvec(:)' & outdirvec(:)',postimevec);
        end
        
        choiceperiodsdur = round(sum(excurperiods(:,2) - excurperiods(:,1)));
        disp(sprintf('%d sec of choiceperiods (%d)',choiceperiodsdur,EXCURPERIODS))
        
        choicevec = logical(list2vec(excurperiods - epstart,tvec - 0.001));
        
        Cdensity    = sum(postfull(Cinds,:),1);
        Ldensity    = sum(postfull(Linds,:),1);
        Rdensity    = sum(postfull(Rinds,:),1);
        tvec_cut    = tvec(choicevec);  % corresponding chopped up times;
        
        % XCorr Left vs. Right %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if plotter && plot_Xcorr_LR
            [Xcorr_LR,tlag] = xcorr(Rdensity(choicevec),Ldensity(choicevec));
            H = figure('units','normalized','outerposition',[.1 .6 .2 .2]);
            plot(tlag/1000,Xcorr_LR,'k-','linewidth',2);
            xlim([-0.5 0.5])
            set(gca,'fontsize',14)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Theta phase histogram of L vs. R density >> minphase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Ldensity; Rdensity;
            LRdensity = Ldensity + Rdensity;
                LRthreshinds = LRdensity >= LR_prop_thresh; 
                if ELIMINATE_LR_EXCLUSIVE_MINPHASE
                    disp('eliminating anything above .9 for minphase est')
                    LRthreshinds = LRthreshinds & (LRdensity < 1); 
                end
                
            if plotter && plot_LRprop_hist 
               H = figure('units','normalized','outerposition',[.2 .06 .2 .2]);
               hist(LRdensity(LRthreshinds(:) & choicevec(:)),20); axis tight
               title('LR density','fontweight','bold','fontsize',12)
            end
            
        % Determine min phase for theta (for LR density)    
        thetaphase;
        thetatimevec;
            tphs = thetaphase(lookup(tvec_cut,thetatimevec - epstart)); % theta phase of each density bin
            LRthreshinds2 = LRthreshinds(lookup(tvec_cut,tvec));
        nbins = 12;  %18;  % 36 bins is 10 deg bins, after Jezek 2011
            phasebins       = -pi:(2*pi/nbins):pi;
            phasebins_c     = phasebins(1:(end-1));       
        [~,I] = histc(tphs,phasebins);

        phasehist = zeros(size(phasebins_c));
        for bbb = 1:nbins
            phasehist(bbb) = sum(LRdensity( (I == bbb) & LRthreshinds2(:) )); 
        end
        
        if PHASECHOICE == 1
            [~,I2] = min(phasehist);  % min
            minphase = phasebins_c(I2);
        elseif PHASECHOICE == 2
            [~,I2] = max(phasehist);  %max
            minphase = phasebins_c(I2);    
        end
        
        % Theta phase histogram of LR density %%%%%%%%%%%%%%%%%%%%%%%%
        if plotter && plot_LRprop_phasehist
            H = figure('units','normalized','outerposition',[.1 .3 .2 .2]);
            B = bar(phasebins_c,phasehist,'histc'); hold on
            set(B,'facecolor','k')
            set(gca,'fontsize',14)
            axis tight
            plot([minphase minphase],[0 max(phasehist)],'r-','linewidth',3); hold on
            title(sprintf('min phase : %0.3f',minphase))
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% Identify theta bins (using minphase) %%%%%%%%%%%%%%%%%%%%%%%%%%%
        [tbins,~,~,~,~,~,tbblocks]  ...
               = thetabinner(thetaphase,thetatimevec,minphase,[epstart epend]);  
            % Filter out bins that don't occur during excurperiods (outbound excursion, and prior to CP)
            filterinds = logical(isExcluded(mean(tbins,2),excurperiods));
            tbins = tbins(filterinds,:);
                tbins_mean = mean(tbins,2);
                totalnumcyc = size(tbins,1);
                cycletimevec = 1:totalnumcyc;
       
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Get Spatial Index (SI) value (0: Left, 1: Right) for each theta bin
        SIvec_orig = nan(1,totalnumcyc);
        Altvec_orig = nan(1,totalnumcyc);
        
        for tb = 2:size(tbins,1)
            bin_a = lookup(tbins(tb,1)-epstart,tvec);  % start index of bin
            bin_b = lookup(tbins(tb,2)-epstart,tvec);  % end index of bin
            Lval = sum(Ldensity(bin_a:bin_b));
            Rval = sum(Rdensity(bin_a:bin_b));
            Cval = sum(Cdensity(bin_a:bin_b));
            LRval = Rval / (Lval + Rval);
            % Filter for at least minimum LR posterior density
            %%%binstarttime = tbins(tb,1)-epstart;
            %%%if binstarttime > 241.676 && binstarttime < 241.68
            %%%    keyboard
            %%%end
            if (Lval + Rval)/(Cval + Rval + Lval) < LR_prop_thresh
                %disp('low thresh')
                continue  
            else
                % Check to make sure cycle did not come after discontinuous theta bin gap
                if tbins(tb,1) == tbins(tb-1,2)   
                    % Store SI value
                    SIvec_orig(tb) = LRval;
                    if tb ~= 1 && ~isnan(SIvec_orig(tb-1))
                        % Calculate Alt speed
                       Altvec_orig(tb) = abs( LRval - SIvec_orig(tb-1) );   % alternation speed  
                    end
                else
                   continue 
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Plot SI in theta bins %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if plotter &&  ~omit_mainplot 
        figure(K);
        AX(8) = subplot(7,1,7);  % go back and plot on top of posteriors
        warning('off','all')
        linkaxes(AX,'x');
        warning('on','all')
        colormap(redblue)
        cmap = flipud(colormap);
        clrscale = (1:(size(colormap,1))) / size(colormap,1);
        
        for tbb = 1:size(tbins,1)
            bin_start = tbins(tbb,1) - epstart;
            bin_end =   tbins(tbb,2) - epstart;
            if all(isExcluded([bin_start bin_end],[plot_start plot_end]))
                if ~isnan(SIvec_orig(tbb)) && any(isExcluded([bin_start bin_end]+epstart,movingperiods))
                    siclr = cmap( lookup(SIvec_orig(tbb),clrscale) , : ) ;
                    % SI color patch
                    patch([bin_start bin_end bin_end bin_start],...
                          [0 0 1 1],siclr,'edgecolor','k','linewidth',2); hold on
                    % Alternation speed (dot)  
                    plot(mean([bin_start bin_end]),Altvec_orig(tbb),'k.','markersize',20); hold on
                end     
            end
        end
        % plot Alt speed trace
        plot(mean(tbins,2)-epstart,Altvec_orig,'k-','linewidth',2); hold on
        
        % plot threshold line
        plot([plot_start plot_end],[ALT_THRESHOLD ALT_THRESHOLD],'k--','linewidth',1)
        set(gca,'fontsize',14)
        ylim([-.1 1.1])
        
        % plot excursions
        for PL = 1:size(excurlist,1)
            if all(isExcluded(excurlist(PL,[1 2])-epstart,[plot_start plot_end]))
                excur_a = excurlist(PL,1) - epstart;
                excur_b = excurlist(PL,2) - epstart;
                patch([excur_a excur_b excur_b excur_a],[1.1 1.1 1.25 1.25],[.85 .85 .85],'edgecolor','none');
            end
        end
        
        set(gca,'xticklabel','')        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
        % Alternation analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % Initialize outputs
        clear out
        out.animalname         = animalinfo{3};
        out.dayep              = dayep;
        out.numperms           = NUMSAMP;
        out.divider1           = '*********************';
        out.altdist_ep         = nan(NUMSAMP,length(ALTDISTBINS)-1);
        out.numaltpers_ep      = nan(NUMSAMP,MAXALTS);
        out.altmeans_ep        = nan(1,NUMSAMP);
        out.altmedians_ep      = nan(1,NUMSAMP); 
        out.divider2           = '*********************';
        out.altperiods         = []; 
        out.tbins              = tbins;  % start and end times of each theta bin
        
        % I. Detect/Permute @ epoch level %%%%%%%%%%%%%%
        
        for PM = 1:NUMSAMP
            
            if mod(PM,100) == 0
                disp(['Permute #: ' num2str(PM)]) 
            end
            
            % Initialize permuted output
            SIvec_perm  = SIvec_orig;               % copy SI values
            Altvec_perm = Altvec_orig;  %
            
            % 1. Permute %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if PM > 1
                
%                 % Iterate through each excursion period
%                 for xc = 1:size(excurlist,1)
%                     
%                     % Identify thetabins with SI values in this excur period
%                     cycinds = logical((tbins_mean > excurlist(xc,1)) & ...
%                                       (tbins_mean < excurlist(xc,2)));
%                     validcycinds = cycinds(:)' & ~isnan(SIvec_orig);
%                     numcyc_xc = sum(validcycinds);
%                     cycinds = find(validcycinds);
%                     
%                     % Permute order
%                     reord = randperm(numcyc_xc);
%                     
%                     % Apply permutation to SI values
%                     SIvec_perm(cycinds) = SIvec_perm(cycinds(reord));
%                     
%                 end
                
                % Iterate through excursion type
                for XTYPE = [1 3 -11]   % 1: L, 2: R, -11: Trackback
                    exinds = excurlist(:,3) == XTYPE;
                    % identify all thetabins that fall in this excursion type
                    excurtype_pers = excurlist(exinds,[1 2]);
                    if ~isempty(excurtype_pers)
                        cycinds = find( isExcluded(tbins_mean,excurtype_pers) );
                        numcyc = length(cycinds);
                        
                        % Resampling method
                        if SAMPMETHOD == 1
                            % Permute order
                            sampord = randperm(numcyc);
                        elseif SAMPMETHOD == 2
                            % Bootstrap
                            sampord = ceil( rand(1,numcyc) * numcyc) ;
                        end
                        
                        % Apply sampling method to SI values
                        SIvec_perm(cycinds) = SIvec_perm(cycinds(sampord));
                    end
                end
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % 2. Calculate Alternation speed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for tb = 2:size(tbins,1)
                if ~isnan(Altvec_perm(tb))   % if it was in original, it can be calculated here
                    Altvec_perm(tb) = abs( SIvec_perm(tb) - SIvec_perm(tb-1) );   % alternation speed
                end
            end
            N = histc(Altvec_perm,ALTDISTBINS);
                N(end-1) = N(end-1) + N(end);
                N(end) = [];
            % install in output
            out.altdist_ep(PM,:)   = N(:);
            out.altmeans_ep(PM)    = nanmean(Altvec_perm);
            out.altmedians_ep(PM)  = nanmedian(Altvec_perm);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % 3. Detect Alternation periods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            altvec      = Altvec_perm > ALT_THRESHOLD;
            altpers     = vec2list(altvec,cycletimevec);
            
            % filter for discontiguous gaps
            for p = size(altpers,1):-1:1
                cyclevec = logical(list2vec(altpers(p,:),cycletimevec));
                cyclelist = tbins(cyclevec,:);
                % detect any gaps within the alt period
                if ~all(  cyclelist(1:(end-1),2) == cyclelist(2:end,1))
                    %disp('gap in alt period detected, deleting')
                    altpers(p,:) = [];
                    continue
                end
            end
            
            % Identify / store alternation periods
                % (for next analyses + plotting below)
            if PM == 1
                Altpers_orig        = altpers;                  % alternation periods in terms of cycle #s
                numaltperiods_orig  = size(Altpers_orig,1);     % total # of alternation periods
                Altpers_orig_xc     = nan(1,numaltperiods_orig); % excursion type of each alt period
                    for vv = 1:numaltperiods_orig
                        % mid time of alternation period
                        altperiods_midtime = tbins_mean(round(mean(altpers(vv,:),2)));                         
                        % excursion #
                        excurnum = find( excurlist(:,1) < altperiods_midtime & ...
                                         excurlist(:,2) > altperiods_midtime );
                        % excursion type
                        Altpers_orig_xc(vv) = excurlist(excurnum,3);
                    end
            end
            
            % Tabulate cycle durations
            altpers_dur = [altpers(:,2) - altpers(:,1)] + 1;     % duration, in theta cycles, of each candidate period
            for mm = 1:MAXALTS
                out.numaltpers_ep(PM,mm) = sum(altpers_dur >= mm);              % periods with at least 4 cycles
            end
            if PM == 1
               Altpers_dur_orig = altpers_dur; 
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
         % II. Detect/Permute @ single altperiod level %%%%%%%%%%%%%%
        
         % Preliminarily, obtain SIvals for each xcursion type for this epoch 
         clear cycinds_xc numcyc_xc sivals_xc
         cycinds_xc{1}      = find( isExcluded(tbins_mean,excurlist(excurlist(:,3) == 1,[1 2])) );  % Right
         cycinds_xc{3}      = find( isExcluded(tbins_mean,excurlist(excurlist(:,3) == 3,[1 2])) );  % Left 
         cycinds_xc{11}     = find( isExcluded(tbins_mean,excurlist(excurlist(:,3) == -11,[1 2])) ); % Trackback
             numcyc_xc{1}   = length(cycinds_xc{1});
             numcyc_xc{3}   = length(cycinds_xc{3});
             numcyc_xc{11}  = length(cycinds_xc{11});
                sivals_xc{1} = SIvec_orig(cycinds_xc{1});
                sivals_xc{3} = SIvec_orig(cycinds_xc{3});
                sivals_xc{11} = SIvec_orig(cycinds_xc{11});
             
         % initialize output
         out.altperiods = nan(numaltperiods_orig,7);  % [ time_a  time_b   tb_a   tb_b  numcyc  xcurtype  pvalue ]
         
         for ll = 1:numaltperiods_orig
             
            % Identify basic information for this alternation period
            tb_a = Altpers_orig(ll,1) -1;  % theta bin # of first cycle
            tb_b = Altpers_orig(ll,2);  % theta bin # of last cycle 
            time_a = tbins(tb_a,1);     % absolute time of start of first cycle
            time_b = tbins(tb_b,2);     % absolute time of end of last cycle
            numcyc = tb_b-tb_a+1;       % number of theta cycles
            xc = Altpers_orig_xc(ll);  % excursion type
            
            if xc == -11
                xc = 11;  % trackbacks >> can't index with negative value
            end
            
            % Identify indices of the surrounding continuous thetabin period ("block period") in which this occurred
            tbc_a = nan;  % "theta bin contig" first bin # of contig theta period (with decoded SI)
            tbc_b = nan;  % "theta bin contig" last bin # of contig theta period (with decoded SI)
            blockval = tbblocks(tb_a);  % dummy val is 1 or 2 -- designates a contiguous theta period
                % initialize
            tbc_a = tb_a;
            tbc_b = tb_b;
            % now find the beginning end by stepping
            while tbc_a > 0
               if tbblocks(tbc_a) == blockval && ~isnan(SIvec_orig(tbc_a)) % same contiguous theta block && SIvec is decoded
                   tbc_a = tbc_a - 1;
               else
                   break
               end
            end
            while tbc_b <= totalnumcyc
               if tbblocks(tbc_b) == blockval && ~isnan(SIvec_orig(tbc_b)) % same contiguous theta block && SIvec is decoded
                   tbc_b = tbc_b + 1;
               else
                   break
               end
            end
            blockdur = tbc_b - tbc_a + 1;  % duration of surrounding block period
            
            altdetect = nan(1,NUMSAMP);  % maximum alternation is 19
            for rrr = 1:NUMSAMP
                
                % Resampling method  (here just resample all the sivals)
                if SAMPMETHOD == 1
                    sivals_resamp = sivals_xc{xc}(randperm(numcyc_xc{xc}));  % Permute order
                elseif SAMPMETHOD == 2
                    sivals_resamp = sivals_xc{xc}(ceil( numcyc_xc{xc} * rand(1,numcyc_xc{xc}) )) ;  % Bootstrap
                end 
                
                % Take out the block period (from the beginning, since simplest)
                si_resamp = sivals_resamp(1:blockdur);
                
                % Detect alternations
                altspeed_resamp = abs(diff(si_resamp));
                altvec_resamp   = altspeed_resamp > ALT_THRESHOLD;
                altpers_resamp  = vec2list(altvec_resamp,1:length(altvec_resamp));
                altpers_dur_resamp = altpers_resamp(:,2) - altpers_resamp(:,1) + 1;
                
                % If produce at least one alternation during this period,
                % of at least duration of the actual alternation
                if any(altpers_dur_resamp >= numcyc)
                    altdetect(rrr) = 1;
                else
                    altdetect(rrr) = 0;
                end
            end
            % Calculate pvalue
            if sum(altdetect) ~= 0
                pval = sum(altdetect) / NUMSAMP;
                %%% testing
                %if pval < 0.001
                %    keyboard
                %end
            else
                pval = - 1 / NUMSAMP;  % negative indicates that the pval is --less than-- the absolute value
            end
            
            
            % install outputs
            out.altperiods(ll,1) = time_a;  % [ time_a  time_b   tb_a   tb_b  numcyc  xcurtype  pvalue ]
            out.altperiods(ll,2) = time_b; 
            out.altperiods(ll,3) = tb_a; 
            out.altperiods(ll,4) = tb_b; 
            out.altperiods(ll,5) = numcyc; 
            out.altperiods(ll,6) = xc;
            out.altperiods(ll,7) = pval; 
            
         end
         
         

         
        

       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % 2. Detect Alternation periods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            altvec      = Altvec_perm > ALT_THRESHOLD;
            altpers     = vec2list(altvec,cycletimevec);
            
            % filter for discontiguous gaps
            for p = size(altpers,1):-1:1
                cyclevec = logical(list2vec(altpers(p,:),cycletimevec));
                cyclelist = tbins(cyclevec,:);
                % detect any gaps within the alt period
                if ~all(  cyclelist(1:(end-1),2) == cyclelist(2:end,1))
                    %disp('gap in alt period detected, deleting')
                    altpers(p,:) = [];
                    continue
                end
            end
            
            % Store original Alternation periods (for plotting below)
            if PM == 1
                Altpers_orig = altpers;
            end
            
            % Tabulate cycle durations
            altpers_dur = [altpers(:,2) - altpers(:,1)] + 1;     % duration, in theta cycles, of each candidate period
            for mm = 1:MAXALTS
                out.numaltpers_ep(PM,mm) = sum(altpers_dur >= mm);              % periods with at least 4 cycles
            end
            if PM == 1
               Altpers_dur_orig = altpers_dur; 
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
     
        
        
        
        
        % 4. Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plot Alt speed hist %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if plotter && plot_altspeed_hist
                H = figure('units','normalized','outerposition',[.1 .8 .2 .2]);
                N = histc(Altvec_orig,ALTDISTBINS);
                N(end-1) = N(end-1) + N(end);
                B2 = bar(bincenterer(ALTDISTBINS),N(1:(end-1)),'histc');
                set(B2,'facecolor','k');
                title('Alt speed hist','fontsize',14,'fontweight','bold')
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Bar graph of alternation periods %%%%%%%%%%
            if plot_altspeed_counts
            H = figure('units','normalized','outerposition',[.9 .6 .2 .4]);
            BR = bar(1:MAXALTS,out.numaltpers_ep(1,:),'histc'); hold on
            set(BR,'facecolor',[255 230 108]/255,'edgecolor','k','linewidth',3)
            % plot mean of permuted
            for rr = 1:MAXALTS
                meanpermval = mean(out.numaltpers_ep(2:end,rr)) ;
                plot([rr-.3 rr+.3],[meanpermval meanpermval],'-','color',[.7 .7 .7],'linewidth',2)
            end
            title({'Alt speed histogram',num2str(out.numaltpers_ep(1,:))},'fontsize',14,'fontweight','bold')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plot Alternation periods on SI patches %%%%%%%%%%%%%%%%%%%
            if plotter && ~omit_mainplot
            figure(K);
            AX(8) = subplot(7,1,7);  % go back and plot on top of posteriors
            for ppp = 1:size(Altpers_orig,1)
                altper_start = tbins(Altpers_orig(ppp,1)-1,1);  % begin of first theta cycle in alt per
                altper_end   = tbins(Altpers_orig(ppp,2),2);  % end of last theta cycle in alt per
                if any(isExcluded([altper_start altper_end]-epstart,[plot_start plot_end]))
                    patch([altper_start altper_end altper_end altper_start] - epstart,...
                        [1.3 1.3 1.5 1.5],[255 230 108]/255 ,'edgecolor','k','linewidth',2)
                end
            end
                % plot p < 0.05 significant alternation periods
            for ap = 1:size(out.altperiods,1)
                pvalue_alt = out.altperiods(ap,end);
                if pvalue_alt < 0.1
                     altper_start = out.altperiods(ap,1);  % begin of first theta cycle in alt per
                     altper_end   = out.altperiods(ap,2);  % end of last theta cycle in alt per             
                     if any(isExcluded([altper_start altper_end]-epstart,[plot_start plot_end]))
                         patch([altper_start altper_end altper_end altper_start] - epstart,...
                             [1.3 1.3 1.5 1.5],'r','edgecolor','k','linewidth',2)
                     end
                     altstring = sprintf('Sig Alt period: %0.1f - %0.1f s (p = %d) N = %d',...
                                            altper_start-epstart,altper_end-epstart,pvalue_alt,NUMSAMP);
                     disp(altstring);
                end
            end
            ylim([-.1 1.6])
            linkaxes(AX,'x');
            set(gca,'fontsize',14)
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Plot null (permuted) distributions 
        
        % (Epoch - wide) plot altspeed distributions
%         if NUMPERMS > 1
%             H = figure('units','normalized','outerposition',[.2 .5 .2 .2]);
%             subplot(1,2,1)
%                 altdist_emp = OUT.altdist(1,:);
%                 altdist_null = mean(OUT.altdist(2:end,:));
%             plot(bincenterer(ALTDISTBINS),altdist_emp/sum(altdist_emp),'r-','linewidth',2); hold on
%             plot(bincenterer(ALTDISTBINS),altdist_null/sum(altdist_null),'k-','linewidth',2)
%             subplot(1,2,2)
%             plot(1:MAXALTS,OUT.numaltpers(1,:),'r-','linewidth',2); hold on
%             plot(1:MAXALTS,mean(OUT.numaltpers(2:end,:),1),'k-','linewidth',2)
%         end
        % (Epoch - wide) plot alt period distributions
        if NUMSAMP > 1
            if plot_altspeed_counts
            H = figure('units','normalized','outerposition',[.8 0 .2 .5]);
            for n = 1:MAXALTS
                subplot(MAXALTS,1,n); hold on
                NN = hist(out.numaltpers_ep(2:end,n),4); hold on
                hist(out.numaltpers_ep(2:end,n),4); hold on
                empval = out.numaltpers_ep(1,n);
                plot([empval empval],[0 max(NN)],'r-','linewidth',4)
                if n == 1
                   title('# alternations by epoch','fontsize',14,'fontweight','bold') 
                end
                xlabel('# alts','fontsize',12)
                ylabel({sprintf('NUMALT = %d',n),'# resamplings'},'fontsize',12)
                set(gca,'fontsize',12)
                axis tight
            end 
            end
        end
            
       
        
            % alternation mean + median
            if 0
                figure
                subplot(1,2,1)
                MM = hist(OUT.altmeans,20);
                hist(OUT.altmeans(2:end),20); hold on
                empval = OUT.altmeans(1);
                plot([empval empval],[0 max(MM)],'r-','linewidth',2)
                subplot(1,2,2)
                MM = hist(OUT.altmedians,20);
                hist(OUT.altmedians(2:end),20); hold on
                empval = OUT.altmedians(1);
                plot([empval empval],[0 max(MM)],'r-','linewidth',2)
            end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        
        NUMALT_DETECT = 3;
        alttimes_long = tbins(Altpers_orig(find(Altpers_dur_orig >= NUMALT_DETECT),:))-epstart;
            alttimes_long_durs = Altpers_dur_orig(Altpers_dur_orig >= NUMALT_DETECT);
            alttimes_long_sig = out.altperiods(Altpers_dur_orig >= NUMALT_DETECT,end) < P_THRESH;
            alttimes_long_pval = out.altperiods(Altpers_dur_orig >= NUMALT_DETECT,end);
            alttimes_long = [alttimes_long    alttimes_long_durs  alttimes_long_pval  alttimes_long_sig];
            
            % p-value scatter
            if plot_pvalue_scatter 
            figure;
            WIDTH = 0.2;
            for tt = 1:length(alttimes_long_pval)
                pvalu = alttimes_long_pval(tt);
                scatter( WIDTH * rand * [1 1] - WIDTH/2,[pvalu pvalu],500,'k','.'); hold on
            end
            xlim([-1 1])
            title(sprintf('p-values of alts >= %d',NUMALT_DETECT),'fontsize',14,'fontweight','bold') 
            end
            end
            
        if ~isempty(alttimes_long)
            disp('long alt times:');
            alttimes_long
        else
            disp('no long alt times')
        end
        
    end
    
   
    
  
        % Place map plot: 2D map %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if plot_placemap
            
            F2 = figure('units','normalized','outerposition',[.7 .2 .2 .4]);
            % plot all positions
            plot(xpos,ypos,'.','markersize',5,'color',[.6 .6 .6]); hold on
            
            % plot CP positions
            %plot(xpos(CPvec),ypos(CPvec),'.','markersize',5,'color',[.4 .8 .4]); hold on
            
            % plot L and R positions
            plot(xpos(Lvec & ~CPvec),ypos(Lvec & ~CPvec),'.','markersize',5,'color',[0 0 1]); hold on
            plot(xpos(Rvec & ~CPvec),ypos(Rvec & ~CPvec),'.','markersize',5,'color',[1 0 0]); hold on
            
            % plot traversal positions
            %plot(xpos(pos_a:pos_b),...
            %    ypos(pos_a:pos_b),'k-','linewidth',5);   %[.65 .65 .65]
            plot(xpos(pos_a:pos_b),...
                ypos(pos_a:pos_b),'-','color',[40 200 40]/255,'linewidth',4);   %[.65 .65 .65]

            %scatter(xpos(pos_a:pos_b),...
            %        ypos(pos_a:pos_b),30,...
            %        'linewidth',1,'markeredgecolor','k','markerfacecolor',[40 200 40]/255);   %[.65 .65 .65]

            % plot head direction arrows
            for ii = pos_a:NUMSAMP_HEADPOS:pos_b   % every 250 ms
                xdir = 5*cos(hdir(ii));
                ydir = 5*sin(hdir(ii));
                % plot "arrow line"
                plot([xpos(ii)  xpos(ii)+xdir],...
                    [ypos(ii)  ypos(ii)+ydir],'linewidth',1.5,'color','k');
                % plot point at base of arrow
                scatter(xpos(ii),ypos(ii),50,'k','.');
            end
            for ii = pos_a:NUMSAMP_HEADPOS:pos_b   % every 250 ms
                xdir = 5*cos(hdir(ii));
                ydir = 5*sin(hdir(ii));                
                % plot point at tip of arrow
                plot(xpos(ii)+xdir,ypos(ii)+ydir,'.','markersize',15,'color','w');
                plot(xpos(ii)+xdir,ypos(ii)+ydir,'o','markersize',6,'color','k','linewidth',1.5);
            end            
            
            lowspeedflag = 0;
            
            % plot p < 0.01 significant alternation periods
            for ap = 1:size(out.altperiods,1)
                if out.altperiods(ap,end) < 0.05
                    altper_start = out.altperiods(ap,1);  % begin of first theta cycle in alt per
                     altper_end   = out.altperiods(ap,2);  % end of last theta cycle in alt per             
                     if ~any(isExcluded([altper_start altper_end]-epstart,[plot_start plot_end]))
                         continue
                     end
                     disp(sprintf('Plotting sig alt on placemap (%0.1f to %0.1f s)',altper_start-epstart,altper_end-epstart));
                    pos_alt_a = lookup(altper_start,postimevec);
                    pos_alt_b = lookup(altper_end,postimevec);                     
                    
                    % detect any low speed
                    lowspeedflag = any(veltrace(pos_alt_a:pos_alt_b) < 4);
                    
                    % plot alternation (enlarged black points)
                    %plot(xpos(pos_alt_a:pos_alt_b),...
                        %ypos(pos_alt_a:pos_alt_b),'.','markersize',15,'color','k');
                    
                    for ii = pos_a:NUMSAMP_HEADPOS:pos_b   % every 250 ms
                        if any(isExcluded(ii,[pos_alt_a pos_alt_b]))  % check if this position is in alternation period
                            xdir = 5*cos(hdir(ii));
                            ydir = 5*sin(hdir(ii));
                            % plot point at base of arrow
                            scatter(xpos(ii),ypos(ii),50,'k','.');
                            % plot "arrow line"
                            plot([xpos(ii)  xpos(ii)+xdir],...
                                [ypos(ii)  ypos(ii)+ydir],'linewidth',1.5,'color','k');
                        end
                    end
                end
            end
                    % second go around to get circles on top of plot
            for ap = 1:size(out.altperiods,1)
                if out.altperiods(ap,end) < 0.05       
                    altper_start = out.altperiods(ap,1);  % begin of first theta cycle in alt per
                     altper_end   = out.altperiods(ap,2);  % end of last theta cycle in alt per             
                     if ~any(isExcluded([altper_start altper_end]-epstart,[plot_start plot_end]))
                         continue
                     end
                     disp(sprintf('Plotting sig alt on placemap (%0.1f to %0.1f s)',altper_start-epstart,altper_end-epstart));
                    pos_alt_a = lookup(altper_start,postimevec);
                    pos_alt_b = lookup(altper_end,postimevec);                        
                    for ii = pos_a:NUMSAMP_HEADPOS:pos_b   % every 250 ms
                        if any(isExcluded(ii,[pos_alt_a pos_alt_b]))  % check if this position is in alternation period   
                            xdir = 5*cos(hdir(ii));
                            ydir = 5*sin(hdir(ii));                            
                            % plot circle at tip of arrow
                            plot(xpos(ii)+xdir,ypos(ii)+ydir,'.','markersize',15,'color',[255 255 20]/255);
                            plot(xpos(ii)+xdir,ypos(ii)+ydir,'o','markersize',6,'color','k','linewidth',1.5);
                            %uistack(ppp,'top');
                        end
                    end
                    
                end
            end
            
            
            
            
            % formatting
            set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[])
            set(gca,'color','w','box','off')
            %set(findobj(gcf, 'type','axes'), 'Visible','off')
            axis tight
            axis square
            
            title(sprintf('%s day %d ep %d (lowspeedflag: %d)',animalinfo{3}(1:3),d,ep,lowspeedflag),...
                  'fontsize',14,'fontweight','bold')
            %plotcount = plotcount + 1;
            
        end    
    
    
    
    
    break
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
% if plot_continuous_old
% 
%         plot_hg = 0;
%         tracekern_ms = [1000];
%         MSBIN = 100;      
%     
%     % Ripple and WG power traces
%     out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'rippletrace', d, ep);
%     riptrace = zscorer(out{d}{ep}{1}.powertrace);
%     riptrace_timevec = out{d}{ep}{1}.eegtimesvec_ref;
%     out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'fastgammatrace', d, ep);
%     wgtrace = zscorer(out{d}{ep}{2}.powertrace);
%     wgtrace_timevec = out{d}{ep}{2}.eegtimesvec_ref;
%     out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'lowgammatrace', d, ep);
%     lgtrace_ca1 = zscorer(out{d}{ep}{1}.powertrace);  % ca1
%     lgtrace_ca3 = zscorer(out{d}{ep}{2}.powertrace);  % ca3 dg
%     lgtrace_timevec = out{d}{ep}{2}.eegtimesvec_ref;
%     if plot_hg
%         out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'highgammatrace', d, ep);
%         hgtrace = zscorer(out{d}{ep}{1}.powertrace);
%     end
%     
%     if ~isempty(tracekern_ms)
%         sig_bins = 1500 * tracekern_ms/1000;
%         gauskern = gaussian(sig_bins,8*sig_bins);
%         wgtrace = smoothvect(wgtrace,gauskern);
%         lgtrace_ca1 = smoothvect(lgtrace_ca1,gauskern);
%         lgtrace_ca3 = smoothvect(lgtrace_ca3,gauskern);
%         if plot_hg
%             hgtrace = smoothvect(hgtrace,gauskern);
%         end
%     end
% 
%     selected_tets = D{1}.selected_tets;
%         maxtet = max(selected_tets);
%     spikethresh = D{1}.spikethresh;
%     %xvecms = D{1}.xvecms;
%     xvecms = round(epstart)*1000:1:round(epend)*1000 - epstart;
%     trajpos = D{1}.trajpos;
%     linposbins = D{1}.linposbins{1};
%     centerarmmax = D{1}.centerarmmax;
%     rightarmmax = D{1}.rightarmmax;
%     
%     figure;
%     
%     % Unit regional patches
%     ax(1) = subplot(6,1,[1 2 3]);
%     celllist = {};
%     % iterate
%     for cc = 1:maxcell
%         tet = adtc_list(cc,3);
%         cellnum = adtc_list(cc,4);
%         region = cellinfo{d}{ep}{tet}{cellnum}.area;
%             [~,~,~,bgclr] = regionfunc(region);
%         type = cellinfo{d}{ep}{tet}{cellnum}.type;
%         celllist = [celllist ; num2str(cc) type region];
%         disp(type)
%         if ~strcmp(type,'principal')
%             bgclr = [1 1 1];
%         end
%         
%         % plot region (CA1, CA3) color patch for region
%         patch([starttime starttime endtime endtime],[cc-0.5 cc+0.5 cc+0.5 cc-0.5],bgclr,'edgecolor','none');
%             hold on
%        
%     end        
%     
%     % First plot ripple patches and wg traces
%     if 0
%         for xx = 2  %:2
%             ymax = 20;
%             if xx == 1
%                 ax(1) = subplot(6,1,[1 2 3]); 
%                 ylevel = maxtet;
%             elseif xx == 2
%                 % obtain all clustered units this epoch
%                 an = find(strcmp(animal_toplot,animals_order));
%                 [adtc_list] = clusteredunits(spikes,cellinfo,an,d,ep,[],1);
%                 maxcell = size(adtc_list,1);
%                 ax(1) = subplot(6,1,[1 2 3]); 
%                 ylevel = maxcell;
%             end
%             % plot 2 SD SWRs that occur in the widnow %%%%%%%%%
%             plotwinvec = logical(list2vec([starttime endtime] + epstart,ripout{d}{ep}.time))';
%             consvec_rip2_toplot = consvec_rip2 & plotwinvec;
%             rip_toplot = vec2list(consvec_rip2_toplot,ripout{d}{ep}.time) - epstart;
%             numrip_toplot = size(rip_toplot,1);
%             % plot all 2 SD ripples in the raster window
%             for rp = 1:numrip_toplot
%                 patch([rip_toplot(rp,1) rip_toplot(rp,2) rip_toplot(rp,2) rip_toplot(rp,1)] ,...
%                     [0 0 (ylevel + ymax) (ylevel + ymax)],...
%                     [1 .9 .9],'edgecolor','none','Parent',ax(1) ); hold on
%             end
%             % plot WG and SWR power traces at bottom%%%%%%%%%%%%%%%%
%             ind1 = lookup(starttime,riptrace_timevec - epstart) ;
%             ind2 = lookup(endtime,riptrace_timevec - epstart) ;
%             ind3 = lookup(starttime,wgtrace_timevec - epstart) ;
%             ind4 = lookup(endtime,wgtrace_timevec - epstart) ;
%             plot((riptrace_timevec(ind1:ind2) - epstart), 4 * -riptrace(ind1:ind2)/2 + ylevel + 10,'-','Color',[.8 .75 .75],'linewidth',2,'Parent',ax(1) ); hold on
%             plot((wgtrace_timevec(ind3:ind4) - epstart), 8 * -wgtrace(ind3:ind4)/2 + ylevel + 10,'-','Color',[0 .3 .9],'linewidth',2,'Parent',ax(1) ); hold on
%             plot((lgtrace_timevec(ind1:ind2) - epstart), 8 * -lgtrace_ca1(ind1:ind2)/2 + ylevel + 10,'-','Color',[1 .1 .1],'linewidth',2,'Parent',ax(1) ); hold on
%             %plot((lgtrace_timevec(ind1:ind2) - epstart), 4 * -lgtrace_ca3(ind1:ind2)/2 + ylevel + 10,'-','Color',[1 .4 .4],'linewidth',2,'Parent',ax(1) ); hold on
%             
%             
%             xlim([starttime endtime])
%             ylim([0  (ylevel + ymax)])
%             set(gca,'ydir','reverse')
%             set(gca,'fontsize',12)
%             set(gca,'tickdir','out');
%             xlabel('Time (ms)') ;   ylabel('unit #');
%         end    
%     end
%     
%  
%     
%     
%                 % spike raster iterate
%             for cc = 1:maxcell
%                 tet = adtc_list(cc,3);
%                 cellnum = adtc_list(cc,4);
%                 if ~isempty(spikes{d}{ep}{tet}{cellnum}.data)
%                     spiketimes = spikes{d}{ep}{tet}{cellnum}.data(:,1);
%                     spiketimes = spiketimes(logical(isExcluded(spiketimes,[starttime endtime]+epstart))) ;  % spikes within window
%                     numunitspk = length(spiketimes);
%                     if numunitspk > 0
%                     for ll = 1:numunitspk
%                         plot([spiketimes(ll) spiketimes(ll)] - epstart,[cc-0.48 cc+0.48],...
%                             'linestyle','-','Color',[0 0 0],'LineWidth',1,'Parent',ax(1) ); hold on
%                     end
%                     end
%                 end
%             end
% 
%     
% 
%     % link axes
%     linkaxes(ax,'x');
% 
% end
% 
% 
