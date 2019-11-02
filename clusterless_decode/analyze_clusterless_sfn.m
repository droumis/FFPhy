datadir = '/opt/data50/kkay/__Decode';

calculate = 0;
plot_pvalues = 1;
plot_behavior = 1;
report_prevalence = 1;

    PHASECHOICE = 1;
    
    if plot_pvalues || plot_behavior
       animals_toplot = {'Bond','Frank','Dave','Government'}; %{'Bond','Frank','Dave','Government'};  %'Bond','Frank','Dave',
    end
    if calculate
        animals_tocalc = {'Bond','Frank','Dave','Government'}; %{'Egypt','Corriander'}; % 'Bond','Frank','Dave','Government','Dave','Corriander','Miles','Eight'};
            manual_days = [];
        MOVINGPERIODS = 2; % 1: vel4, 2: nonimmobile05 (half sec)
        EXCURPERIODS = 1;  % 1: Excursions + <CPbuff // 2: Head out + <CP buff
        % Spatial index
        LR_prop_thresh = 0.1;   % minimum posterior density that is either in L or R (versus C)        
        % Alternation detection parameters
        ALT_THRESHOLD       = 0.5;
        NUMSAMP            = 10000 + 1;  % # of resampling methods
            SAMPMETHOD = 1;   % 1: permutation, 2: bootstramp
        ALTDISTBINS         = 0:0.05:1;
        MAXALTS             = 3;    
        P_THRESH            = 0.001;
    end 
    
    % Study figures %%%%%%%%%%%%
    plot_Xcorr_LR = 0;
    plot_LRprop_hist = 0;
    plot_LRprop_phasehist = 0;
    plot_altspeed_hist = 0;

    if plot_behavior || plot_pvalues
       BEHAVE_P_THRESH = 0.05;
       LENGTH_TOPLOT = 2:3;
    end
    if plot_pvalues
        pbins = 0:.01:.35;
    end
    
    
    % Data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if calculate
        
        cd(datadir)
        
        % loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander',...
            'Iceland','Justice','Kapital','Laplace','Eight','Ten','Conley','Miles'};
        
        clear out
        
        for aa = 1:length(animals_tocalc)
            
            animalname = animals_tocalc{aa};
            animalpref = animalname(1:3);
            animalinfo = animaldef(animalname);
            daydir = getdaydir(animalname);
            an = find(strcmp(animalname,animals_order));
            task = loaddatastruct(animalinfo{2},animalinfo{3},'task');
            if isempty(manual_days)
                days = 1:length(task);
            else
                days = manual_days;
            end
            dayeps = [];
            for dd = days
                if ~isempty(task{dd})
                    eps = wtrackeps(task,dd);
                    dayeps = [dayeps ; dd * ones(length(eps),1)  eps(:)];
                end
            end
            
            % Initialize outputs
            out.date               = date;
            out.animalname         = animalinfo{3};
            out.dayeps             = dayeps;
            out.numperms           = NUMSAMP;
            out.ALTDISTBINS        = ALTDISTBINS;
            out.divider1           = '%%%%%%%%%%%%%%%%%%%%%%%%%%%%';
            out.altdist_ep         = {}; %nan(NUMSAMP,length(ALTDISTBINS)-1);
            out.numaltpers_ep      = {}; %nan(NUMSAMP,MAXALTS);
            out.altmeans_ep        = {}; %nan(1,NUMSAMP);
            out.altmedians_ep      = {}; %nan(1,NUMSAMP);
            out.divider2           = '%%%%%%%%%%%%%%%%%%%%%%%%%%%%';
            out.altperiods         = {};
            out.descript           = '1: day 2: ep 3: time_a, 4: time_b, 5: tb_a, 6: tb_b, 7: numcyc, 8: xcurtype, 9-11: speed, accel, angspeed, 12: pval';

            for de = 1:size(dayeps,1)
            
                d = dayeps(de,1);
                ep = dayeps(de,2);
                
                disp(sprintf('%s: d %d ep %d',animalinfo{3}(1:3),d,ep))
                
                % load Decode file
                if ~exist('P','var') || isempty(P) || ~strcmp(P.animalname,animalname) || ...
                        ~all(P.dayep == [d ep])
                    % Decode file %%%%%%%%%%%%%%
                    filename = sprintf('%s_%d_%d_fullepoch_1_outpos.mat',animalpref,d,ep);
                    cd(datadir);
                    filedir = dir(filename);
                    if isempty(filedir)
                        disp(sprintf('not finding %s, skipping',filename))
                        continue
                    end
                    load(filedir.name,'P');                    
                    tvec = (1:size(P.posteriors{1},2))/1000 - 0.001;   % 1 ms timevec
                    % Positional %%%%%%%%%%%%%%%%%%%
                    pos = loaddatastruct(animalinfo{2},animalinfo{3},'pos',P.dayep(1));
                    linpos = loaddatastruct(animalinfo{2},animalinfo{3},'linpos',P.dayep(1));
                    lindist = linpos{d}{ep}.statematrix.lindist;
                    CPbuff = choicepointer(linpos{d}{ep});
                    CPprevec = lindist < CPbuff;
                    postimevec = pos{d}{ep}.data(:,1);
                    veltrace = pos{d}{ep}.data(:,5);
                    epstart = postimevec(1);
                    epend = postimevec(end);
                    headang = loaddatastruct(animalinfo{2}, animalinfo{3}, 'headang', d);
                    angveltrace = abs(headang{d}{ep}.angvel)  * 180 / pi;    % convert to degrees
                    % tets, clust spikes
                    spikes = loaddatastruct(animalinfo{2},animalinfo{3},'spikes',d);
                    cellinfo = loaddatastruct(animalinfo{2},animalinfo{3},'cellinfo');
                    [adtc_list] = clusteredunits(spikes,cellinfo,an,d,ep,[],1);
                    maxcell = size(adtc_list,1);
                    selected_tets = unique(P.selected_tets);
                    numtets = length(selected_tets);
                    
                    % LFP %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    thetabins = loaddatastruct(animalinfo{2},animalinfo{3},'thetabins',d);
                    thetatet = thetabins{d}.thetatet;
                    if strcmp(animalpref,'Dav') && d == 3
                        disp('Dav thetatet = 1, d 3 ep 2');
                        thetatet = 1;
                        %elseif strcmp(animpref,'Bond') && d == 6
                        %    disp('Bon thetatet = 30, d 6 ep 4');
                        %    thetatet = 30;
                    end
                    theta = loadeegstruct(animalinfo{2},animalinfo{3},'theta',d,ep,thetatet);
                        thetatrace = theta{d}{ep}{thetatet}.data(:,1);
                        thetaphase = double(theta{d}{ep}{thetatet}.data(:,2))/10000;  % hilbert phase
                        thetatimevec = geteegtimes(theta{d}{ep}{thetatet});                    
                end
            
                % Identify Moving periods
                timefilterscript
                if MOVINGPERIODS == 1
                    [dummy,~] = evaluatetimefilter(animalinfo{2},animalinfo{3},{ vel4 },[d ep]);
                elseif MOVINGPERIODS == 2
                    [dummy,~] = evaluatetimefilter(animalinfo{2},animalinfo{3},{ nonimmobile05 },[d ep]);
                elseif MOVINGPERIODS == 3
                    [dummy,~] = evaluatetimefilter(animalinfo{2},animalinfo{3},{ nonimmobile1 },[d ep]);
                end
                movingperiods = dummy{d}{ep};   %
                movingvec = logical(list2vec(movingperiods,postimevec));
               
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Process posterior into image
                validxbins = P.validxbins{1};
                postfull    = P.posteriors{1}(validxbins,:);
                
                % Process position data
                [posvec_all] = posvecmaker(lindist,linpos{d}{ep});
                [posvec_all,outersep] = posvecmaker_compress(posvec_all);

                % convert to RGB
                yvec = 1:size(postfull,1);
                Cinds = yvec <= CPbuff;                          % Center arm
                Linds = (yvec > CPbuff) & (yvec <= outersep);    % Left arm
                Rinds = ~Cinds & ~Linds;                         % Right arm
                
                excurperiods = [];
                
                % Identify excursions
                if EXCURPERIODS == 1
                    % Excursion + before CP
                    excursions = loaddatastruct(animalinfo{2},animalinfo{3},'excursions',P.dayep(1));
                    excurlist = excursions{d}{ep}.excurlist;
                    inds = ismember(excurlist(:,3),[1 3 -11 23 32]); % Outbound trajs (1,3) + Center well trackbacks + R>>L and L>>R
                    excurlist = excurlist(inds,1:3);
                    excurvec = list2vec(excurlist(:,[1 2]),postimevec);
                    intersectvec = CPprevec & excurvec & movingvec;
                    excurperiods = vec2list(intersectvec,postimevec);
                else
                    keyboard
                end
                
                choiceperiodsdur = round(sum(excurperiods(:,2) - excurperiods(:,1)));
                disp(sprintf('%d sec of excurperiods (%d)',choiceperiodsdur,EXCURPERIODS))
                
                choicevec = logical(list2vec(excurperiods - postimevec(1),tvec));
                
                Cdensity    = sum(postfull(Cinds,:),1);
                Ldensity    = sum(postfull(Linds,:),1);
                Rdensity    = sum(postfull(Rinds,:),1);
                tvec_cut    = tvec(choicevec);  % corresponding chopped up times;
                
                % XCorr Left vs. Right %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if plot_Xcorr_LR
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
                
                if plot_LRprop_hist
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
                if plot_LRprop_phasehist
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
                    if (Lval + Rval)/(Cval + Rval + Lval) < LR_prop_thresh
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
                
            % Alternation analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % I. Detect/Permute @ epoch level %%%%%%%%%%%%%%
            
            for PM = 1:NUMSAMP
                
                if mod(PM,100) == 0
                    disp(['Permute #: ' num2str(PM)])
                end
                
                % 1. Permute %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Initialize permuted output
                SIvec_perm  = SIvec_orig;   % copy SI values
                Altvec_perm = Altvec_orig;  %
                
                % Iterate through excursion type
                if PM > 1
                    % must identify non-NaN cycles prior to permutation
                    nonnancyc = ~isnan(SIvec_orig) ;                    
                    for XTYPE = [1 3 -11]   % 1: L, 2: R, -11: Trackback
                        exinds = excurlist(:,3) == XTYPE;
                        % identify all thetabins that fall in this excursion type
                        excurtype_pers = excurlist(exinds,[1 2]);
                        if ~isempty(excurtype_pers)
                            excurinds = logical(isExcluded(tbins_mean,excurtype_pers));
                            cycinds = find(excurinds & nonnancyc(:));
%                            cycinds = find( isExcluded(tbins_mean,excurtype_pers) );
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
                            SIvec_perm(cycinds) = SIvec_orig(cycinds(sampord));
                        end
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                % 2. Calculate Alternation speed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                for tb = 2:size(tbins,1)
                    if ~isnan(Altvec_orig(tb))   % if it was in original, it can be calculated here
                        Altvec_perm(tb) = abs( SIvec_perm(tb) - SIvec_perm(tb-1) );   % alternation speed
                    end
                end
                N = histc(Altvec_perm,ALTDISTBINS);
                N(end-1) = N(end-1) + N(end);
                N(end) = [];
                % install in output
                out.altdist_ep{de}(PM,:)   = N(:);
                out.altmeans_ep{de}(PM)    = nanmean(Altvec_perm);
                out.altmedians_ep{de}(PM)  = nanmedian(Altvec_perm);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % 3. Detect Alternation periods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                altvec      = Altvec_perm > ALT_THRESHOLD;
                altpers     = vec2list(altvec,cycletimevec);
                    % eliminate alternation periods lasting only 1 cycle
                	altpers(  (altpers(:,2) - altpers(:,1)) == 0 , : ) = [];
                    
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
                        if isempty(excurnum)
                            Altpers_orig_xc(vv) = nan;
                            disp('***************** excurtype not found, setting to nan')
                        else
                            Altpers_orig_xc(vv) = excurlist(excurnum,3);
                        end
                    end
                end
                
                % Tabulate cycle durations
                altpers_dur = [altpers(:,2) - altpers(:,1)] + 1;     % duration, in theta cycles, of each candidate period
                for mm = 1:MAXALTS
                    out.numaltpers_ep{de}(PM,mm) = sum(altpers_dur >= mm);              % periods with at least 4 cycles
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
         out.altperiods{de} = nan(numaltperiods_orig,10);  % [ time_a  time_b   tb_a   tb_b  numcyc  xcurtype  speed  accel  angspeed  pvalue ]
         
         % Iterate through each alternation period
         for ll = 1:numaltperiods_orig
             
            % Identify basic information for this alternation period
            tb_a = Altpers_orig(ll,1);  % theta bin # of first cycle
            tb_b = Altpers_orig(ll,2);  % theta bin # of last cycle 
            time_a = tbins(tb_a,1);     % absolute time of start of first cycle
            time_b = tbins(tb_b,2);     % absolute time of end of last cycle
            numcyc = tb_b-tb_a+1;       % number of theta cycles
            xc = Altpers_orig_xc(ll);  % excursion type
            
            if xc == -11
                xc = 11;  % trackbacks >> can't index with negative value
            end
            if isnan(xc)
                disp('xcursion type cant be found, skipping')
                continue
            end
            
            % Identify behavior (speed, accel, angular speed)
                % take mean
            pos_a           = lookup(time_a,postimevec);
            pos_b           = lookup(time_b,postimevec);
            speed           = mean(veltrace(pos_a:pos_b));         % cm/s
            accel           = mean(diff(veltrace(pos_a:pos_b)));  % cm/s^2
            angspeed        = mean(angveltrace(pos_a:pos_b));    % deg / s
            lowspeedflag    = any(veltrace(pos_a:pos_b) < 4);
            
            % Permute each event %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % First need to identify indices of the surrounding continuous thetabin period ("block period") in which this occurred
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
            blockdur = tbc_b - tbc_a;  % duration of surrounding block period
            
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
            else
                pval = - 1 / NUMSAMP;  % negative indicates that the pval is --less than-- the absolute value
            end
            %%%%%%%%%%%%%  end of permute code  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % install outputs
            out.altperiods{de}(ll,1) = d;  % [ day ep  time_a  time_b   tb_a   tb_b  numcyc  xcurtype  speed accel angspeed pvalue ]
            out.altperiods{de}(ll,2) = ep;             
            out.altperiods{de}(ll,3) = time_a - epstart;   
            out.altperiods{de}(ll,4) = time_b - epstart; 
            out.altperiods{de}(ll,5) = tb_a; 
            out.altperiods{de}(ll,6) = tb_b; 
            out.altperiods{de}(ll,7) = numcyc; 
            out.altperiods{de}(ll,8) = xc;
            out.altperiods{de}(ll,9) = speed;
            out.altperiods{de}(ll,10) = accel;
            out.altperiods{de}(ll,11) = angspeed; 
            out.altperiods{de}(ll,12) = lowspeedflag; 
            out.altperiods{de}(ll,13) = pval; 
            
         end
         
     
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % 2. Detect Alternation periods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            altvec      = Altvec_perm > ALT_THRESHOLD;
            altpers     = vec2list(altvec,cycletimevec);
                    % eliminate alternation periods lasting only 1 cycle
                	altpers(  (altpers(:,2) - altpers(:,1)) == 0 , : ) = [];
            
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
                out.numaltpers_ep{de}(PM,mm) = sum(altpers_dur >= mm);              % periods with at least 4 cycles
            end
            if PM == 1
               Altpers_dur_orig = altpers_dur; 
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %  Printed reportage of alt periods >= 3 cycles
            NUM_LONGALT_REPORT = 3;
            alttimes_long           = tbins(Altpers_orig(find(Altpers_dur_orig >= NUM_LONGALT_REPORT),:))-epstart;
            alttimes_long_durs      = Altpers_dur_orig(Altpers_dur_orig >= NUM_LONGALT_REPORT);
            alttimes_long_sig       = out.altperiods{de}(Altpers_dur_orig >= NUM_LONGALT_REPORT,end) < P_THRESH;
            alttimes_long_pval      = out.altperiods{de}(Altpers_dur_orig >= NUM_LONGALT_REPORT,end);
            alttimes_long           = [alttimes_long    alttimes_long_durs  alttimes_long_pval  alttimes_long_sig];
            
            if ~isempty(alttimes_long)
                disp('long alt times:');
                alttimes_long
            else
                disp('no long alt times')
            end
            
            
            
            end
            
            cd(datadir);
            filename = sprintf('Analysis_%s',animalpref);
            save(filename,'out','-v7.3');
            
        end
        
        
    end
    
    
    
    if plot_pvalues
        
        % Collect p-values
        pvalues_collect = {};
        for L = 2:3  % 2 is exactly 2 alternations, 3 is >= 3
            for M = 1:2   % 1: non-low speed,  2: includes low-speed
                pvalues_collect{L}{M} = [];
                cd(datadir)
                for aa = 1:length(animals_toplot)
                    animalname = animals_toplot{aa};
                        animalinfo = animaldef(animalname);
                    filename = dir(sprintf('Analysis_%s.mat',animalname(1:3)));
                    load(filename.name,'out')
                    for de = 1:length(out.altperiods)
                        if ~isempty(out.altperiods{de})
                            if L == 2
                                cycinds = out.altperiods{de}(:,7) == 2;
                            elseif L == 3
                                cycinds = out.altperiods{de}(:,7) >= 3;
                            end
                            if M == 1
                                speedinds = out.altperiods{de}(:,12) == 0;
                            elseif M == 2
                                speedinds = out.altperiods{de}(:,12) == 1;
                            end
                            pvalues = abs(out.altperiods{de}(cycinds & speedinds,end));
                            pvalues_collect{L}{M} = [pvalues_collect{L}{M} ; pvalues(:)];
                        end
                    end
                end
            end
        end
        
        %%% p-value scatter
%           K = figure('units','normalized','outerposition',[.1 .06 .3 .7]); hold on
%         WIDTH = 0.2;
%         for L = LENGTH_TOPLOT
%             for M = 1:2  % M   1: moving only  2: includes low speed
%                 if L == 2
%                     ptclr = 'k';
%                     ptstyle = 'o';
%                     ptsize = 50;
%                 elseif L == 3
%                     ptclr = 'k';
%                     ptsize = 500;
%                     ptstyle = '.';
%                 end
%                 if M == 1
%                     ptclr = 'k';
%                 elseif M == 2
%                     ptclr = 'r';
%                 end
%                 for tt = 1:length(pvalues_collect{L}{M})
%                     onep = pvalues_collect{L}{M}(tt);
%                     scatter( WIDTH * rand * [1 1] - WIDTH/2,[onep onep],ptsize,ptclr,ptstyle,'linewidth',2); hold on
%                 end
%             end
%         end
%         xlim([-1 1])
%         animstring = mat2str(cell2mat(animals_toplot));
%         title(sprintf('%s: p vals, alt pers',animstring),'fontsize',14,'fontweight','bold')
        
          % p-value histogram
          if 0
          K = figure('units','normalized','outerposition',[.1 .06 .3 .5]); hold on
        for L = LENGTH_TOPLOT
            subplot(2,1,L-1);
            NNN = [];
            for M = 1:2  % M   1: moving only  2: includes low speed
                if M == 1
                    faceclr = 'k';
                elseif M == 2
                    faceclr = 'r';
                end
                
                NN = histc(pvalues_collect{L}{M},pbins);
                NN(end-1) = NN(end-1) + NN(end);
                NNN = [NNN ; NN(:)'];
            end
            
            B = bar(bincenterer(pbins),NNN(:,1:(end-1))','stacked'); hold on
                %set(B,'facecolor',faceclr);
            xlabel('P value','fontsize',14)
            ylabel('# of alternation periods','fontsize',14)
            
            set(gca,'fontsize',14)

            maxvalue = 10*ceil(max(sum(NNN,1))/10);
            ylim([0 maxvalue])
            set(gca,'ytick',0:20:maxvalue)
            
            % p < 0.05 line
            plot([0.05 0.05],[0 maxvalue],'k--','linewidth',2)
            
            % title
            animstring = mat2str(cell2mat(animals_toplot));
            title(sprintf('%s: P values, alt pers',animstring),'fontsize',14,'fontweight','bold')   
            

        
        end
        
          end
          % p-value histogram ->> log version
          K = figure('units','normalized','outerposition',[.1 .06 .3 .5]); hold on
        for L = LENGTH_TOPLOT
            subplot(2,1,L-1);
            NNN = [];
            collectvals = [];
            for M = 1:2  % M   1: moving only  2: includes low speed
                if M == 1
                    faceclr = 'k';
                elseif M == 2
                    faceclr = 'r';
                end
                
                logbins = 0:.1:(4+.1);
                neglogvals = -log10(abs(pvalues_collect{L}{M}));
                
                NN = histc(neglogvals,logbins);
                NN(end-1) = NN(end-1) + NN(end);
                NNN = [NNN ; NN(:)'];
                
                collectvals = [collectvals ; neglogvals(:)];
            end
            
            B = bar(bincenterer(logbins),NNN(:,1:(end-1))','stacked'); hold on
                %set(B,'facecolor',faceclr);
            xlabel('-log10(P value)','fontsize',14)
            ylabel('# of alternation periods','fontsize',14)
            
            set(gca,'fontsize',14)

            xlim([0 max(logbins)])
            maxvalue = 10*ceil(max(sum(NNN,1))/10);
            if maxvalue > 20
                set(gca,'ytick',0:10:maxvalue)
            else
                maxvalue = maxvalue - 5;
                set(gca,'ytick',0:5:maxvalue)
            end
            ylim([0 maxvalue])
            % p < 0.05 line
            plot(-log10([0.05 0.05]),[0 maxvalue],'k--','linewidth',2)
            
            % title
            animstring = mat2str(cell2mat(animals_toplot));
            title(sprintf('%s:  P values, alt pers',animstring),'fontsize',14,'fontweight','bold')   
            

        
        end
                
        
        
            
    end
    
    
    if plot_behavior
        
        % Collect behavior values
        behaviorvals = {};
        dayeptimes = {};
        for M = 1:2
            for L = 2:3  % 2 is exactly 2 alternations, 3 is >= 3
                behaviorvals{L}{M} = [];
                dayeptimes{L}{M} = [];
                for aa = 1:length(animals_toplot)
                    animalname = animals_toplot{aa};
                    filename = dir(sprintf('Analysis_%s.mat',animalname(1:3)));
                    load(filename.name,'out')
                    for de = 1:length(out.altperiods)
                        if ~isempty(out.altperiods{de})
                            pinds = out.altperiods{de}(:,end) < BEHAVE_P_THRESH;
                            if L == 2
                                cycinds = out.altperiods{de}(:,7) == 2;
                            elseif L == 3
                                cycinds = out.altperiods{de}(:,7) >= 3;
                            end
                            if M == 1
                                speedinds = out.altperiods{de}(:,12) == 0;
                            elseif M == 2
                                speedinds = out.altperiods{de}(:,12) == 1;
                            end
                            % behavior values
                            behav_vals = out.altperiods{de}(pinds & cycinds & speedinds,[ 9 10 11] );
                            behaviorvals{L}{M} = [behaviorvals{L}{M} ; behav_vals];
                            % day ep times
                            dayeptimes{L}{M} = [dayeptimes{L}{M} ; out.altperiods{de}(pinds & cycinds & speedinds,[1 2 3 4 7])];
                        end
                    end
                end
            end
        end
        
        % behavior scatter plot
        K = figure('units','normalized','outerposition',[.5 .06 .3 .5]); hold on
        for L = LENGTH_TOPLOT 
            subplot(2,1,L-1);
            for M = 1:2
                
                if M == 1   % M = 1: moving only
                    if L == 2
                        ptsize = 500;
                        ptstyle = '.';
                        ptclr = 'k'; %[.85 .85 .85];
                    elseif L == 3
                        ptsize = 500;
                        ptstyle = '.';
                        ptclr = 'k';
                    end
                elseif M == 2  % M = 2: low speed too
                    if L == 2
                        ptsize = 50;
                        ptstyle = 'o';
                        ptclr = 'k'; %[.85 .85 .85];
                    elseif L == 3
                        ptsize = 50;
                        ptstyle = 'o';
                        ptclr = 'k';
                    end
                end
                for tt = 1:size(behaviorvals{L}{M},1)
                    one_b = behaviorvals{L}{M}(tt,:);
                    scatter(one_b(1),one_b(3),ptsize,ptclr,ptstyle,'linewidth',2); hold on
                    %scatter3(one_b(1),one_b(2),one_b(3),ptsize,ptclr,ptstyle,'linewidth',2); hold on
                end
            end
            xlabel('Speed (cm/s)','fontsize',14)
            ylabel('Angular speed (deg/s)','fontsize',14)
            set(gca,'ytick',0:2:10)
            %zlabel('accel','fontweight','bold','fontsize',12)
            grid on
            xlim([0 60])
            ylim([0 10])
            set(gca,'fontsize',14)
        end
        
        animstring = mat2str(cell2mat(animals_toplot));
        subplot(2,1,1);
        title(sprintf('%s: Behav vals, alt pers',animstring),'fontsize',14,'fontweight','bold')
        
    end
    
    dayeptimes{3}{1}
    
    
if report_prevalence
    for ee = 1:length(out.altperiods)
       d = out.altperiods{ee}(1,1);
       ep = out.altperiods{ee}(1,2);
       excursions = loaddatastruct(animalinfo{2},animalinfo{3},'excursions',[d ep]);
                    excurlist = excursions{d}{ep}.excurlist;
                    inds = ismember(excurlist(:,3),[1 3 -11]); % outbound trajs (1,3) + outbound trackbacks
                    excurlist = excurlist(inds,1:3);       
                    numoutbound = size(excurlist,1);
                    
       numsigaltpers = sum(out.altperiods{ee}(:,13) < 0.05);
       
       disp(sprintf('# out excur: %d, # sigaltpers: %d -- %0.2f altper / excur',numoutbound,numsigaltpers,numsigaltpers/numoutbound))
       
    end
end
    
    break
    
   
  
