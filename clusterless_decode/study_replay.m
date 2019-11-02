% Cluster decoding on W-track -- psuedo 1D "armdists" (Wu--Foster-2014 Fig 1)

datadir_pseudo = '/opt/data13/kkay/Superclustdecode_data';
datadir_dir = '/opt/data13/kkay/Superclustdecode_dir_data';

calc_cand = 1;   % identifies candidate events and saves in P variable (Wu-Foster-2014)
    if calc_cand
       targetfile = ''; %'Bon_4_2.mat';  % leave blank if want to process all files 
       plotcheck = 0;
            ploteeg = 0;
            if ploteeg
                tracekern_ms = [];
                clear tracespec
                tracespec{1} = {'rippletrace',1,[.7 .7 .7]} ;   % tracefile name + TFgroup
                tracespec{2} = {'fastgammatrace',2,[.4 .6 1]} ;
                tracespec{3} = {'lowgammatrace',1,[.8 .1 .1]}   ;
            end
       BIAS_CENTER = 1;         % 1: center well as reference / 2: center junct as reference
            CANDTHRESH = 2;     % spiking SD threshold above mean
            MINBIN = 3;         % minimum number of decoding bins with significant decode 
            MINCOVERAGE = 30;   % cm
            MINCORR = 0.3;      % weighted correlation minimum
            MINDIR = 0;       
            MINDIRBIAS = 0;
    end

plot_events = 0;   % calls plot_replay_events
    if plot_events
       animal_toplot = 'Bond';
            animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
       eventname = 'ripples';
       dayep = [4 6];

       loaddata = 1;
            % other
       plot_remoteW = 0;
       tracekern_ms = 10;       
       omitraster = 0;
       omitposterior = 0;
       plot_nonlocalbins = 1;    
    end

if calc_cand
    
        cd(datadir_dir)
        if isempty(targetfile)
            filenames = dir('Fra_8*.mat');
        else
            filenames = dir(targetfile);
        end     
        
    % load each file 
    for rr = 1:length(filenames)
       
        % load dir file   'D'
        cd(datadir_dir)
        if isempty(targetfile)
            filenames = dir('Fra_8*.mat');
        else
            filenames = dir(targetfile);
        end        
        filename = filenames(rr).name;
        load(filename,'D')
        
        % load pseudo file  'P'
        cd(datadir_pseudo)
        filename = filenames(rr).name;
        load(filename,'P')

        animalname = P.animalname;
            animalinfo = animaldef(animalname);
        d = P.dayep(1);
        ep = P.dayep(2);
            remoteep = P.remoteep;
        eps = [ep remoteep];
        
        disp(sprintf('%s %d %d replay',animalname(1:3),d,ep))

        spikes = loaddatastruct(animalinfo{2},animalinfo{3},'spikes',d);
        pos = loaddatastruct(animalinfo{2},animalinfo{3},'pos',d);
            epstart = pos{d}{ep}.data(1,1);
            epend = pos{d}{ep}.data(end,1);
            postimevec = pos{d}{ep}.data(:,1);
            mstimevec = epstart:.001:epend;
        timefilterscript    % 
        outperiods = evaluatetimefilter(animalinfo{2},animalinfo{3},{ripvel},[d ep]);  
            ripvelperiods = outperiods{d}{ep};
           
            
        sdf = [];  % spike density function
            
        % Determine candidate events for each encoding-decoding epoch (LOCAL and REMOTE)

        replayoutput = cell(1,2); 
        replay_nondur_count = [0 0];
        cand_count = [nan nan];
        
        for nn = 1:length(eps)  % iterate through LOCAL and REMOTE
            
            candstarts = [];
            candends = [];
            canddurations = [];
            numcandevents = nan;
            arminds = cell(1,3);
            
            % basic values
            binsize_t = P.binvec_c(2) - P.binvec_c(1);  % time bin size
            binsize_p = P.linposbins{nn}(2) - P.linposbins{nn}(1);  % position bin size            
            
            % first identify linpos bins for each arm
            arminds = cell(1,3);  % 1: center, 2: right, 3: left
                lastcenterbin = P.lastcenterbin(nn);
                lastrightbin = P.lastrightbin(nn);
            arminds{1} = 1:lastcenterbin ;                            % center inds
            arminds{2} = (lastcenterbin+1):lastrightbin ;           % right inds
            arminds{3} = (lastrightbin+1):length(P.linposbins{nn}) ;  % left inds
            
            % also identify threshold probability
            numposbins = length(P.linposbins{nn});
            minprob = 5 * (1/numposbins) ;
            
            % get place cell spike trains
            detc = P.detc{nn};
            numcells = size(detc,1);
            spk = [];
            for cc = 1:numcells
                tet = detc(cc,3);
                cellnum = detc(cc,4);
                if ~isempty(spikes{d}{ep}{tet}{cellnum}.data)
                    spk = [spk ; spikes{d}{ep}{tet}{cellnum}.data(:,1)]; % spike times
                end
            end
            [sdf,~] = histc(spk,mstimevec);
                gauskern = gaussian(15,15 * 8);  % 15 ms smoothing
            sdf = smoothvect(sdf,gauskern);
            
            if 0
                % check spiking function
               figure; plot(mstimevec/1000,sdf);
            end
            baseline = mean(sdf);
            thresh = CANDTHRESH * std(sdf);
            mindur = 0;
            tmpevent = extractevents(sdf, thresh, baseline, 0, mindur, 0)';
                startind = tmpevent(:,1);
                endind = tmpevent(:,2);
                    candstarts = epstart + startind / 1000;
                    candends = epstart + endind / 1000;
            
            % Filter for events that occur at low speeds (< 4 cm/s)
            inds = isExcluded(candstarts,ripvelperiods);
            candstarts(~inds) = [];
            candends(~inds) = [];
            canddurations = candends - candstarts;
                numcandevents = length(canddurations);
                disp(sprintf('%d events',numcandevents));
                cand_count(nn) = numcandevents;
                
            % Now, determine whether candidate events are REPLAY events   
            replayoutput{nn} = [];
            replay_nondur_count(nn) = 0;
            for c = 1:numcandevents

                replayfound = 0;  
                
                % decode bin inds
                a = lookup(candstarts(c),P.binvec_c);
                b = lookup(candends(c),P.binvec_c);

                    numdecbins = b - a + 1;
                    binvec = 1:numdecbins;
 
                    % now obtain "MAP" (maximum a priori) function
                        
                        MAP = [];
                        MAPPOS = [];
                        
                        POST = P.posteriors{nn}(:,a:b);  % posterior probability matrix for the event
                        MAP = max(POST,[],1);  % maximum probability in each decoding time bin
                        MAPPOS = max(POST,[],2);
                        if 1
                            % smooth over time
                            kern = gaussian(1, 1 * 8);           % smooth MAP w/ Gaussian of size 1 time bin (15 ms)
                            MAP = smoothvect(MAP,kern);
                            % smooth over position
                            kern = gaussian(2, 2 * 8);  % smooth MAP-pos w/ Gaussian of size 2 position bins bin (2 cm)
                            MAPPOS = smoothvect(MAPPOS,kern);
                        else
                            % don't smooth
                            MAP = map;
                        end

                    % Re-define candidate MAP event from the peak of MAP to where
                    % the MAP crosses below threshold
                    [~,maxbinind] = max(MAP);
                    threshinds = MAP > minprob;
                        if all(threshinds == 0)
                            continue
                        end
                    threshlist = vec2list(threshinds,binvec);
                    cand_a = []; cand_b = [];
                    for zz = 1:size(threshlist,1)
                        if isExcluded(maxbinind,threshlist(zz,:)) == 1
                            cand_a = binvec(threshlist(zz,1));
                            cand_b = binvec(threshlist(zz,2));
                            break
                        end
                    end
                    
                    if isempty(cand_a)
                        continue
                    end
                    
                        % cut the posterior and map
                    POST = POST(:,cand_a:cand_b);  % filter the posterior image to obtain the image for the MAP event
                    MAP = MAP(cand_a:cand_b);
                        % re-define indices
                    a1 = a + cand_a - 1 ;   
                    b1 = a + cand_b - 1 ;
                        stepvec = a1:b1;
                            
                    % Now determine 3 properties of the candidate MAP event:
                    % i. duration   (# of time bins in event)
                    % ii. arm coverage
                    % iii. weighted correlation
                    
                    % i. duration
                    numcandbins = size(POST,2);
                        binvec2 = 1:numcandbins;  % shorter bin vector
                        
                    % skip if not sufficient duration
                    if numcandbins < MINBIN
                        continue
                    end
                    
                    % ii. arm coverage
                    armcoverage = binsize_p * sum(max(POST,[],2) > minprob);
                    
                    % skip if not sufficient arm coverage
                    if armcoverage < MINCOVERAGE
                        continue
                    end
                    
                    % iii. weighted correlation
                    [ wcorr, bestcorr, besttraj ] = weightedcorr(POST',arminds);
                    
                    if bestcorr < MINCORR
                        continue
                    end
                    
                    replay_nondur_count = replay_nondur_count + 1;

                    %%%%%%%%%%%%%%%%%%%%% DIRECTIONAL REPLAY DETECTION %%%%%%%%%
                    
                    % (Review plot check)
                    if plotcheck
                        H = figure('units','normalized','outerposition',[0 0 .24 .94]);
                        if 0
                        % Plot image
                        subplot(9,1,[1 2 3])
                        linposbins = P.linposbins{nn} ;
                        imagesc(binvec2,linposbins,POST); hold on
                        ylabel('Linearized position','FontSize',12);
                        colormap(gray); colormap(flipud(colormap))
                        caxis([0 0.1]);
                        ylim([-.01 max(linposbins)]);
                        set(gca,'tickdir','out');
                        set(gca,'fontsize',12)
                        set(gca,'ydir','normal')
                        %xlabel('Bin #')
                        xlim([binvec2(1)-0.5 binvec2(end)+0.5])
                        % plot lines demarcating positions of arms
                        % between Center and Right
                        plot([0 binvec2(end)+0.5],[P.lastcenterbin(nn) P.lastcenterbin(nn)],'--','Color','k','linewidth',2)
                        % between Right and Left
                        plot([0 binvec2(end)+0.5],[P.lastrightbin(nn) P.lastrightbin(nn)],'--','Color','k','linewidth',2)
                        
                        % plot animals position
                        if nn == 1
                            j = lookup(candstarts(c),postimevec);
                            k = lookup(candends(c),postimevec);
                            winpos = P.armdists( j:k );  % armdists_cut (linear) positions over the course of the ripple
                            plot( (postimevec(j:k)-candstarts(c)) / binsize_t,winpos,'x','linewidth',6,'Color',[.8 .8 .8]); hold on
                        end
                        end
                        
                        if 0
                        subplot(9,1,[1 2 3])
                        string1 = sprintf('%s day %d ep %d',animalname(1:3),d,ep);
                        string2 = sprintf('Event time: %0.1f sec (CR: %0.1f CL: %0.1f RL:%0.1f)',candstarts(c)-epstart,wcorr);
                        title({string1,string2},'fontsize',14,'fontweight','bold')
                        end
                        
                        % Plot MAP below
                        subplot(9,1,4);
                        plot([binvec2(1) binvec2(end)],[minprob minprob],'-','Color',[.4 1 .4],'linewidth',3); hold on
                        plot(binvec2,MAP,'b-','linewidth',5); hold on
                        xlim([binvec2(1)-0.5 binvec2(end)+0.5])
                        ylabel('MAP Function','FontSize',12);
                        set(gca,'fontsize',12)                        
                        
                        %pause
                        %close all 
                        
                    end
                    
                    % Calculate directionality of each time bin of the
                    % candidate event ('dirtrace')
                        % ingredients
                    D; cand_a; cand_b;
                        numpos = D.lastbin(nn);
                        offset = numpos;
                        % initialize output
                    dirtrace = nan(1,numcandbins);  % initialize output    
                    for tbin = 1:numcandbins
                        binnum = stepvec(tbin);
                        num = 0; den = 0;
                        for pbin = 1:numpos
                            num = num + abs( D.posteriors{nn}(pbin,binnum) - D.posteriors{nn}(pbin+offset,binnum) );
                            den = den + D.posteriors{nn}(pbin,binnum) + D.posteriors{nn}(pbin+offset,binnum);
                        end
                        dirtrace(tbin) = num / den;
                    end
                    if plotcheck
                        % plot check
                        subplot(9,1,5);
                        plot([binvec2(1) binvec2(end)],[0.3 0.3],'-','Color',[.4 1 .4],'linewidth',3);         hold on                
                        plot(binvec2,dirtrace,'k-','linewidth',5); hold on
                        plot(binvec2,dirtrace,'ok','markersize',10,'linewidth',4); 
                        xlim([binvec2(1)-0.5 binvec2(end)+0.5])
                        ylabel('Directionality','FontSize',12);
                        set(gca,'fontsize',12)                        
                        
                    end
                    
                    % Calculate directional "bias"
                    biastrace = nan(1,numcandbins);
                    for tbin = 1:numcandbins
                        binnum = stepvec(tbin);
                        num = 0; den = 0;
                        for pbin = 1:numpos
                            % At the point, check if we're referencing
                            % center junct (Wu--Foster-2014) or home well
                                % in Wu-style, bias > 0 : facing toward center junction
                                %              bias < 0 : facing away from center junction
                                % in Frank-style, bias > 0 : facing OUTBOUND (away from home well)
                                %                 bias < 0 : facing INBOUND (toward home well)
                            if BIAS_CENTER
                                if pbin <= D.lastcenterbin(nn)
                                    B_SIGN = +1;  % no need to reverse on the home arm
                                else
                                    B_SIGN = -1;  % on the outer arms, reverse 
                                end
                            end
                            num = num + B_SIGN * (D.posteriors{nn}(pbin,binnum) - D.posteriors{nn}(pbin+offset,binnum)) ;
                            den = den + D.posteriors{nn}(pbin,binnum) + D.posteriors{nn}(pbin+offset,binnum);
                        end
                        biastrace(tbin) = num / den;
                    end                    
                    if plotcheck
                        % plot check
                        subplot(9,1,6);
                            plot([binvec2(1) binvec2(end)],[0.3 0.3],'-','Color',[.4 1 .4],'linewidth',3); hold on
                            plot([binvec2(1) binvec2(end)],[0 0],'k--','linewidth',1.5);
                            plot([binvec2(1) binvec2(end)],[-0.3 -0.3],'-','Color',[.4 1 .4],'linewidth',3);
                        plot(binvec2,biastrace,'-','Color',[.3 .3 .3],'linewidth',5); hold on
                        plot(binvec2,biastrace,'o','Color',[.3 .3 .3],'markersize',10,'linewidth',4);
                        ylabel('Dir Bias','FontSize',12);
                        set(gca,'fontsize',12)
                        xlim([binvec2(1)-0.5 binvec2(end)+0.5])
                        ylim([-1 1])
                    end
                    
                    
                    if plotcheck
                        % plot check, directional posterior
                        if nn == 1
                            linposbins = D.linposbins{nn} ;
                        elseif nn == 2
                            % quick fix
                            linposbins = 1:D.lastbin(nn);
                        end
                        
                        % if you want to use the center-well as the reference,
                        %       then you need to transpose between the positional islands in D.posteriors
                        if BIAS_CENTER
                            
                            BLUE = nan(numpos,numcandbins);
                            RED = nan(numpos,numcandbins);
                            for pbin = 1:numpos
                                if pbin <= D.lastcenterbin(nn)
                                    % in center arm
                                    BLUE(pbin,:) = D.posteriors{nn}(pbin+offset,a1:b1);    % BLUE is Frank-style INBOUND posteriors (facing toward center well)
                                    RED(pbin,:) = D.posteriors{nn}(pbin,a1:b1);            % RED is Frank-style OUTBOUND posteriors (facing away from center well)
                                else
                                    % in outer arm
                                    BLUE(pbin,:) = D.posteriors{nn}(pbin,a1:b1);           %  "  same as above
                                    RED(pbin,:) = D.posteriors{nn}(pbin+offset,a1:b1);
                                end
                            end
                        else
                            % here is w/ center junction as reference (Wu--Foster)
                            BLUE = D.posteriors{nn}(p1:p2,a1:b1);              % BLUE is Wu-style INBOUND posteriors (facing towards center junct)
                            RED  = D.posteriors{nn}((p2+1):(p2*2),a1:b1);      % RED is Wu-style OUTBOUND posteriors (facing away from center junct)
                        end
                        
                        % now plot
                        scalemax = max([BLUE(:) ; RED(:)]);
                        BLUE = BLUE/scalemax;  % normalize images
                        RED = RED/scalemax;
                       
                        BRimage = ones(size(BLUE));  % initialize
                        BRimage = repmat(BRimage,[1 1 3]);
                        % Blue-red inversion procedure
                        for xx = 1:size(BRimage,1)
                            for yy = 1:size(BRimage,2)
                                redval = RED(xx,yy);
                                blueval = BLUE(xx,yy);
                                BRimage(xx,yy,1) = BRimage(xx,yy,1) - blueval;
                                BRimage(xx,yy,2) = BRimage(xx,yy,2) - redval - blueval;
                                BRimage(xx,yy,3) = BRimage(xx,yy,3) - redval;
                            end
                        end
                        % fix for rounding errors
                        BRimage(BRimage > 1) = 1;
                        BRimage(BRimage < 0) = 0;
                        
                        
                        %BRimage = im2uint8(BRimage * 255);
                        
                        %dirimage = cat(3,RED,GREEN,BLUE);
                        subplot(9,1,[1 2 3]);
                        image(BRimage); hold on
                        
                        %imagesc(binvec2,linposbins,POST); hold on
                        ylabel('Linearized position','FontSize',12);
                        colormap(gray); colormap(flipud(colormap))
                        caxis([0 scalemax]);
                        ylim([-.01 max(P.linposbins{nn})]);
                        set(gca,'tickdir','out');
                        set(gca,'fontsize',12)
                        set(gca,'ydir','normal')
                        %xlabel('Bin #')
                        xlim([binvec2(1)-0.5 binvec2(end)+0.5])
                        % plot lines demarcating positions of arms
                        % between Center and Right
                        plot([0 binvec2(end)+0.5],[P.lastcenterbin(nn) P.lastcenterbin(nn)],'--','Color','k','linewidth',2)
                        % between Right and Left
                        plot([0 binvec2(end)+0.5],[P.lastrightbin(nn) P.lastrightbin(nn)],'--','Color','k','linewidth',2)
                        
                        % plot animals position
                        if nn == 1
                            j = lookup(candstarts(c),postimevec);
                            k = lookup(candends(c),postimevec);
                            winpos = P.armdists( j:k );  % armdists_cut (linear) positions over the course of the ripple
                            plot( (postimevec(j:k)-candstarts(c)) / binsize_t,winpos,'x','linewidth',6,'Color',[.8 .8 .8]); hold on
                        end
                        
                        if 1
                            %keyboard
                        else
                            pause
                            close all
                            continue
                        end
                    end
                        
                    % Identify continuous segments (of at least min # of time bins)
                    %    of the same sign of BIAS (+ is outbound, - is inbound) within the candidate event           
                    biasperiods{1} = []; % outbound
                    biasperiods{2} = []; % inbound
                        % first fill in gaps in posterior (no clustered spiking) with neighboring
                        % bins' bias signs, if both the neighboring bins (on either side) are the same
                    nanbins = find(isnan(biastrace));
                    biastrace2 = biastrace;
                    for nb = nanbins
                        if nb ~= 1 && nb ~= numcandbins
                            if sign(biastrace(nb-1)) == sign(biastrace(nb+1))
                                neighborsign = sign(biastrace(nb+1));
                                biastrace2(nb) = .001 * neighborsign;  % use tiny value as a marker, .001
                            end
                        end
                    end
                        % next, identify "bias periods" = continuous segments of the same sign
                    biasperiods{1} = vec2list(biastrace2 > 0,1:numcandbins);
                    biasperiods{2} = vec2list(biastrace2 <= 0,1:numcandbins);

                    % Now iterate through each candidate bias period
                    %       each candidate bias period will be evaluated for replay
                    for btype = 1:2
                       for bp = size(biasperiods{btype},1):-1:1
                           
                           for passnum = 1:2  % passnum == 1 is the entire directional period
                                              % passnum == 2 is, if the entire candidate bias period period fails replay criteria,
                                              % a subset of its times will be attempted
                                              %     (from peak
                                              %     directionality to times
                                              %     above the threshold)
                            if passnum == 1
                                bin1 = biasperiods{btype}(bp,1);
                                bin2 = biasperiods{btype}(bp,2);
                            elseif passnum == 2
                                % redefine replay from a subset of times within the directional period, from peak
                                [~,peakbin] = max(biastrace(bin1:bin2));
                                bin1star = []; bin2star = [];
                                biaspers = vec2list(biastrace > MINDIRBIAS,1:numcandbins) ;
                                for ww = 1:size(biaspers,1)
                                   peakbin2 = (peakbin - 1 + bin1);
                                   % identify the bias period that contains the peak bias time bin
                                   if any(peakbin2 == biaspers(ww,:)) || ...
                                           (peakbin2 > biaspers(ww,1) && peakbin < biaspers(ww,2))
                                        bin1star = biaspers(ww,1);
                                        bin2star = biaspers(ww,2);
                                        %disp('replay boundaries re-defined')
                                        break
                                   end
                                end
                                if ~isempty(bin1star)
                                    bin1 = bin1star;
                                    bin2 = bin2star;
                                else
                                   break 
                                end
                            end
                            
                            % now calculate replay statistics!
                            
                            % i. duration   (# of time bins in event)
                            dur = (bin2 - bin1 + 1) ;
                            if dur < MINBIN;
                                %disp('MINBIN HALT!')
                                break
                            end
                            
                            % ii. arm coverage (references to non-directional posterior)
                            cvg = binsize_p * sum(max(POST(:,bin1:bin2),[],2) > minprob);
                            if 0
                                if cvg < MINCOVERAGE
                                    disp('directional MINCOVERAGE not long enough')
                                    break
                                end
                            end
                            
                            % iii. weighted correlation (references non-directional posterior)
                            [ wcorr, bestcorr, besttraj, besttrajstring] = weightedcorr(POST(:,bin1:bin2)',arminds);
                            if 0
                                if bestcorr < MINCORR
                                    disp('directional corr not high enough')
                                    if passnum == 1
                                        continue
                                    elseif passnum == 2
                                        break
                                    end
                                end
                            end
                            
                            % iv. calculate mean directionality
                            meanD = mean(dirtrace(bin1:bin2));
                            if meanD < MINDIR
                                disp('HALT!')
                                if passnum == 1
                                    continue
                                elseif passnum == 2
                                    break
                                end
                            end
                            
                            % v. calculate mean directional bias
                            meanDB = nanmean(biastrace(bin1:bin2));
                            if abs(meanDB) < MINDIRBIAS 
                                disp('HALT2!')
                                if passnum == 1
                                    continue
                                elseif passnum == 2
                                    break
                                end
                            end
                            
                            % vi.  direction of motion  (from the sign of the weighted correlation)
                                % +1: traveling AWAY FROM the center well
                                % -1: traveling TOWARD the center well
                            motiondir = sign(wcorr(besttraj));  
                            
                            % temporary test 
                            if besttraj == 3
                                %disp('besttraj == 3, skipping')
                                break
                            end
                            
                            % vii. Forward vs. Reverse
                            replayT = nan;  % +1: forward, -1: reverse
                            if BIAS_CENTER
                                if motiondir == sign(meanDB)
                                    replayT = 1;
                                else
                                    replayT = -1;
                                end
                            else
                                disp('still have yet to code center-junction based replay type')
                                disp('in this case, direction needs to be analyzed differently on each arm')
                                keyboard
                            end
                                
                            % identify replay start and end time in epoch
                                r1 = a1 - 1 + bin1;
                                r2 = a1 - 1 + bin2;
                                
                                rstart = P.binvec_c(r1);
                                rend = P.binvec_c(r2);
                                
                                % sanity check
                                if length(r1:r2) ~= dur
                                    disp('something is wrong, these lengths should match')
                                end
                                                  
                                
                            % Store [rstart, rend, numbins, armcoverage, weighted corr, besttraj, F/R]
                            replay = [rstart, rend, dur, cvg, bestcorr, besttraj, replayT];
                            replayoutput{nn} = [replayoutput{nn} ; replay];
                            replayfound = 1;
                            
                            % (plot check)  plot a line indicating detected replay
                            if plotcheck
                                subplot(9,1,[1 2 3]);
                                if btype == 1
                                    replayclr = 'r';
                                elseif btype == 2
                                    replayclr = 'b';
                                end
                                plot([bin1 bin2],[10 10],'-','Color',replayclr,'linewidth',4); hold on
                                plot([bin1:bin2],10 * ones(1,dur),'x','Color',replayclr,'markersize',15,'linewidth',4); hold on
                            end
                            
                            if plotcheck
                                %pause
                                %close all
                            end
                            
                            %disp('replay identified, breaking')
                            break % if already identified a replay event, don't do a sub-period pass (i.e. passnum == 2)
                            
                           end
                       end
                    end
                    
                    if ploteeg && replayfound
                        subplot(9,1,[7 8 9])
                        
                        numbands = length(tracespec);
                        ztrace = cell(1,numbands);
                        timevecs = cell(1,numbands);
                        for g = 1:numbands
                            out = loadtracestruct(animalinfo{2}, animalinfo{3},tracespec{g}{1}, d, ep);
                            TF = tracespec{g}{2};
                            if 1
                                [ztrace{g},~,~] = zscoretrace(out{d}{ep}{TF},tracekern_ms);
                            elseif 1
                                ztrace{g} = zscorer(out{d}{ep}{TF}.powertrace);
                            end
                            timevecs{g} = out{d}{ep}{TF}.eegtimesvec_ref;
                            clear out
                        end
                        
                        t1 = P.binvec_c(a1) - binsize_t/2;
                        t2 = P.binvec_c(b1) + binsize_t/2;
                        for bb = 1:numbands
                            ind1 = lookup(t1,timevecs{bb});
                            ind2 = lookup(t2,timevecs{bb});
                                firsttime = timevecs{bb}(ind1);
                            plot(1000 * (timevecs{bb}(ind1:ind2)-firsttime),...
                                  ztrace{bb}(ind1:ind2),'Color',tracespec{bb}{3},'linewidth',3); hold on
                             axis tight
                             ylim([-2 10])
                        end
                        ylabel('Ripple/gamma Power (z)','FontSize',12);
                        set(gca,'fontsize',12)                            
                        
                        
                        subplot(9,1,[1 2 3])
                        if replayT == 1
                            replaytypestring = 'FORWARD';
                        elseif replayT == -1
                            replaytypestring = 'REVERSE';
                        end
                        string1 = sprintf('%s day %d ep %d',animalname(1:3),d,ep);
                        string2 = sprintf('Event time: %0.1f sec (%s,%s)',...
                                candstarts(c)-epstart,replaytypestring,besttrajstring);
                        string3 = sprintf('(CR: %0.1f CL: %0.1f RL:%0.1f)',wcorr);
                        title({string1,string2,string3},'fontsize',14,'fontweight','bold')                        
                    
                        
                        clear out
                        
                        pause
                        close all
                    end
                    
                    close all
            end
            
            
            
            if nn == 1
                estring = 'Local';
            elseif nn == 2
                estring = 'Remote';
            end
            
            if ~isempty(replayoutput{nn})
                numforward = sum(replayoutput{nn}(:,end) == 1) ;
                numreverse = sum(replayoutput{nn}(:,end) == -1) ;
            else
                numforward = 0;
                numreverse = 0;
            end
            
            disp(sprintf('%s %d %d : %s For: %d Rev: %d',...
                         animalname(1:3),d,ep,estring,numforward,numreverse))
            
        end
        
        
        % update and save P file
        P.replay = replayoutput;
        P.cand_count = cand_count;
        P.replay_nondur_count = replay_nondur_count;
        P.rep_descript = 'start,end,numbins,armcov,wcorr,besttraj,localremote,F/R';
        cd(datadir_pseudo)
        save(filename,'P','-v7.3')
        
        clear P D spikes
        
    end
    
    
end

    

if plot_events

    plot_replay_events

end


    


              % extra code for arm segmented MAP function (Wu--Foster)
                        %                     % arm segment version, not using for now
                        %                     MAP = cell(1,3);
                        %                     for ar = 1:3
                        %                         post = P.posteriors{nn}(arminds{nn}{ar},a:b);  % posterior probability matrix for the event, in this arm
                        %                         map = max(post,[],1);  % maximum probability in each decoding time bin
                        %                         % smooth map functions
                        %                         if 1
                        %                             mapkern = gaussian(1, 1 * 8);
                        %                             mapsmooth = smoothvect(map,mapkern);   % smooth MAP w/ 15 ms Gaussian
                        %                             MAP{ar} = mapsmooth;
                        %                         else
                        %                             MAP{ar} = map;
                        %                         end
                        %                     end
                        %                     % (plot check)  MAP functions
                        %                     if 1
                        %                         figure
                        %                         for ar2 = 1:3
                        %                             subplot(3,1,ar2);
                        %                             xvec = (1:numdecbins)*binsize_t*1000;
                        %                             plot(xvec,MAP{ar2},'b-','linewidth',5); hold on
                        %                             % also plot minimum probability level (absolute level)
                        %                             minprob = 5 * (1/length(P.linposbins{nn})) ;
                        %                             plot([xvec(1) xvec(end)],[minprob minprob],'r-','linewidth',3);
                        %                         end
                        %                         keyboard
                        %                         %close all
                        %                     end
    




