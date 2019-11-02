
% Find continuous Theta-associated alternations in Pos decodes 
allan = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander',...
                'Eight','Ten','Conley','Miles'};
mostan = {'Government','Egypt','Chapati','Dave','Frank','Bond','Corriander',...
                'Eight','Ten','Conley','Miles'};
core4   = {'Bond','Frank','Dave','Government'};
anim1st = {'Government','Egypt','Chapati','Dave','Higgs','Frank'};         
anim2nd = {'Bond','Corriander','Eight','Ten','Conley','Miles'};

animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander',...
    'Iceland','Justice','Kapital','Laplace','Eight','Ten','Conley','Miles'};

calculate                   = 0;
plot_pvalues                = 1;
plot_behavior               = 1;
report_longest_altpers      = 1;
    P_thresh_long_altpers   = 0.001;
report_prevalence           = 0;

    DECODE_RUN = 1;   % 0: Uniform, 1: Random walk, 100: (something old?)
    EXCURTYPE_toplot = [ 1 ] ; % 1: Pro, 2: Ret, 3: In L, 4: In R   
    
if plot_pvalues || plot_behavior || report_longest_altpers
   animals_toplot = mostan; %mostan; %{'Government'}; %allan; %{'Bond','Frank','Dave','Government'}; %{'Bond','Frank','Dave','Government'};  %'Bond','Frank','Dave',
end

if calculate
    
    animals_tocalc =  allan; %{'Egypt','Corriander'}; % 'Bond','Frank','Dave','Government','Dave','Corriander','Miles','Eight'};
        manual_dayep = [];

    % Analysis targets %
    PRECHOICE_ONLY      = 1;        % 0: full excur, 1: CP only
    EXCURTYPE_tocalc    = [1];      % 1: Pro, 2: Ret, 3: In L, 4: In R
    PHASECHOICE         = 1;        % 0: min firing phase
                                    % 1: min LR density  ** Choose this
                                    % 2: max LR density
    MOVINGPERIODS       = 2;        % Periods to analyze: 1: vel4, 2: nonimmobile05 (half sec)
    LR_prop_thresh      = 0.1;      % Minimum LR density (vs. C) to calculate LR value        

    % Alt period detection params %
    ALT_THRESHOLD       = 0.1; 
    NUMSAMP             = 10000 + 1;    % # of resamples
        SAMPMETHOD = 1;                 % 1: permutation, 2: bootstrap
    MAXALTS             = 5;    
    P_THRESH            = 0.05;

    % Theta %
    nbins = 12;  %18;  % 36 bins is 10 deg bins, after Jezek 2011
        phasebins       = -pi:(2*pi/nbins):pi;
        phasebins_c     = phasebins(1:(end-1));
end 

% Saved file directories %
if DECODE_RUN == 1
    decodedir = '/opt/data50/kkay/__Decode/LR_decode_20_4_Randwalk_new';    
    thetabindir = '/opt/data50/kkay/__Decode/LR_decode_20_4_Randwalk_new';
    cd(decodedir) 
elseif DECODE_RUN == 100
    decodedir = '/opt/data50/kkay/__Decode/Dir_decode_20_4_Uniform';   LR_decode_20_4_Randwalk_new     
    thetabindir = '/opt/data50/kkay/__Decode/Dir_decode_20_4_Uniform/Thetabins_LR';
    cd(decodedir)
end

% Internal study figures %%%%%%%%%%%%
plot_Xcorr_LR = 0;
plot_LRprop_hist = 0;
plot_LRprop_phasehist = 0;
plot_altspeed_hist = 0;

if plot_behavior || plot_pvalues
   BEHAVE_P_THRESH = 0.05;
   LENGTH_TOPLOT = 3; % [2 3]
end
if plot_pvalues
    pbins = 0:.01:.35;
end
    
    
% Data loading %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if calculate
    
    cd(decodedir)
    
    animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander',...
        'Iceland','Justice','Kapital','Laplace','Eight','Ten','Conley','Miles'};
    
    for XT = EXCURTYPE_tocalc
        
        EXCURTYPE = XT;
    
        out = [];
    
    for aa = 1:length(animals_tocalc)
        
        animalname = animals_tocalc{aa};
        animalpref = animalname(1:3);
        animalinfo = animaldef(animalname);
        daydir = getdaydir(animalname);
        an = find(strcmp(animalname,animals_order));
        task = loaddatastruct(animalinfo{2},animalinfo{3},'task');
        
        if ~isempty(manual_dayep)
            dayeps = manual_dayep;
        else
            days = 1:length(task);
            dayeps = [];
            for dd = days
                if ~isempty(task{dd})
                    eps = wtrackeps(task,dd);
                    dayeps = [dayeps ; dd * ones(length(eps),1)  eps(:)];
                end
            end
        end
        
        % Initialize outputs
        out.date               = date;
        out.animalname         = animalinfo{3};
        out.dayeps             = dayeps;
        out.numperms           = NUMSAMP;
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
            
            % Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~exist('P','var') || isempty(P) || ~strcmp(P.animalname,animalname) || ...
                    ~all(P.dayep == [d ep]) || EXCURTYPE ~= P.EXTYPE
                
                % Decode file %%%
                filename = sprintf('%s_%d_%d_Pos_%d.mat',animalpref,d,ep,EXCURTYPE);
                cd(decodedir);
                filedir = dir(filename);
                if isempty(filedir)
                    disp(sprintf('not finding %s, skipping',filename))
                    continue
                end
                load(filedir.name,'P');
                
                % Positional %
                pos = loaddatastruct(animalinfo{2},animalinfo{3},'pos',P.dayep(1));
                linpos = loaddatastruct(animalinfo{2},animalinfo{3},'linpos',P.dayep(1));
                lindist = linpos{d}{ep}.statematrix.lindist;
                CP_buff = choicepointer(linpos{d}{ep});
                CP_bin = ceil(CP_buff);
                segind = linpos{d}{ep}.statematrix.segmentIndex;
                seg2 = (segind == 2) ; seg3 = (segind == 3) ;
                seg4 = (segind == 4) ; seg5 = (segind == 5) ;
                Carm = (lindist <= CP_buff)    ;       % C arm pos times
                Larm = ~Carm & (seg2 | seg3)   ;       % L arm pos times
                Rarm = ~Carm & (seg4 | seg5)   ;       % R arm pos times
                postimevec = pos{d}{ep}.data(:,1);
                veltrace = pos{d}{ep}.data(:,5);
                epstart = postimevec(1);
                epend = postimevec(end);
                headang = loaddatastruct(animalinfo{2}, animalinfo{3}, 'headang', d);
                angveltrace = abs(headang{d}{ep}.angvel)  * 180 / pi;    % convert to degrees
                %[posvec_all] = posvecmaker(lindist,linpos{d}{ep});
                %[posvec_all,outersep] = posvecmaker_compress(posvec_all);
                yvec    = P.yvec_pos;
                outersep = P.CLR_lengths(1) + P.CLR_lengths(2);
                
                %                     % Head direction data %
                %                     headdir = linpos{d}{ep}.statematrix.segmentHeadDirection(:,1);   % head direction relative to center well -- values > 0 are outbound
                %                     postimevec_nonnan = postimevec(~isnan(headdir));
                %                     headdir_nonnan = headdir(~isnan(headdir));
                %                     headdir2 = interp1(postimevec_nonnan,headdir_nonnan,postimevec,'linear');
                %                     outvec = (headdir2 >= 0);
                %                     invec  = (headdir2 < 0);
                %                     disp(num2str(sum(invec & Larm)/30))
                %                     keyboard
                
                % By Excursion type, identify Choice Times + C,L,R pos bins %%%%
                Choicevec = [];
                if EXCURTYPE == 1       % Pro
                    exflags = [1 3 -11];
                    Choicevec = Carm;   % Pro Choice Times   : Center arm
                    Cin = yvec <= CP_bin;                          % Center arm
                    Lin = (yvec > CP_bin) & (yvec <= outersep);    % Left arm
                    Rin = ~Cin & ~Lin;                         % Right arm
                elseif EXCURTYPE == 2   % Ret
                    exflags = [2 4 -11];
                    Choicevec = Carm;   % Ret "Choice" Times : Center arm
                    Cin = yvec <= CP_bin;                          % Center arm
                    Lin = (yvec > CP_bin) & (yvec <= outersep);    % Left arm
                    Rin = ~Cin & ~Lin;                         % Right arm
                elseif EXCURTYPE == 3   % In L
                    exflags = [-33 4 32];
                    Choicevec = Larm;   % L inbound "Choice" Times : Left arm
                    Cin = (yvec > CP_bin) & (yvec <= outersep) ;    % Left arm
                    Lin = yvec > outersep;                          % Right arm (here the "left side")
                    Rin = ~Cin & ~Lin;                              % Right arm
                elseif EXCURTYPE == 4   % In R
                    exflags = [-22 2 23];
                    Choicevec = Rarm;   % R inbound "Choice" Times: Right arm
                    Cin = yvec > outersep;                          % Center arm
                    Lin = yvec <= CP_bin;    % Left arm
                    Rin = ~Cin & ~Lin;                         % Right arm
                end
                
                % Identify Excursions (Left, Right, Trackback) %%%
                inds = ismember(P.excurlist(:,3),exflags);
                excurlist   = P.excurlist(inds,1:3);
                numexcur = size(excurlist,1);
                
                % Moving periods %
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
                
                % Tets + clust spikes
                spikes = loaddatastruct(animalinfo{2},animalinfo{3},'spikes',d);
                cellinfo = loaddatastruct(animalinfo{2},animalinfo{3},'cellinfo');
                [adtc_list] = clusteredunits(spikes,cellinfo,an,d,ep,[],1);
                maxcell = size(adtc_list,1);
                selected_tets = unique(P.selected_tets);
                numtets = length(selected_tets);
                
                % Theta %
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
                tph = thetabins{d}.tph{ep};
                tvec_tph = thetabins{d}.timevec_tph{ep};
                %thetaphase = double(theta{d}{ep}{thetatet}.data(:,2))/10000;  % hilbert phase
                %thetatimevec = geteegtimes(theta{d}{ep}{thetatet});
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % I. Construct Tvec and LRdensity, i.e. decode concatenated across excursions
            % (Also, histogram LR density to get theta phase at which LR density is minimal (minphase)) %%%%%%%%%%%%%%%%%%%%%%
            
            Tvec = [];
            CLRdec = [];        % [row: times // col: L: 1, R: 2 ]
            Phasehist_LR = zeros(size(phasebins_c));
            Analysisperiods = []  ; % all analysis pers in epoch
            Analysis_totaldur = 0;
            
            % Get total analysis vector + periods for this epoch
            Analysisvec = false(size(postimevec));
            Tvec_collect = [];
            for x = 1:numexcur
                % Define analysis periods %   (Pre-choice point, Excursion times, Movement)
                excurvec    = list2vec(excurlist(x,[1 2]),postimevec);
                
                if PRECHOICE_ONLY
                    intersectvec    = Choicevec & excurvec & movingvec;   % *** Prior to Choice, This Excursion, Movement ****
                else
                    intersectvec    =  excurvec & movingvec;   % *** This Excursion, Movement ****
                end
                
                Analysisvec(intersectvec) = 1;
                % Collect (overlapping) decoding times
                Tvec_collect    =  [Tvec_collect ; P.Tvecs{x}(:)];
            end
            Tvec = unique(Tvec_collect);
            CLRdec = nan(length(Tvec),3);  % [ times x CLR ]
            Analysisperiods        = vec2list( Analysisvec,postimevec );
            Analysis_totaldur   = Analysis_totaldur + round( sum( Analysisperiods(:,2) - Analysisperiods(:,1) ));
            
            %   Iterate excursions to collect
            for xx = 1:numexcur
                
                % Obtain Posterior (C,L,R components)
                post_exc    = P.Posts{xx};   % Posterior during excursion
                tvec_exc    = P.Tvecs{xx};   % Time vector of bins during excursion
                Cdensity    = sum(post_exc(Cin,:),1);
                Ldensity    = sum(post_exc(Lin,:),1);
                Rdensity    = sum(post_exc(Rin,:),1);
                LRdensity       = Ldensity + Rdensity;
                LRthreshinds    = LRdensity >= LR_prop_thresh;
                
                tphs = tph(lookup(tvec_exc,tvec_tph - epstart));   % Theta phase of each cut density bin
                
                % Finally, install onto epoch-master Decode variable
                inds_ep = lookup(tvec_exc,Tvec);
                CLRdec(inds_ep,:)   =  [ Cdensity(:)  Ldensity(:)  Rdensity(:) ];
                
                % Also, add LR density to accumulator phase histogram
                [~,I] = histc(tphs,phasebins);
                for bbb = 1:nbins
                    %LRdensity_addition = sum(LRdensity( (I' == bbb) & LRthreshinds_cut(:) )); % Correct phase bin & meets threshold
                    LRdensity_addition = sum(LRdensity( (I' == bbb) )); % Correct phase bin & meets threshold
                    Phasehist_LR(bbb) = Phasehist_LR(bbb) + LRdensity_addition; % Add to total
                end
                
                % (plot check) Hist of LR densities during this excursion's analysis periods
                if plot_LRprop_hist
                    H = figure('units','normalized','outerposition',[.2 .06 .2 .2]);
                    hist(LRdensity,20); axis tight
                    title('LR density','fontweight','bold','fontsize',12)
                    keyboard
                end
                
            end
            
            disp(sprintf('%d s of epoch analyzed (EXCURTYPE: %d)',Analysis_totaldur,EXCURTYPE))
            
            if Analysis_totaldur == 0
                disp('hmm, skipping')
                continue
            end
            
            %%% II. Divide epoch into Theta bins using LR minphase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % (a) Select phase + divide theta cycles %%%%%%%%%%%%%%%%%%%%%%%
            
            if PHASECHOICE == 0
                tbins = thetabins{d}.thetabins{ep};  % min firing phase (this is not preferable, since LFP variable)
                tbblocks = thetabins{d}.tbblocks{ep};
                minphase = thetabins{d}.minphase;
            elseif PHASECHOICE == 1
                [~,I2] = min(Phasehist_LR);  % min
                minphase = phasebins_c(I2);
                [tbins, ~, ~, ~, ~, ~, tbblocks]  ...
                    = thetabinner(tph,tvec_tph,minphase,[epstart epend]);
            elseif PHASECHOICE == 2
                [~,I2] = max(Phasehist_LR);  % max
                minphase = phasebins_c(I2);
                [tbins, ~, ~, ~, ~, ~, tbblocks]  ...
                    = thetabinner(tph,tvec_tph,minphase,[epstart epend]);
            end
            
            % (Optional plot check) Theta phase histogram of LR density (ep-wide) %%%%%
            if 0 %plot_LRprop_phasehist
                H = figure('units','normalized','outerposition',[.1 .3 .2 .2]);
                B = bar(phasebins_c,Phasehist_LR,'histc'); hold on
                set(B,'facecolor','k')
                set(gca,'fontsize',14)
                axis tight
                plot([minphase minphase],[0 max(Phasehist_LR)],'r-','linewidth',3); hold on
                title(sprintf('min phase : %0.3f',minphase))
                keyboard
            end
            
            % (b) Store this in a file (for plotting later)
            clear tbins_LR;
            tbins_LR.date = date;
            tbins_LR.PHASECHOICE = PHASECHOICE;
            tbins_LR.minphase = minphase;
            tbins_LR.tbins = tbins;
            tbins_LR.tbblocks = tbblocks;
            cd(thetabindir);
            filename_tb = sprintf('Thetabins_LR_%d_%s_%d_%d.mat',EXCURTYPE,animalname(1:3),d,ep);
            save(filename_tb,'tbins_LR','-v7.3');
            cd(decodedir)
            
            % (c) Filter for Theta bins that occur within Analysis periods %%%%%%%%%%%%%%%
                filterinds = logical(  isExcluded( mean(tbins,2) , Analysisperiods )  ); % Get theta bins that occur analysis periods
            Tbins = tbins(filterinds,:);
            Tbins_mean = mean(tbins,2);
            Numcyc = size(Tbins,1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%% III. Calculate LR Decode + Alt speed + Indicator vectors / Theta bin  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Tbins;                              % [ starttime endtime] of theta cycle
            Tbins_mean      = mean(Tbins,2);    % mid-time of theta cycle
            Contig_tb       = [];               % 111122221111: Contig indicator vector for theta bins
            LR_tb           = nan(Numcyc,1);    % R/(R+L) proportion / theta cycle
            Alt_tb          = nan(Numcyc,1);    % LR alternation speed / theta cycle
            Excur_tb    = [];               % Actual excursion-type (L, R, Trackback) the rat was on during the Theta bin
            
            % Calculate LR index (0: Left, 1: Right)
            % Calculate Alt speed
            % Construct Indicator vectors
            
            % Iterate each theta cycle
            for tb = 1:Numcyc
                
                bin_a = lookup(Tbins(tb,1) - epstart,Tvec);  % Start index of bin
                bin_b = lookup(Tbins(tb,2) - epstart,Tvec);  % End index of bin
                midtime = mean(Tbins(tb,:));
                %disp(num2str(midtime-epstart))
                
                % (Provisional - fix later) Identify the excursion-type of this bin
                excur = nan;
                for xxxx = size(excurlist,1):-1:1
                    if sum(isExcluded(midtime,excurlist(xxxx,[1 2]))) > 0
                        excur = [excur ; excurlist(xxxx,3)];
                        break
                    end
                end
                if isnan(excur)
                    excur = nan;
                    disp('nan excur -- should be rare')
                %elseif length(unique(excur)) > 1
                %    keyboard
                %else
                %    excur = unique(excur);
                end

                Excur_tb = [Excur_tb ; excur];
                
                %                  % Check discrepancy due to decoding bin size (e.g. 20 ms)
                %                     decodedur = Tvec(bin_b) - Tvec(bin_a);
                %                     actualdur = Tbins(tb,2) - Tbins(tb,1);
                %                     diffdur = abs(decodedur - actualdur) ;
                %                     disp(num2str(diffdur));
                %                     if diffdur > 0.008
                %                         keyboard
                %                     end
                
                % Post density of C, L, R arms
                Cval = sum(CLRdec(bin_a:bin_b,1));     % Center arm density of this theta bin
                Lval = sum(CLRdec(bin_a:bin_b,2));     % Left " "
                Rval = sum(CLRdec(bin_a:bin_b,3));     % Right " "
                
                if any(isnan([Cval Lval Rval]))
                    %keyboard
                    disp('*********nan CLR val (should only happen occasionally) ****************')
                    disp(num2str([Cval Lval Rval]))
                end
                
                % Calc LR proportion
                LR_prop = (Lval + Rval)/(Cval + Rval + Lval);
                
                % Calc LR (if exceed min LR density)
                if LR_prop > LR_prop_thresh
                    
                    LRFLAG = 2;
                    
                    if LRFLAG == 1
                        % Option 1: Proportion
                        LRval =  Rval / (Lval + Rval);
                    elseif LRFLAG == 2
                        % Option 2: Sign
                        LRval =  Rval / (Lval + Rval);
                        LRval =  sign(LRval - 0.5);  % +1: Right, -1: Left
                    end
                    
                    LR_tb(tb) =  LRval;
                else
                    LR_tb(tb) = nan;
                end
                
                if tb > 1
                    % Check to make sure cycle did not come after discontinuous theta bin gap
                    if Tbins(tb,1) == Tbins(tb-1,2)  % Contiguous
                        
                        % Calc Alt speed %
                        if ~isnan(LR_tb(tb-1))
                            Alt_tb(tb) = abs( LR_tb(tb) - LR_tb(tb-1) );
                        end
                        
                        % Mark Contiguity %
                        if Contig_tb(end) == 1
                            Contig_tb = [Contig_tb ; 1];  % keep
                        else
                            Contig_tb = [Contig_tb ; 2];   % keep
                        end
                        
                    else           % Discontiguous
                        
                        % Mark Discontiguity %
                        if Contig_tb(end) == 1
                            Contig_tb = [Contig_tb ; 2];  % switch
                        else
                            Contig_tb = [Contig_tb ; 1];   % switch
                        end
                    end
                elseif tb == 1
                    Contig_tb = [Contig_tb ; 1];  % Initialize w/ 1
                    continue
                end
                
            end
            
            if 0  %(optional) % check alternations
                if LRFLAG == 1
                    figure;
                    plot(Tbins_mean-epstart,Alt_tb,'k.'); hold on
                    plot(Tbins_mean-epstart,Alt_tb,'k-'); hold on
                    ylim([-0.1 2.1])
                elseif LRFLAG == 2
                    plot(Tbins_mean-epstart,Alt_tb,'r.'); hold on
                    plot(Tbins_mean-epstart,Alt_tb,'r-'); hold on
                    ylim([-0.1 2.1])
                end
                %break
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % IV. Detect Alternation Periods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Initialize output
            Altpers = [];
            Altpers_xc = [];      % excursion type of the alt period
            Altpers_dur = [];
            numaltperiods = nan;  % total # of alternation periods this ep
            
            % Initialize while-loop state variables
            cy = 2;        % pointer index
            altcycdur = 0; % alt per dur
            per_a = 1;     % alt per start
            per_b = nan;   % alt per end
            
            while cy ~= Numcyc
                
                if Contig_tb(cy) == Contig_tb(cy-1)         % Contiguous
                    
                    if Alt_tb(cy) > ALT_THRESHOLD    % this theta bin qualifies as alternation
                        altcycdur = altcycdur + 1;   % alternation duration (in cycles)
                        per_b = cy;                  % update pointer
                    else            % not reach alternation thresh
                        % Store alternation period
                        if altcycdur > 2
                            Altpers = [Altpers ; per_a  per_b ];
                        end
                        % Reset pointers
                        altcycdur = 0;
                        per_a = cy;
                        per_b = nan;
                    end
                    
                else                                        % Discontiguous
                    
                    %Reset pointers
                    altcycdur = 0;
                    per_a = cy;
                    per_b = nan;
                    
                end
                
                cy = cy + 1;
                
            end
            
            numaltperiods  = size(Altpers,1);
            
            if numaltperiods == 0
                disp('no alt periods detected')
                continue
            end
            
            % Identify excursion types for alt periods
            Altpers_xc     = nan(numaltperiods,1);
            for vv = 1:numaltperiods
                extype = unique( Excur_tb( Altpers(vv,1) : Altpers(vv,2) )  );
                if length(extype) ~= 1
                    disp('ALTPERIOD STRADDLES TWO EXCUR: ASSUMING TRACKBACK')
                    [~,III] = max(abs(extype));
                    extype = extype(III);
                    %keyboard
                end
                Altpers_xc(vv) = extype;
            end
            
            if isempty(Altpers)
                continue
            end
            
            % Tabulate cycle durations
            Altpers_dur = [ Altpers(:,2) - Altpers(:,1) ] + 1;            % Duration, in theta cycles, of each candidate period
            for mm = 1:MAXALTS
                out.numaltpers_ep{de}(1,mm) = sum(Altpers_dur >= mm);   % Periods with at least mm cycles
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % V. Permute single alt periods to get P-values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Initialize output
            %     [ time_a   time_b   tb_a   tb_b   numcyc   xcurtype   speed   accel   angspeed   pvalue ]
            out.altperiods{de} = nan(numaltperiods,10);
            
            % (i) Preliminary: obtain LRvals for each excursion type (within this epoch)
            cycinds_xc = {};    % indices of theta cycles of this excursion type
            numcyc_xc = {};     % # of theta cycles of this excursion type
            LRvals_xc = {};     % LR values for each of these theta cyccles in this excursion type
            for XF_flag = exflags   % iterate through each excursion flag value (only positive!)
                    xf_inds = excurlist(:,3) == XF_flag;
                    XF = abs(XF_flag);
                cycinds_xc{XF}      = find( isExcluded(Tbins_mean,excurlist(xf_inds,[1 2])) );    % Right
                numcyc_xc{XF}       = length(cycinds_xc{XF});
                LRvals_xc{XF}       = LR_tb(cycinds_xc{XF});
            end
            
            % Iterate each alt period
            for ll = 1:numaltperiods
                
                % Basic info %
                tb_a    = Altpers(ll,1);        % theta bin # of first cycle
                tb_b    = Altpers(ll,2);        % theta bin # of last cycle
                time_a  = Tbins(tb_a,1);        % absolute time of start of first cycle
                time_b  = Tbins(tb_b,2);        % absolute time of end of last cycle
                numcyc  = tb_b - tb_a + 1;      % number of theta cycles
                xf      = abs(Altpers_xc(ll));  % excursion flag (abs makes it indexable)
                
                if isnan(xf)
                    disp('xcursion type cant be found')
                    keyboard
                    %continue
                end
                
                disp(sprintf('%s: d%d ep%d xc%d altper #%d (%d)',animalinfo{3}(1:3),d,ep,xf,ll,numcyc))
                
                % Behavior (calc mean of speed, accel, angular speed) %
                pos_a           = lookup(time_a,postimevec);
                pos_b           = lookup(time_b,postimevec);
                speed           = mean(veltrace(pos_a:pos_b));          % cm/s
                accel           = mean(diff(veltrace(pos_a:pos_b)));    % cm/s^2
                angspeed        = mean(angveltrace(pos_a:pos_b));       % deg / s
                lowspeedflag    = any(veltrace(pos_a:pos_b) < 4);       % overlaps with low-speed time
                
                % Permutation procedure (i, ii) %%%%%
                
                % (i) Identify indices of the surrounding theta bin contiguous period (or "block period")
                %       in which this alternation period occurred
                %       strategy is to use a while loop and step backwards and
                %       forwards
                tbc_a = [];  % first bin # of surrounding contig theta period (with decoded LR)
                tbc_b = [];  % last bin # of surrounding contig theta period (with decoded LR)
                
                % Initialize %
                tbc_a = tb_a;  % first bin # of surrounding contig theta period (with decoded LR)
                tbc_b = tb_b;  % last bin # of surrounding contig theta period (with decoded LR)
                blockval = tbblocks(tb_a);  % dummy val is 1 or 2: designates a contiguous theta period
                
                % Find beginning by stepping backwards %
                while tbc_a > 0
                    if tbblocks(tbc_a) == blockval && ~isnan(LR_tb(tbc_a)) % Same contig theta block && LR is decoded
                        tbc_a = tbc_a - 1;
                    else
                        break  % terminate
                    end
                end
                
                % Find end by stepping forward %
                while tbc_b <= Numcyc
                    if tbblocks(tbc_b) == blockval && ~isnan(LR_tb(tbc_b)) % Same contig theta block && LR is decoded
                        tbc_b = tbc_b + 1;
                    else
                        break  % terminate
                    end
                end
                
                blockdur = tbc_b - tbc_a + 1;  % Duration (in cycles) of surrounding block period
                
                % (ii) Permute
                altdetect = nan(1,NUMSAMP);
                for rrr = 1:NUMSAMP
                    
                    % Resample %
                    %   (here we just go ahead and resample all the LR vals in the epoch)
                    if SAMPMETHOD == 1
                        LRvals_resamp = LRvals_xc{xf}(randperm(numcyc_xc{xf}));  % Permutation (order)
                    elseif SAMPMETHOD == 2
                        LRvals_resamp = LRvals_xc{xf}(ceil( numcyc_xc{xf} * rand(1,numcyc_xc{xf}) )) ;  % Bootstrap
                    end
                    
                    % Extract out a block of the same size as surrounding contig block period
                    %       (do from the the beginning as this is simplest)
                    if blockdur < numcyc_xc{xf}         % essentially always the case
                        LR_resamp = LRvals_resamp(1:blockdur);
                    elseif numcyc_xc{xf} == 0
                        disp('Zero numcyc_xc -- should NOT happen but skipping this alternation')
                        keyboard
                        continue
                    elseif blockdur > numcyc_xc{xf}  % sometimes (as in trackback) this will happen --
                        LR_resamp = LRvals_resamp;
                        if abs(xf) < 10
                            disp('non-trackback excursion type has a lack of permute cycles? unexpected')
                            keyboard
                        end
                    end
                    
                    
                    
                    % Detect alternations
                    altspeed_resamp = abs(diff(LR_resamp));
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
                %%%%%%%  end of permute code  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % install outputs
                out.altperiods{de}(ll,1) = d;  % [ day ep  time_a  time_b   tb_a   tb_b  numcyc  xcurtype  speed accel angspeed pvalue ]
                out.altperiods{de}(ll,2) = ep;
                out.altperiods{de}(ll,3) = time_a - epstart;
                out.altperiods{de}(ll,4) = time_b - epstart;
                out.altperiods{de}(ll,5) = tb_a;
                out.altperiods{de}(ll,6) = tb_b;
                out.altperiods{de}(ll,7) = numcyc;
                out.altperiods{de}(ll,8) = xf;
                out.altperiods{de}(ll,9) = speed;
                out.altperiods{de}(ll,10) = accel;
                out.altperiods{de}(ll,11) = angspeed;
                out.altperiods{de}(ll,12) = lowspeedflag;
                out.altperiods{de}(ll,13) = pval;
                
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %  Printed reportage of alt periods >= 3 cycles
            NUM_LONGALT_REPORT = 3;
            alttimes_long           = Tbins(Altpers(find(Altpers_dur >= NUM_LONGALT_REPORT),:))-epstart;
            alttimes_long_durs      = Altpers_dur(Altpers_dur >= NUM_LONGALT_REPORT);
            alttimes_long_sig       = out.altperiods{de}(Altpers_dur >= NUM_LONGALT_REPORT,end) < P_THRESH;
            alttimes_long_pval      = out.altperiods{de}(Altpers_dur >= NUM_LONGALT_REPORT,end);
            alttimes_long           = [alttimes_long    alttimes_long_durs  alttimes_long_pval  alttimes_long_sig];
            
            if ~isempty(alttimes_long)
                disp('long alt times:');
                alttimes_long
            else
                disp('no long alt times')
            end
            
            
        end
        
        cd(decodedir);
        filename = sprintf('Analysis_%d_%s',EXCURTYPE,animalpref);
        save(filename,'out','-v7.3');
        
    end
    
    end
    
end
    
    
    
    if plot_pvalues
        
        % Collect p-values
        pvalues_collect = {};
        for L = 2:3  % 2 is exactly 2 alternations, 3 is >= 3
            for M = 1:2   % 1: non-low speed,  2: includes low-speed
                pvalues_collect{L}{M} = [];
                for aa = 1:length(animals_toplot)
                    animalname = animals_toplot{aa};
                        animalinfo = animaldef(animalname);
                    cd(decodedir)
                    filename = dir(sprintf('Analysis_%d_%s.mat',EXCURTYPE_toplot,animalname(1:3)));
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
        
          % p-value histogram (linear version)
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
            for L = 3  % 2 is exactly 2 alternations, 3 is >= 3
                behaviorvals{L}{M} = [];
                dayeptimes{L}{M} = [];
                for aa = 1:length(animals_toplot)
                    animalname = animals_toplot{aa};
                    cd(decodedir)
                    filename = dir(sprintf('Analysis_%d_%s.mat',EXCURTYPE_toplot,animalname(1:3)));
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
            set(gca,'xtick',0:5:40)
            %zlabel('accel','fontweight','bold','fontsize',12)
            grid on
            xlim([0 40])
            ylim([0 10])
            set(gca,'fontsize',14)
        end
        
        animstring = mat2str(cell2mat(animals_toplot));
        subplot(2,1,1);
        title(sprintf('%s: Behav vals, alt pers',animstring),'fontsize',14,'fontweight','bold')
        
    end
   
    
    
    
    
    
    if report_longest_altpers

        
        % Collect behavior values
        altpers_all  = [];  % 13 columns (same as out.altperiods above)
        %                     1   2     3       4        5      6     7       8     9     10    11       12             13       14  
        %                   [ a   d     ep    time_a  time_b   tb_a   tb_b  durcyc  xc  speed  accel  angspeed    lowspeedflag pvalue ]

        for M = 1:2
            for aa = 1:length(animals_toplot)
                animalname = animals_toplot{aa};
                an = find(strcmp(animalname,animals_order));
                cd(decodedir)
                filename = dir(sprintf('Analysis_%d_%s.mat',EXCURTYPE_toplot,animalname(1:3)));
                load(filename.name,'out')
                for de = 1:length(out.altperiods)
                    if ~isempty(out.altperiods{de}) && ~isempty(out.altperiods{de})
                        pinds = out.altperiods{de}(:,end) < P_thresh_long_altpers;
                        cycinds = out.altperiods{de}(:,7) >= 3;
                        if M == 1
                            speedinds = out.altperiods{de}(:,12) == 0;
                        elseif M == 2
                            speedinds = out.altperiods{de}(:,12) == 1;
                        end
                        % behavior values
                        allinfo = out.altperiods{de}(pinds & cycinds & speedinds, :);
                            numalts = size(allinfo,1);
                        altpers_all  = [altpers_all  ;  an * ones(numalts,1) allinfo];
                    end
                end
            end
        end
          
        altpers_all = sortrows(altpers_all,-8);
        
        % Report [an day numcyc]
        disp('[an day numcyc]')
        disp(num2str(altpers_all(:,[1 2 8])));
        
        % Report [an day ep numcyc time_a time_b]
        disp('[an day ep numcyc time_a time_b]')
        disp(num2str(altpers_all(:,[1 2 3 8 4 5])));
        
    end
    
    
    
    
    
    
    
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
    
   
  
