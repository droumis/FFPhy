function kk_clusterdecode(directoryname,fileprefix,dayeps,animalname,varargin)
%

         
animalinfo = animaldef(animalname);
   daydir = getdaydir(animalname);
blank = '';
tetfilter = '';
savedir = '/opt/data13/kkay/___Superfourcode_data';
remoteW = 0;
extratime = 500;  % ms of extra time before and after ripple to decode
plot_powertraces = 0;
plot_infunction = 0;
xdel = 1;
xkern = 2;  % # of xdel bins to smooth occupancy / firing rate maps
winsize = 0.5;
min_activecells = 3;
manual_period = [];

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'savedir'
            savedir = varargin{option+1};    
        case 'plot_infunction'
            plot_infunction = varargin{option+1};
        case 'remoteW'
            remoteW = varargin{option+1};
        case 'sigma_transmat'
            sigma_transmat = varargin{option+1};            
        case 'cellchoice'
            cellchoice = varargin{option+1};
        case 'extratime'
            extratime = varargin{option+1};        
        case 'winsize'
            winsize = varargin{option+1};        
        case 'manual_period'
            manual_period = varargin{option+1};    
        case 'plot_powertraces'
            plot_powertraces = varargin{option+1};      
        case 'exclude_inters'
            exclude_inters = varargin{option+1};
    end
end

task = loaddatastruct(directoryname,fileprefix,'task');

epochfilter = epochmaker('runW_rip');

% Identify all day eps
dayeps_all = evaluatefilter(task,epochfilter);

if isempty(dayeps)
    dayeps = dayeps_all;
end

% Iterate through epochs to decode

for de = 1:size(dayeps,1)
    
    day = dayeps(de,1);
    ep = dayeps(de,2);  % "local" ep
    
    disp(sprintf('%s day %d ep %d',animalname(1:3),day,ep))
    
    clear E D
    
    % Load data
    linpos = loaddatastruct(directoryname,fileprefix,'linpos',day);
    pos = loaddatastruct(directoryname,fileprefix,'pos',day);
    cellinfo = loaddatastruct(directoryname,fileprefix,'cellinfo',day);
    trajencode = loaddatastruct(directoryname,fileprefix,'trajencode',day);
    spikes =  loaddatastruct(directoryname,fileprefix,'spikes',day);
    if plot_powertraces
        out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'rippletrace', day, ep);
        riptrace = zscorer(out{day}{ep}{1}.powertrace);
        riptrace_timevec = out{day}{ep}{1}.eegtimesvec_ref;
        out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'fastgammatrace', day, ep);
        wgtrace = zscorer(out{day}{ep}{2}.powertrace);
        wgtrace_timevec = out{day}{ep}{2}.eegtimesvec_ref;
    end
    
   % (if specified) Identify the remote W-track epoch # 
   if remoteW
      task = loaddatastruct(directoryname,fileprefix,'task',day);
      eps = dayeps_all(dayeps_all(:,1)==day,2)';
      currenv = task{day}{ep}.environment;
      remoteeps = [];
      for xx = eps
         if ~strcmp(currenv,task{day}{xx}.environment)
             remoteeps = [remoteeps xx];
         end
      end
      remoteep = max(remoteeps);  %take the latest epoch (most efficient running behavior)
      epochs_enc = [ep remoteep];
      clear task
   else
       epochs_enc = ep;
       remoteep = [];
   end
  
   num_encodeeps = length(epochs_enc);   % [local epoch #, remote epoch #]

   % From each encoding epoch, collect basic epoch data
   postimevec = {};
   numpossamps = [nan nan];
   starttime = [nan nan];
   endtime = [nan nan]; 
   for zz = 1:num_encodeeps
        ep2 = epochs_enc(zz);
        postimevec{zz} = linpos{day}{ep2}.statematrix.time;
        numpossamps(zz) = length(linpos{day}{ep2}.statematrix.time);
        starttime(zz) = postimevec{zz}(1);
        endtime(zz) = postimevec{zz}(end);
   end   
   
   
   % Identify cells to encode-decode with (for local and remote W-track epochs)
   detc{1} = [];
   detc{2} = [];
   if cellchoice == 0
       % 0: cellinfo - based  (not using)
       cellfilter = 'isequal($type, ''principal'') && ( isequal($area, ''CA1'') || isequal($area, ''CA3'') || isequal($area, ''CA2'') )';  
       all_detc = evaluatefilter(cellinfo,cellfilter);
            cellinds{1} = (all_detc(:,1) == day) & ((all_detc(:,2) == ep) | (all_detc(:,2) == remoteep));
            cellinds{2} = (all_detc(:,1) == day) & ((all_detc(:,2) == ep) | (all_detc(:,2) == remoteep));
        detc = {};
        detc{1} = all_detc(cellinds{1},:);
        detc{2} = all_detc(cellinds{2},:);
            numunits(1) = size( detc{1},1) ;   
            numunits(2) = size( detc{2},1) ;  
   elseif cellchoice == 1
       % 1: Only CA1, CA2 P, CA3 units, and in addition checks whether place-active
       %        for each encoding epoch (LOCAL and, if available, REMOTE W-track epoch), 
            animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
            an = find(strcmp(animalname,animals_order));

                % i. load all principals (Right and Left hemispheres) for this animal
            cellfilter = 'isequal($type, ''principal'') && ( isequal($area, ''CA1'') || isequal($area, ''CA3'') || isequal($area, ''CA2'') )';  
            pr_detc = evaluatefilter(cellinfo,cellfilter); 
                    % identify those active in the present epoch
                pr_detc = pr_detc(  (pr_detc(:,1) == day) & (pr_detc(:,2) == ep),:);
            
                % ii. exclude CA2 N units
           load('/opt/data13/kkay/Superlin_data/Classtable_09-Apr-2015.mat');
                ca2n_adtc = sortrows([ classtable6_adtc{2}{2} ; classtable6_adtc{2}{3} ],[1 2 3 4]);
                ca2n_dtc = ca2n_adtc(ca2n_adtc(:,1) == an,2:4);
                if ~isempty(ca2n_dtc)
                   for qq = size(pr_detc,1):-1:1
                        if rowfind(pr_detc(qq,[1 3 4]),ca2n_dtc)
                            disp('ca2n unit ignored')
                            pr_detc(qq,:) = [];
                        end
                   end
                end

           dtc_clust = unique(pr_detc(:,[1 3 4]),'rows');     
                numclust = size(dtc_clust,1);
           adtc_clust = [an * ones(numclust,1)   dtc_clust  ];
           
                % iii. Use activematrix (record of whether unit was place-active in W-track epochs), generated by save_activeW_units
                   %    To filter only for place-active units
           load('/opt/data13/kkay/__WG/ACTIVEMATRIX.mat','ACT')
           placethresh = 2;   % in Hz
                col = find(ACT.ratethresh_field == placethresh) + 1;
           
           % now, for the LOCAL and REMOTE W-track epochs, respectively,
           %    identify the subset of these clustered units that were place-active
           numunits = [nan nan];
           for oo = 1:num_encodeeps
                ep2 = epochs_enc(oo);
                for cc2 = 1:size(adtc_clust,1)
                    rind = rowfind(adtc_clust(cc2,:),ACT.adtc);
                    if rind > 0
                        entry = ACT.actmatrix{rind};
                        rind2 = find(entry(:,1) == ep2);
                        if ~isempty(rind2)
                            act_flag = entry(rind2,col);   % place-active flag -- 1 if active, 0 if not active
                        else
                            act_flag = 0;  % also, if not clustered, then 0 
                        end
                    else
                        disp('cannot find unit in ACT')
                        keyboard                        
                    end
                    % if place-active, add to corresponding unit list
                    if act_flag == 1
                        detc{oo} = [detc{oo} ; day ep adtc_clust(cc2,[3 4]) ];  % if place-active, then add to list
                    end
                end
                numunits(oo) = size( detc{oo},1) ;
           end
           disp(sprintf('%d units clustered, %d W1 place-active, %d W2 place-active',numclust,numunits(1),numunits(2)))
           if numunits(1) < 12 && numunits(2) < 12
               disp('too few units')
               continue
           end
   end
       
   
   % Initialize outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
      % positional variables
   xbins = cell(1,2);
   armdists = cell(1,2);        % two elements since up to two different W-tracks
   validxbins = cell(1,2);   
   occ = cell(1,2);
   occ_sm = cell(1,2);
   centerarmmax = [nan nan];
   rightarmmax = [nan nan];
   maxarmdist = [nan nan];
     lastcenterbin = [nan nan];
     lastrightbin = [nan nan];
     % spike data
   spk_e = cell(1,2);
    
   % internal use, position variables (not to save in output file later)
   armdists2 = cell(1,2);
   armdists_cat = cell(1,2);
      a_cut = cell(1,2);            % these are the boundaries of each positional "island"  (for separate smoothing)
      b_cut = cell(1,2);
      c_cut = cell(1,2);
      d_cut = cell(1,2);
      e_cut = cell(1,2);
      f_cut = cell(1,2);
  encodeperiods = cell(1,2);    

   %%% Iterate through each encoding epoch to collect basic data %%%%%%%%%%%%%%%%%%%%%%%%      
   for nn = 1:num_encodeeps
       
       epenc = epochs_enc(nn);
       
       disp(['Encoding : Day ',num2str(day), ', Epoch ',num2str(epenc)])
       
       % Identify 2 SD SWRs (here, to use to exclude from encoding model)
       ripout = kk_getconstimes(animalinfo{2},animalinfo{3},[day epenc ],'ripplescons',1,...
           'consensus_numtets',3,'minthresh',2,...
           'exclusion_dur',0,'minvelocity',0,'maxvelocity',4);
       consvec_rip2 = ripout{day}{epenc}.cons;
       consvectimes_rip2 = ripout{day}{epenc}.time;
       periodtimes_rip2 = vec2list(consvec_rip2,consvectimes_rip2);
       
       % Basic epoch data
       postimevec_epenc = linpos{day}{epenc}.statematrix.time;
            numpossamps_encep = length(postimevec_epenc);
            starttime_epenc = postimevec_epenc(1);
            endtime_epenc = postimevec_epenc(end);
       
       % If later using modelnum 4, then get high velocity periods
       if modelnum == 4
           timefilterscript;
           tfout = evaluatetimefilter(animalinfo{2},animalinfo{3},{vel10},[day epenc]);
           hispeedperiods = tfout{day}{epenc};
       end               
            
       % Identify W-track segment indices
       lindist = linpos{day}{epenc}.statematrix.lindist;
       seg1 = linpos{day}{epenc}.statematrix.segmentIndex==1;
       seg2 = linpos{day}{epenc}.statematrix.segmentIndex==2;
       seg3 = linpos{day}{epenc}.statematrix.segmentIndex==3;
       seg4 = linpos{day}{epenc}.statematrix.segmentIndex==4;
       seg5 = linpos{day}{epenc}.statematrix.segmentIndex==5;
       
       % Initialize position variables : arm dists %%%%%%%%%%%%%%%%%%%
       armdists{nn}  = nan(1,numpossamps_encep) ;   % (horizontal vector)   +200 / +400 / +600 linear distance from center well     
            centerarmmax(nn) =    max( lindist(seg1) ) ;                                % maximum linear distance of center arm
            rightarmmax(nn)  =  max( lindist(seg4 | seg5)     -  centerarmmax(nn) ) ;   % maximum distance from the center junct to right well
            maxarmdist(nn)   =  max( lindist(seg2 | seg3)  )  +  rightarmmax(nn)  ;     % maximum armdist (end of the 1D plot)
                lastcenterbin(nn) = ceil(centerarmmax(nn));
                lastrightbin(nn) = ceil(rightarmmax(nn))+ceil(centerarmmax(nn));
            armdists{nn}(seg1)        = lindist(seg1)                               +   200;  % >200: center arm
            armdists{nn}(seg4 | seg5) = lindist(seg4 | seg5)  - centerarmmax(nn)  +   400;  % >400: right arm
            armdists{nn}(seg2 | seg3) = lindist(seg2 | seg3)  - centerarmmax(nn)  +   600;  % >600: left arm
                armdists{nn} = armdists{nn}(:);

                
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       %%%% Construct positional space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %xbins{nn} = min(armdists{nn}):xdel:(max(armdists{nn}) + xdel);        %
       xbins{nn} = 0:xdel:1000;        %
       
       %indices for each positional "island" a to b, c to d, e to f
       a_cut{nn} = lookup(min(armdists{nn}(seg1)),xbins{nn},-1);
       b_cut{nn} = lookup(max(armdists{nn}(seg1)),xbins{nn},-1);
       c_cut{nn} = lookup(min(armdists{nn}(seg4 | seg5)),xbins{nn},-1);
       d_cut{nn} = lookup(max(armdists{nn}(seg4 | seg5)),xbins{nn},-1);
       e_cut{nn} = lookup(min(armdists{nn}(seg2 | seg3)),xbins{nn},-1);
       f_cut{nn} = lookup(max(armdists{nn}(seg2 | seg3)),xbins{nn},-1);
       
       validxbins{nn}([a_cut{nn}:b_cut{nn} ...
                       c_cut{nn}:d_cut{nn} ...
                       e_cut{nn}:f_cut{nn} ]) = 1 ;   % indicates which xbins indices were actually occupied by animal
            disp(sprintf('%d bins valid',sum(validxbins{nn})))
            validxbins{nn} = logical(validxbins{nn});   
            
                
       % armdists_cat (used for plotting and output later -- same as armdists but with no 200-400-600 buffer) 
       armdists_cat{nn}              = nan(1,numpossamps_encep);
       armdists_cat{nn}(seg1)        = lindist(seg1)                               ;  % 
       armdists_cat{nn}(seg4 | seg5) = lindist(seg4 | seg5)    ;  %
       armdists_cat{nn}(seg2 | seg3) = lindist(seg2 | seg3)   +   rightarmmax(nn) ;  % 
        
       % Calculate occupancy %%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        % if eliminating high speed periods w/ modelnum 4
        if modelnum == 4
            validposinds = logical(isExcluded(postimevec{nn},hispeedperiods));
           armdists_filt = armdists{nn}( validposinds );
           occ{nn} = histc(armdists_filt,xbins{nn}) * 1/29.97;
        else
            occ{nn} = histc(armdists{nn},xbins{nn}) * 1/29.97;
        end       
       
            % sanity check for stitching procedure
            if 0
               figure;
                    % island plot
               subplot(1,2,1)
               plot(1:1001,occ{nn});
                    % stitch plot
               subplot(1,2,2)
               cutocc = occ{nn}(validxbins{nn});
               plot(1:length(cutocc),cutocc); hold on
                        % stitch dividers
                        plot([lastcenterbin(nn) lastcenterbin(nn)],[0 max(cutocc)],'k--')
                        plot([lastrightbin(nn)  lastrightbin(nn)],[0 max(cutocc)],'k--')
                        plot([lastbin(nn)       lastbin(nn)],[0 max(cutocc)],'k--')
            end
       occ_sm{nn} = smoothvect(occ{nn},gaussian(xkern,8*xkern));
            % ** Later, exclude bins with very low occupancy and replace with NaN
            lowocc = occ_sm{nn} < xdel*.1 ;  % Maggie used threshold of 0.1 s / cm    
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

       %%%%%%  Create encoding model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       spk_e{nn} = cell(1,numunits(nn));        % encoding spike times
       unitmap{nn} = cell(1,numunits(nn));      % place maps
       
       for cc = 1:numunits(nn)
           
           tet = detc{nn}(cc,3);
           cellnum = detc{nn}(cc,4);
                     
           if ~isempty(spikes{day}{epenc}{tet}{cellnum}.data)
               spiketimes = spikes{day}{epenc}{tet}{cellnum}.data(:,1);
           else
               continue
           end
           
           %  Spikes in experimenter-transcribed (in notebook) epoch
           inds_epenc =  (  spiketimes  >  starttime_epenc  )  &  ( spiketimes < endtime_epenc );
           %  Spikes in SWRs    Spikes that occur within SWR (2 SD) periods
           inds_swr = isExcluded(spiketimes, periodtimes_rip2)  ;
           
           % Determine choice of encoding spikes :  "model number"
           if modelnum == 1
               % basic case
               inds_e = inds_epenc;
               encodevec = ones(1,length(consvectimes_rip2));
           elseif modelnum == 2
               % exclusion case               
               % inds4 = isExcluded(filedata.params(:,1)/10000,gfperiods);
               inds_e = inds_epenc & inds_trajencode;
               %encodevec = []; % would need to code this 
               %disp(sprintf('%d of %d spikes excluded from encoding',sum(inds_d)-sum(inds_e),sum(inds_d)))
           elseif modelnum == 3
               % non-SWR
               inds_pre = inds_epenc;
               inds_e = inds_epenc & ~inds_swr;
                    % also, formulate encoding periods in a period format (necessary
                    % to create accurate occupancy map below)
                    encodeperiods{nn} = vec2list(~consvec_rip2,consvectimes_rip2);
               if 0
                    disp(sprintf('%d of %d spikes excluded from encoding',sum(inds_pre)-sum(inds_e),sum(inds_pre)))      
               end
           elseif modelnum == 4
               % head speed is > 10 cm/s
               inds_hispeed = logical(isExcluded(spiketimes,hispeedperiods));
               inds_pre = inds_epenc;
               inds_e = inds_epenc & inds_hispeed;
                    % also, formulate encoding periods in a period format (necessary
                    % to create accurate occupancy map below)
                    encodeperiods{nn} = vec2list(~consvec_rip2,consvectimes_rip2);
               if 0
                    disp(sprintf('%d of %d spikes excluded from encoding',sum(inds_pre)-sum(inds_e),sum(inds_pre)))      
               end               
           else
               keyboard
           end
           
           % Spike times (recording clock) for Encoding
           spk_e{nn}{cc} = spiketimes(inds_e);  
           
           % Calculate firing rate map for each unit
           posinds = lookup(spiketimes(inds_e),postimevec{nn});
                spk_armdists = armdists{nn}(posinds); 
                    spk_bincounts = histc(spk_armdists,xbins{nn});
                    spk_bincounts_sm = smoothvect(spk_bincounts,...
                                                  gaussian(xkern,8*xkern));
                         spk_bincounts_sm = spk_bincounts_sm(:);
           unitmap{nn}{cc} = spk_bincounts_sm ./ occ_sm{nn} + 0.05;   

       end
       
   end
   
   % To check, plot unit place maps
   if 0
       for yy = 2
           for hh = 1:length(unitmap{yy})
               figure ; plot(unitmap{yy}{hh})
               title(sprintf('%d %d %d %d',detc{yy}(hh,:)),'fontsize',14,'fontweight','bold')
               pause
               close all
           end
       end
   end
   
   
   %%% Second, collect decoding data %%%%%%%%%%%%%%%%%%%%%%%%   
   

   
   
   ep;  % the decoding epoch #
   
   clear linpos
                  
   % Decode for each encoding epoch  (since unit common lists might be different,
                            %         doing two separate encoding epochs)
   binvec = [];             % edges of decoding time bins
   binvec_c = [];           % centers of decoding time bins
   
   S = cell(1,2);             % [ <times> <bin #> <cellnumber> ]
   spikecounts = cell(1,2);   %  row: cell#   x  column: bin, spike count 
   cellcounts = cell(1,2);
   activebins = cell(1,2);

   % Set up bin time vector
   if isempty(manual_period)
       startdec = starttime(1);
       enddec = ceil(endtime(1)+winsize);
   else
       startdec = starttime(1) + manual_period(1);
       enddec = starttime(1) + manual_period(2);
   end
   binvec = startdec:winsize:enddec;
        numbins = length(binvec) - 1 ;               % number of decoding bins
        binvec_c = binvec(1:(end-1)) + winsize/2 ;   % bin center times
      
       % Go through each encode-decode list's respective unit list
        %   and then calculate the binned spike counts
       for nnnn = 1:num_encodeeps
           
           S{nnnn} = [];
           spikecounts{nnnn} = [];
           
           for cc2 = 1:numunits(nnnn)
               
               tet2 = detc{nnnn}(cc2,3);
               cellnum2 = detc{nnnn}(cc2,4);
               spk = spikes{day}{ep}{tet2}{cellnum2}.data;   % decode all clustered spikes
               if ~isempty(spk)
                   spk_d = spk(:,1);
               else
                   spk_d = [];
               end
               
               % Assign each spike into a bin
               [~,spikebins] = histc(spk_d, binvec);
               validinds = spikebins > 0;           % where index matrix is non-zero
               spikebins = spikebins(validinds);    % filter for spikes contained in the epoch
               numspk = length(spikebins);
               
               % Concatenate to celldata ( [ <spike times> <bin #> <cellnumber> ] ) and
               if numspk > 0
                   tmpcelldata = [ spk_d(validinds)   spikebins(:)   cc2*ones(numspk,1)  ];
                   S{nnnn} = [S{nnnn} ; tmpcelldata];
               end
               
               % Register spikecounts: # of spikes from each cell in each bin (cell # x bin #)
               spikecount_unit = zeros(1,numbins);
               for i = 1:length(spikebins)
                   spikecount_unit(spikebins(i)) = spikecount_unit(spikebins(i)) + 1;
               end
               
               spikecounts{nnnn} = [spikecounts{nnnn} ; spikecount_unit];
               
           end
           
           % Sort unit data
           S{nnnn} = sortrows(S{nnnn},1);
           
           % For each time bin, count # of distinct cells contributing at least one spike
           cellcounts{nnnn} = sum((spikecounts{nnnn} > 0),1);
           
           % Find all decoding bins with enough cells firing to count as 'active'
           activebins{nnnn} = find(cellcounts{nnnn} >= min_activecells);   % by bin index
           
       end

    %%% Sixth, decode (Bayes' Rule) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

    posteriors = cell(1,2);
    xbins_cut = cell(1,2);
    numlinpos = [nan nan];
    numcutlinpos = [nan nan];
    
    for nnn = 1:num_encodeeps
   
        % define these to make code cleaner
        numlinpos(nnn) = length(xbins{nnn});     % total # of spatial locations (traj-dist combination)
        
        % initialize output
        posteriors{nnn} = nan(numlinpos(nnn),numbins);
        
        for bb = 1:numbins
            
            % inds of spikes in this bin
            inds = ( S{nnn}(:,2) == bb ) ;
            
            % Identify active ("a") encoding cells in this bin
            cellnums_a = unique( S{nnn}(inds,3) )' ;
                numactive = length(cellnums_a);
                
            % If no clustered spiking in this bin, then posterior set to 0s
            if isempty(cellnums_a)
                posteriors{nnn}(:,bb) = zeros(numlinpos(nnn),1);   % just report a series of 0s if there are no active spikes this bin
                continue
            end
            
            % Iterate through each cell, calculating their individual likelihood functions
            L = nan(numactive,numlinpos(nnn)); 
            for vv = 1:numactive
                ccc = cellnums_a(vv);  % cell number in detc list
                bin_numspk = spikecounts{nnn}(ccc,bb);
                % iterate over spatial locations
                for pp = 1:numlinpos(nnn)
                    % retrieve estimated firing rate at that spatial location
                    fr = unitmap{nnn}{ccc}(pp) + 0.001;  % (adding a small value since actual Poisson dist is asymptotic)
                    lambda = fr * winsize;
                    % calculate the Poisson probability of firing the bin's # of spikes there
                    L(vv,pp) = poisscdf(bin_numspk,lambda) - poisscdf(bin_numspk-1,lambda);
                end
                % normalize each unit's likelihood function
                L(vv,:) = L(vv,:)/sum(L(vv,:));
            end
            % multiply the individual units' likelihoods together (independence assumption)
            LL = prod(L,1);
            % normalize to obtain posterior, then install
            posteriors{nnn}(:,bb) = LL' / sum(LL);
                  
            % troubleshoot
            if 0
                H = figure
                hold on
                % plot only the posterior + position data of the current traj
                plot(L(:,inds)','--','linewidth',2)  % individual likelihoods
                plot(posteriors(bb,inds),'-k','linewidth',3) % posterior
                pause
                close(H)
            end
            
            % print out progress
            if bb == round(numbins/100)
                disp('1%')
            elseif bb == round(numbins/10)
                disp('10%')
            elseif bb == round(numbins/2)
                disp('50%')
            elseif bb == round(numbins*.75)
                disp('75%')
            end
            
        end
        
        % eliminate invalid position bins %%%%%%%%%
       
        if 0
            posteriors{nnn}((d_cut{nnn}+1):(e_cut{nnn}-1),:) = [];
            posteriors{nnn}((b_cut{nnn}+1):(c_cut{nnn}-1),:) = [];       
        else
           posteriors{nnn} = posteriors{nnn}(validxbins{nnn},:); 
        end
        
        numcutlinpos(nnn) = size(posteriors{nnn},1);
        
        
    end
    
        % save output
        cd(savedir)
        disp('done calculating full epoch posteriors')
        savefilename = sprintf('%s_%d_%d',animalname(1:3),day,ep);
        P = struct;
        P.animalname = animalname;
        P.dayep = [day ep];
        P.remoteep = remoteep;
        P.detc = detc;          
        P.binvec_c = binvec_c;
        P.armdists = armdists_cat{1};   % actual animal positions
        P.linposbins{1} = (1:numcutlinpos(1)) * xdel ;
        if remoteW && length(xbins_cut) > 1
            P.linposbins{2} = (1:numcutlinpos(2)) * xdel ;
        end
        P.lastcenterbin = lastcenterbin;
        P.lastrightbin = lastrightbin;
        P.lastbin = numcutlinpos ;
        P.activebins = activebins;
        P.posteriors = posteriors;

        save(savefilename,'P')
        
        clear P
            
    
    
end


