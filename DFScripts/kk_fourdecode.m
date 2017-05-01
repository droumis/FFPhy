function kk_fourdecode(directoryname,fileprefix,dayeps,animalname,varargin)
%         
animalinfo = animaldef(animalname);
   daydir = animalinfo{2};
savedir = '/opt/data13/kkay/___Superfourcode_data/Decode';
plot_infunction = 0;
min_activecells = 3;
place_thresh_fr = 2;
place_thresh_size = 6;
%winsize = .125;  % obsolete
numCodes = 5;

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'savedir'
            savedir = varargin{option+1};    
        %case 'winsize'
        %    winsize = varargin{option+1};  
        case 'remoteW'
            remoteW = varargin{option+1};
        case 'sigma_transmat'
            sigma_transmat = varargin{option+1};            
        case 'cellchoice'
            cellchoice = varargin{option+1};
    end
end

task = loaddatastruct(directoryname,fileprefix,'task');

epochfilter = epochmaker('runW_rip');

% Identify all day eps
dayeps_all = evaluatefilter(task,epochfilter);

if isempty(dayeps)
    dayeps = dayeps_all;
end

clear task

% Iterate through epochs to decode
for de = 1:size(dayeps,1)
    
    d = dayeps(de,1);
    ep = dayeps(de,2);  % "local" ep
    
    disp(sprintf('%s day %d ep %d',animalname(1:3),d,ep))
    
    clear E D
    
    % Load data
    task = loaddatastruct(directoryname,fileprefix,'task',d);
    linpos = loaddatastruct(directoryname,fileprefix,'linpos',d);
    pos = loaddatastruct(directoryname,fileprefix,'pos',d);
    cellinfo = loaddatastruct(directoryname,fileprefix,'cellinfo',d);
    spikes =  loaddatastruct(directoryname,fileprefix,'spikes',d);
        superlinfields_dir = '/opt/data13/kkay/Superlin_data/';
        filename = dir(sprintf('%sSuperlin_%s.*',superlinfields_dir,animalname(1:3)));
        load([superlinfields_dir filename(end).name],'superlin')
    
   % (if specified) Identify the remote W-track epoch # 
   if remoteW
      task = loaddatastruct(directoryname,fileprefix,'task',d);
      eps = dayeps_all(dayeps_all(:,1)==d,2)';
      currenv = task{d}{ep}.environment;
      remoteeps = [];
      for xx = eps
         if ~strcmp(currenv,task{d}{xx}.environment)
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
   epstart = [nan nan];
   epend = [nan nan]; 
   for zz = 1:num_encodeeps
        ep2 = epochs_enc(zz);
        postimevec{zz} = linpos{d}{ep2}.statematrix.time;
        numpossamps(zz) = length(linpos{d}{ep2}.statematrix.time);
        epstart(zz) = postimevec{zz}(1);
        epend(zz) = postimevec{zz}(end);
   end   
   
   % Identify cells to encode-decode with (for local and remote W-track epochs)
   detc{1} = [];
   detc{2} = [];
   
   if cellchoice == 0
   
       keyboard
   
   elseif cellchoice == 1
       
       animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander',...
           'Iceland','Justice','Kapital','Laplace','Eight','Ten','Conley','Miles'};
       an = find(strcmp(animalname,animals_order));
       
       % 1: Only CA1, CA2 P, CA3 units, and in addition checks whether place-active
       %        for each encoding epoch (LOCAL and, if available, REMOTE W-track epoch),
       
       load('/opt/data13/kkay/Unitgroups/Unitgroup_Principal_16_animals.mat','Principal_16_animals');
       load('/opt/data13/kkay/Unitgroups/Unitgroup_Principal_CA1_16_animals.mat','Principal_CA1_16_animals');

       % identify those active in the present epoch
       ADTC = Principal_16_animals;
            selectinds = (ADTC(:,1) == an) & (ADTC(:,2) == d)  ;
       adtc_all = ADTC(selectinds,[1 2 3 4]);
            
       % ii. exclude CA2 N units
       if 0
           load('/opt/data13/kkay/Superlin_data/Classtable_09-Apr-2015.mat');
           ca2n_adtc = sortrows([ classtable6_adtc{2}{2} ; classtable6_adtc{2}{3} ],[1 2 3 4]);
           ca2n_dtc = ca2n_adtc(ca2n_adtc(:,1) == an,2:4);
           if ~isempty(ca2n_dtc)
               for qq = size(adtc,1):-1:1
                   if rowfind(adtc(qq,[1 3 4]),ca2n_dtc)
                       disp('ca2n unit ignored')
                       adtc(qq,:) = [];
                   end
               end
           end
       end
       
      
       % Now, for the LOCAL and REMOTE W-track epochs, respectively,
       %    identify the subset of these clustered units that were place-active
       
       detc{1} = cell(1,numCodes);  % {1:loc / 2:rem}{code}
       detc{2} = cell(1,numCodes);  % 
       
       for ee = 1:num_encodeeps 
           
           ep2 = epochs_enc(ee);
           
           for cc2 = 1:size(adtc_all,1)
               
               tet = adtc_all(cc2,3);
               cellnum = adtc_all(cc2,4);
               
               PF_flag = 0;
               
               SLSTATE = 2;
               [ind_sl,regnum] = superlinfind([d ep2 tet cellnum],superlin,SLSTATE);
               
               if ind_sl == 0
                   disp('************fourdecode: cannot find this unit************')
                   continue
               end
                        
               dataentry = superlin.data{regnum}.detc_output{SLSTATE}(ind_sl);
               
               % Analyze pair types
               for Code = 1:numCodes  % Spatial type
                   
                   cendir_flag = 0;
                   
                   if Code == 1           % Pro Full
                       trajs = [3 1];       % [L,R]
                   elseif Code == 2       % Ret Full
                       trajs = [4 2];       % [L,R]
                   elseif Code == 3       % R Dir
                       trajs = [1 2];       % [O,I]
                   elseif Code == 4       % L Dir
                       trajs = [3 4];       % [O,I]
                   elseif Code == 5       % Cen Dir
                       trajs = [];         
                       cendir_flag = 1;     % [O,I] in Center arm
                   end
                   
                   % Place field detection
                   PF_flag = [ 0 0 ];  % flag if detect at least 1 place field
                                       % below, will require this to analyze the epoch
                   
                   if ~cendir_flag    % Conventional Spatial pairs
                       
                       for ttt = 1:2  % trajs # 1 and # 2 (S = 1 or 2:  L vs. R // S = 3 or 4 OR  In vs Out)
                           
                           traj = trajs(ttt);
                           
                           if traj <= length(dataentry.trajdata) && ...
                                   ~isempty(dataentry.trajdata{traj})
                               
                               datamat = dataentry.trajdata{traj};   % 1st col: location, 2nd col: raw occ, 3rd col: raw spk count, 5: smoothed place map
                               xvec = datamat(:,1);
                               occ = datamat(:,2);
                               spkcnt = datamat(:,3);
                               placemap = datamat(:,5);
                               
                               % Place field detection
                               placevec = placemap > place_thresh_fr;  % 2 Hz minimum
                               placelist = vec2list(placevec,xvec);  % start and end indices of each place field
                               if ~isempty(placelist)
                                   placefieldsizes = placelist(:,2) - placelist(:,1);  % size of each field
                                   if any(placefieldsizes >= place_thresh_size)  % 6 cm minimum
                                       PF_flag(ttt) = 1;  % flag that p.f. detected
                                   end
                               end
                           else
                               % if unavailable data for a traj, then
                               % disregard entirely
                               PF_flag = [0 0];
                               break
                           end
                           
                       end
                       
                   elseif cendir_flag == 1   % Center arm Directional
                       
                        for oi = 1:2  % Out, In
                           
                            if oi == 1          % 1: OUT
                               trajs = [1 3]; 
                            elseif oi == 2      % 2: IN
                               trajs = [2 4];
                            end
                            
                            for traj = trajs
                                
                                if traj <= length(dataentry.trajdata) && ...
                                        ~isempty(dataentry.trajdata{traj})
                                    
                                    datamat = dataentry.trajdata{traj};   % 1st col: location, 2nd col: raw occ, 3rd col: raw spk count, 5: smoothed place map
                                    xvec = datamat(:,1);
                                    occ = datamat(:,2);
                                    spkcnt = datamat(:,3);
                                    placemap = datamat(:,5);
                                    
                                    % Identify choice point (CP) linear distance
                                    CPbuff = choicepointer(linpos{d}{ep2});
                                    
                                    % Place field detection  (here, only in center arm)
                                    centervec = (xvec < CPbuff);  
                                    placevec = (placemap > place_thresh_fr) & centervec;  % 2 Hz minimum
                                    placelist = vec2list(placevec,xvec);  % start and end indices of each place field
                                    if ~isempty(placelist)
                                        placefieldsizes = placelist(:,2) - placelist(:,1);  % size of each field
                                        if any(placefieldsizes >= place_thresh_size)  % 6 cm minimum
                                            PF_flag(oi) = 1;  % flag that p.f. detected
                                        end
                                    end
                                    
                                else
                                    % if unavailable data for a traj, then
                                    % disregard entirely
                                    PF_flag = [0 0];
                                    break
                                end
                                
                            end
                           
                       end                
                       
                   end
                   
                   % if place-active, add to corresponding unit list
                   if any(PF_flag)
                       detc{ee}{Code} = [detc{ee}{Code} ; d ep2 tet cellnum ];  % if place-active, then add to list
                   end
                   
               end
        
           end
       end
       
       numclust = size(adtc_all,1);
       numunits{1} = [size(detc{1}{1},1)  size(detc{1}{2},1)  size(detc{1}{3},1)  size(detc{1}{4},1)  size(detc{1}{5},1)];
       numunits{2} = [size(detc{2}{1},1)  size(detc{2}{2},1)  size(detc{2}{3},1)  size(detc{2}{4},1)  size(detc{2}{5},1)];
       
       disp(sprintf('%d clust units: %d %d %d %d W1 place-active, %d %d %d %d W2 place-active',numclust,numunits{1},numunits{2}))
       
   end
       
   
   % Initialize outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   PM = {}; % {nn}{Code}   TWO BIN (L and R) PLACE MAP
  
   %%% Iterate through each encoding epoch to collect basic data %%%%%%%%%%%%%%%%%%%%%%%%      
   for nn = 1:num_encodeeps
       
       epenc = epochs_enc(nn);
       
       % Identify choice point (CP) linear distance
       CPbuff = choicepointer(linpos{d}{epenc});
            
       disp(['Encoding : Day ',num2str(d), ', Epoch ',num2str(epenc)])
       
       %%%%%%  Create encoding model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       for C = 1:numCodes  % CODE:  Pro, Retro, Right, Left
          
           cendir_flag = 0;
           
           if C == 1           % Pro Full
               trajs = [3 1];       % [L,R]
           elseif C == 2       % Ret Full
               trajs = [4 2];       % [L,R]
           elseif C == 3       % R Dir
               trajs = [1 2];       % [O,I]
           elseif C == 4       % L Dir
               trajs = [3 4];       % [O,I]
           elseif C == 5       % Cen Dir
               trajs = [];
               cendir_flag = 1;     % [O,I] in Center arm
           end
           
            PM{nn}{C} = nan(numunits{nn}(C),2);      % initialize
           
            for cc = 1:numunits{nn}(C)
           
                tet = detc{nn}{C}(cc,3);
                cellnum = detc{nn}{C}(cc,4);
                
                % 3. Spatial index %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                SLSTATE = 2;
                [ind_sl,regnum] = superlinfind([d ep tet cellnum],superlin,SLSTATE);
                if ind_sl == 0
                    disp('CANNOT FIND IN SUPERLIN (fourdecode)')
                    continue
                end
                
                dataentry = superlin.data{regnum}.detc_output{SLSTATE}(ind_sl);
                
                occnorm = [nan nan];  % epi 1, epi 2
                
                if ~cendir_flag
                    
                    for tt = 1:2   % episode pair    PRO RETRO:   1: L, 2: R  for     RIGHT LEFT:  1: OUT 2:IN
                        
                        tr = trajs(tt);
                        
                        if tr <= length(dataentry.trajdata) && ...
                                ~isempty(dataentry.trajdata{tr})
                            
                            datamat = dataentry.trajdata{tr};   % 1st col: location, 2nd col: raw occ, 3rd col: raw spk count, 5: smoothed place map
                            
                            % what we want    
                            occ = datamat(:,2);
                            spkcnt = datamat(:,3);
                            
                            % Tabulation of spikes and time occupied
                            totalspk = sum(spkcnt);
                            totaltime = sum(occ);
                            occnorm(tt) = totalspk / totaltime;
                        else
                            disp('this should not happen')
                            keyboard
                        end
                    end
                
                elseif cendir_flag == 1
                    
                     for oi = 1:2   % cendir (here episode pair is outbound inbound) 
                        
                         if oi == 1
                             trajs = [1 3];   % outbounds
                         elseif oi == 2
                             trajs = [2 4];   % inbounds
                         end
                         
                         totalspk = 0; 
                         totaltime = 0; 
                         
                         for tr = trajs
                             
                             if tr <= length(dataentry.trajdata) && ...
                                     ~isempty(dataentry.trajdata{tr})
                                 
                                 datamat = dataentry.trajdata{tr};   % 1st col: location, 2nd col: raw occ, 3rd col: raw spk count, 5: smoothed place map
                                    xvec = datamat(:,1);
                                 
                                 % Since Cen Dir: get Center-arm data 
                                 pos_a = lookup(0,xvec);
                                 pos_b = lookup(CPbuff,xvec);
                                 occ = datamat(pos_a:pos_b,2);
                                 spkcnt = datamat(pos_a:pos_b,3);
                                    totalspk = totalspk + sum(spkcnt);
                                    totaltime = totaltime + sum(occ);
                                 
                             else
                                 disp('this should not happen')
                                 keyboard
                             end
                             
                         end

                         % 1-bin place field
                         occnorm(oi) = totalspk / totaltime;
                         
                     end
                    
                end
                
                % Spatial index is simply occupancy normalized firing rate
                PM{nn}{C}(cc,:) = occnorm ;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end
                       
       end
            
   end

   
   %%% Second, collect decoding data %%%%%%%%%%%%%%%%%%%%%%%%   
  
   ep;  % the decoding epoch #
   
   clear linpos
                  
   % Decode for each encoding epoch  (since unit common lists might be different,
                            %         doing two separate encoding epochs)
                            
    %%% Load thetabins (previously calculated) %%%%%%%%%%%%%%%%5
    tbins = loaddatastruct(directoryname,fileprefix,'thetabins',d);
        thetatet = tbins{d}.thetatet;
        thetabins = tbins{d}.thetabins{ep};
        thetabins_c = mean(thetabins,2);
        contigflag = tbins{d}.contigflag{ep};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S = {};              % {ee}{C} [ <times> <bin #> <cellnumber> ]
    spikecounts = {};    % {ee}{C}     row: cell#   x  column: bin, spike count
    cellcounts = {};     % {ee}{C}
    activebins = {};     % {ee}{C}

   % Set up decoding bins
   if 0
       %    startdec = epstart(1);
       %    enddec = ceil(epend(1) + winsize);
       %    binvec = startdec:winsize:enddec;
       %    binvec_c = binvec(1:(end-1)) + winsize/2 ;   % bin center times
   elseif 1
       binlist      = thetabins ;             % edges of decoding time bins
       binlist_c    = thetabins_c ;           % centers of decoding time bins
       numbins      = size(binlist,1) ;       % number of decoding bins
       bindurs      = thetabins(:,2) - thetabins(:,1);
   end
   
        % Go through each encode-decode list's respective unit list
        %   and then calculate the binned spike counts
        for nnnn = 1:num_encodeeps
            
            for CC = 1:numCodes   % Code

                % initialize
                S{nnnn}{CC}             = [];    % [ spktimes binnums cellnum ]
                spikecounts{nnnn}{CC}   = [];    % thetabin vector w/ # of spikes from the cell
                
                for u = 1:numunits{nnnn}(CC)  % units
                    
                    tet2 = detc{nnnn}{CC}(u,3);
                    cellnum2 = detc{nnnn}{CC}(u,4);
                    if cellnum2 < length(spikes{d}{ep}{tet2}) ...
                            &&  ~isempty(spikes{d}{ep}{tet2}{cellnum2})
                        spk = spikes{d}{ep}{tet2}{cellnum2}.data;   % decode all clustered spikes
                        if ~isempty(spk)
                            spk_d = spk(:,1);
                        else
                            spk_d = [];
                        end
                    else
                       spk_d = [];  % this can happen if cell not clustered in current epoch but had PF in remote
                    end
                    
                    % Assign each spike into a bin
                    
                    if 1
                        % Reverse procedure for finding each spike's
                        % decoding bin
                        %tic
                        numspk0 = length(spk_d);
                        spikebins = nan(numspk0,1);
                        for bn = 1:numbins
                            spk_inds = logical(isExcluded(spk_d,binlist(bn,:)));
                            spikebins(spk_inds) = bn;
                        end
                        %toc
                        % remove spikes that weren't binned
                            validinds = ~isnan(spikebins);
                        spikebins(~validinds) = [];
                    elseif 0
                        % OBSOLETE %
                        [~,spikebins] = histc(spk_d, binvec);
                        validinds = spikebins > 0;           % where index matrix is non-zero
                        spikebins(~validinds) = [];    % filter for spikes contained in the epoch
                    end
                    
                    % Concatenate to celldata ( [  <spike times>  <bin #>  <cellnumber>  ] )
                    numspk = length(spikebins);
                    if numspk > 0
                        tmpcelldata = [ spk_d(validinds)   spikebins(:)   u*ones(numspk,1)  ];
                        S{nnnn}{CC} = [S{nnnn}{CC} ; tmpcelldata];
                    end
                    
                    % Register spikecounts: # of spikes from each cell in each bin (cell # x bin #)
                    spikecount_unit = zeros(1,numbins);
                    for i = 1:length(spikebins)
                        spikecount_unit(spikebins(i)) = spikecount_unit(spikebins(i)) + 1;
                    end
                    
                    spikecounts{nnnn}{CC} = [spikecounts{nnnn}{CC} ; spikecount_unit];
                    
                end
                
                if ~isempty(S{nnnn}{CC})
                
                    % Sort unit data
                    S{nnnn}{CC} = sortrows(S{nnnn}{CC},1);
                    
                    % For each time bin, count # of distinct cells contributing at least one spike
                    cellcounts{nnnn}{CC} = sum((spikecounts{nnnn}{CC} > 0),1);
                    
                    % Find all decoding bins with enough cells firing to count as 'active'
                    activebins{nnnn}{CC} = find(cellcounts{nnnn}{CC} >= min_activecells);   % by bin index
                    
                else
                    
                    S{nnnn}{CC} = [];
                    cellcounts{nnnn}{CC} = [];
                    activebins{nnnn}{CC} = [];
                    
                end
            end
        end
        
        
    %%% Sixth, decode (Bayes' Rule) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

    posteriors = {};
    
    for nnn = 1:num_encodeeps
     
        for CCC = 1:numCodes

            % initialize output
            posteriors{nnn}{CCC} = nan(numbins,2);
            
            if isempty(S{nnn}{CCC})
                continue
            end
            
            for bb = 1:numbins
                
                % Inds of spikes in this bin
                inds = ( S{nnn}{CCC}(:,2) == bb ) ;
                
                % Duration of the bin
                bindur = bindurs(bb);
                
                % Identify active ("a") encoding cells in this bin
                cellnums_a = unique( S{nnn}{CCC}(inds,3) )' ;
                numactive = length(cellnums_a);
                
                % If no clustered spiking in this bin, then posterior set to Nan
                if isempty(cellnums_a)
                    posteriors{nnn}{CCC}(bb,:) = [nan nan];   % just report a series of 0s if there are no active spikes this bin
                    continue
                end
                
                % Iterate through each cell, calculating their individual likelihood functions
                L = nan(numactive,2);
                for vv = 1:numactive
                    uuu = cellnums_a(vv);  % cell number in detc list
                    bin_numspk = spikecounts{nnn}{CCC}(uuu,bb);
                    % iterate over spatial indices
                    for pp = 1:2
                        % retrieve estimated firing rate at that spatial location
                        fr = PM{nnn}{CCC}(uuu,pp) + 0.001;  % (adding a small value since actual Poisson dist is asymptotic)
                        lambda = fr * bindur;
                        % calculate the Poisson probability of firing the bin's # of spikes there
                        L(vv,pp) = poisscdf(bin_numspk,lambda) - poisscdf(bin_numspk-1,lambda);
                    end
                    % normalize each unit's likelihood function
                    L(vv,:) = L(vv,:)/sum(L(vv,:));
                end
                % multiply the individual units' likelihoods together (independence assumption)
                LL = prod(L,1);
                % normalize to obtain posterior, then install
                posteriors{nnn}{CCC}(bb,:) = LL' / sum(LL);
                
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
        end
            
    end
        
    disp('done calculating full epoch posteriors')

% save output
cd(savedir)
savefilename = sprintf('%s_%d_%d',animalname(1:3),d,ep);
P = struct;
P.date = date;
P.animalname = animalname;
P.dayep = [d ep];
P.remoteep = remoteep;
P.detc = detc;
% time stuff
P.thetatet = thetatet;
P.binlist = binlist;
P.binvec_c = binlist_c;
P.numbins = numbins;
P.bindurs = bindurs;
P.activebins = activebins;
% decoded output
P.divider = '%%%%%%%%%%%%%%';
P.posteriors = posteriors;
%P.winsize = winsize;
P.min_activecells = min_activecells;
P.place_thresh_fr = place_thresh_fr;
P.place_thresh_size = place_thresh_size;
P.descript = 'second {} cell is remote ep, if available';

save(savefilename,'P')

clear P
    
    
end

    
end


