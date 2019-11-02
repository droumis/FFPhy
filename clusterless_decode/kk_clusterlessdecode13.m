function kk_clusterlessdecode12(directoryname,fileprefix,dayeps,animalname,varargin)
  
% outpos 

animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander',...
                 'Iceland','Justice','Kapital','Laplace','Eight','Ten','Conley','Miles'} ;
animalinfo = animaldef(animalname);
   daydir = getdaydir(animalname);
blank = '';
tetfilter = '';
savedir = '/opt/data13/kkay/Superposteriors_data';
decodemode = 1;  % default is all epoch
modelnum = 1;
remoteW = 0;
spikethresh = 0;
mdel = 10;  % in uV
xdel = 1;  % in cm
dt = .001;
sigma_randomwalk = 1;
plot_infunction = 0;
TETSET = 5;
nummsbin = 1; 5;
dec_start = 0; 72;
epsi = eps(0);
SKIPSAVED = 0;
CELLMAX = 0;
Priormodel = 1;  % 1: Random walk, 2: Empirical
DCoffset = 1e-9;
manualstoptime = [];

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'savedir'
            savedir = varargin{option+1};    
        case 'nummsbin'
            nummsbin = varargin{option+1};
        case 'dec_start'
            dec_start = varargin{option+1};
        case 'plot_infunction'
            plot_infunction = varargin{option+1};
        case 'decodemode'
            decodemode = varargin{option+1};            
        case 'modelnum'
            modelnum = varargin{option+1};  
        case 'remoteW'
            remoteW = varargin{option+1};
        case 'dt'
            dt = varargin{option+1};        
        case 'mdel'
            mdel = varargin{option+1};
        case 'xdel'
            xdel = varargin{option+1}; 
        case 'sxker'
            sxker = varargin{option+1}; 
        case 'smker'
            smker = varargin{option+1};  
        case 'sigma_randomwalk'
            sigma_randomwalk = varargin{option+1};            
        case 'TETSET'
            TETSET = varargin{option+1};
        case 'spikethresh'
            spikethresh = varargin{option+1};
        case 'SKIPSAVED'
            SKIPSAVED = varargin{option+1};                  
        case 'CELLMAX'
            CELLMAX = varargin{option+1};       
        case 'DCoffset'
            DCoffset = varargin{option+1};
        case 'Priormodel'
            Priormodel = varargin{option+1};
        case 'manualstoptime'
            manualstoptime = varargin{option+1};
    end
end

task = loaddatastruct(directoryname,fileprefix,'task');
    epochfilter = epochmaker('runW_rip');

% identify day eps
dayeps_all = evaluatefilter(task,epochfilter);
if isempty(dayeps)
    dayeps = dayeps_all;
end
clear task

% Iterate through epochs to decode
for de = size(dayeps,1):-1:1
    
    d = dayeps(de,1);
    ep = dayeps(de,2);  % "local" ep
    
   % Load data
   linpos = loaddatastruct(directoryname,fileprefix,'linpos',d);
   pos = loaddatastruct(directoryname,fileprefix,'pos',d);
   tetinfo = loaddatastruct(directoryname,fileprefix,'tetinfo',d);  
        
   % Identify tetrodes
   if TETSET == 5
       disp('selecting tets that have clustered non-N principal units')
       % CA1, CA2 P, CA3 units
       load('/opt/data13/kkay/Unitgroups/Unitgroup_Principal_16_animals.mat','Principal_16_animals');
       adtc = Principal_16_animals;
       % ii. identify non CA2 N units
       load('/opt/data13/kkay/Superlin_data/Classtable_09-Apr-2015.mat');
       an = find(strcmp(animalname,animals_order));
       
            if SKIPSAVED
                cd('/opt/data50/kkay/__Decode')
               filename = sprintf('%s_%d_%d_fullepoch_1_outpos.mat',animalname(1:3),d,ep);
               if ~isempty(dir(filename))
                   disp(sprintf('********** %s found, skipping',filename));
                  continue 
               end
            end
       
       adtc_ca2n = sortrows([ classtable6_adtc{2}{2} ; classtable6_adtc{2}{3} ],[1 2 3 4]);
       adtc_nonN = adtc(~ismember(adtc,adtc_ca2n,'rows'),:);
            animalind = adtc_nonN(:,1) == an;
            
            % only day with max cells
            if CELLMAX
                adtc_list = adtc_nonN(animalind,:);
                days = unique(adtc_list(:,2))';
                maxday = 0; cellmax = 0;
                for dy = days
                    numcellsday = sum(adtc_list(:,2) == dy);
                    if numcellsday >= cellmax
                        maxday = dy;
                        cellmax = numcellsday;
                    end
                end
                if d ~= maxday
                    disp(sprintf('(CELLMAX) Ignoring day %d : %d cells',d,sum(adtc_list(:,2) == d)));
                    continue
                else
                    disp(sprintf('(CELLMAX) Choosing day %d : %d cells',d,cellmax));
                end
            end
            
            dayind    = adtc_nonN(:,2) == d;
          selected_tets = unique(adtc_nonN(animalind & dayind,3));
            selected_tets = selected_tets(:)';
          maxtet = max(selected_tets);
    end           
       
   % Initialize outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
      % positional variables
   xbins = [];
   armdists = [];        % two elements since up to two different W-tracks
   centerarmmax = nan;
   rightarmmax = nan;
      % prior model %
   stateM = [];
      % likelihood variables %
   mark0 = [];
   Xnum = [];
   Lint = [];
   numspikes_e = [];
   posind_spike_e = [];
      % internal use (not to save in output file later)
   num_xbins = nan;
   armdists2 = []; armdists_cat = [];
      C_cut = []; L_cut = [];R_cut = [];

  encodeperiods = [];    

   %%% First, iterate through each encoding epoch to collect basic data %%%%%%%%%%%%%%%%%%%%%%%%      
       
       disp(['Encoding : Day ',num2str(d), ', Epoch ',num2str(ep)])
       
       % Basic epoch data
       postimevec = linpos{d}{ep}.statematrix.time;
            epstart = postimevec(1);
            epend = postimevec(end);
       
       % *********at some point use posvecmaker.m to do the below code     
            
       % Identify W-track segment indices
       lindist = linpos{d}{ep}.statematrix.lindist;
       seg1 = linpos{d}{ep}.statematrix.segmentIndex==1;
       seg2 = linpos{d}{ep}.statematrix.segmentIndex==2;
       seg3 = linpos{d}{ep}.statematrix.segmentIndex==3;
       seg4 = linpos{d}{ep}.statematrix.segmentIndex==4;
       seg5 = linpos{d}{ep}.statematrix.segmentIndex==5;
       CPbuff = choicepointer(linpos{d}{ep});
          maxC = ceil(CPbuff);
          
       % Headdir
       headdir = linpos{d}{ep}.statematrix.segmentHeadDirection(:,1);   % head direction relative to center well -- values > 0 are outbound
            postimevec_nonnan = postimevec(~isnan(headdir));
            headdir_nonnan = headdir(~isnan(headdir));
            headdir2 = interp1(postimevec_nonnan,headdir_nonnan,postimevec,'linear');          
       
       % Initialize pospos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       io  = nan(size(lindist)) ;   % outpos: W-track linearized (and concatenated, with R on top) position
       io2 = nan(size(lindist)); 
       
                outinds = (headdir2 >= 0);   % 
                ininds =  (headdir2 <  0);   %             

                % (sanity check) Plot binary in/out w/ head direction 
                if 1
                    sttime = pos{d}{ep}.data(1,1);
                    figure;       hold on
                        plot(pos{d}{ep}.data(outinds,1)-sttime,headdir2(outinds),'.','markersize',10,'Color',[1 0 0]);   % out: red
                        plot(pos{d}{ep}.data(ininds,1)-sttime,headdir2(ininds),'.','markersize',10,'Color',[0 0 1]);     %  in: blue                   
                    keyboard
                end  
        
                    % Outbound vs. Inbound, Center-Left-Right
                    io( ininds )    =   1;     % 1: in 
                    io( outinds )   =   2;    % 2: out  
                    
                io = io(:)   ;
                io2 = io(:)' ;            
                                   
       %%%% Construct positional mark space  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       xbins =  [1 2] ; %1 : xdel : maxR;    % bin vector
       num_xbins = 2;
                                
      %%%%% Uniform prior (used at the start of every decode) %%%%%%%%%%%
      uniformprior = [0.5 ; 0.5];
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

      %  Construct the transition matrix from animal's behavior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      stateM;
      
       %%%%%%  Create encoding model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       % initialize epoch cell, will fill with tetrode data
       
       mark0 = cell(1,maxtet);
       numspikes_e = nan(1,maxtet);
       posind_spike_e = cell(1,maxtet);
       
       for tet = selected_tets
           
           filedata = loadparamsfile(daydir,d,tet);  % spike data from day directory
           
           disp(filedata.paramnames{2})
           disp(filedata.paramnames{3})
           disp(filedata.paramnames{4})
           disp(filedata.paramnames{5})
           
           % Identify spikes for encoding 
           
           %  Spikes in experimenter-transcribed (in notebook) epoch
           inds_epenc =  (  filedata.params(:,1)/10000  >=  epstart  )  &  ( filedata.params(:,1)/10000 <= epend );
           %  Spikes that has spikethresh (uV) voltage amplitude in at least one channel
           inds_thresh = any(filedata.params(:,2:5) > spikethresh, 2) ;
           %  Spikes at >4 cm/s and out-bound heading
           timefilterscript
           [dummy,~] = evaluatetimefilter(animalinfo{2},animalinfo{3},{ vel12 },[d ep]);
                vel4_periods = dummy{d}{ep};
                %outdir_periods = vec2list((headdir2 >= 0),postimevec);
           inds_moving = isExcluded(filedata.params(:,1)/10000,vel4_periods)  ;
           %inds_outdir = isExcluded(filedata.params(:,1)/10000,outdir_periods)  ;
           
           % encoding periods
           vel4vec  = list2vec(vel4_periods,postimevec);
           %outvec   = list2vec(outdir_periods,postimevec);
           encodeperiods = vec2list(vel4vec,postimevec);
           
           % Determine choice of encoding spikes :  "model number"
           if modelnum == 1
               inds_e = inds_epenc & inds_thresh & inds_moving;  
           elseif modelnum == 2
               % identify non-FS spikes
               maxwidths = filedata.params(:,6)/10; 
               inds_nonFS = maxwidths > .2;
               % spikes
               inds_e = inds_epenc & inds_thresh & inds_moving  & inds_nonFS;              
           else
               keyboard
           end
           
            disp(sprintf('%d spikes',sum(inds_e)))             

           % Spike times (recording clock) for E and D
           spktimes_e{tet} = filedata.params(inds_e,1)/10000;    % exclusive, to ENCODE
                numspikes_e(tet) = length(spktimes_e{tet});
                
           % Amplitude mark vector (4 channel amplitudes) 
           mark0{tet} = filedata.params(inds_e,2:5);
           
           % Obtain position indices for each spike 
           posbin = (postimevec(2)-postimevec(1));
           centerpostimes = [postimevec ; postimevec(end) + posbin]  - posbin/2;
           [        ~ ,        posind_spike_e{tet} ]  =  histc(spktimes_e{tet},centerpostimes)  ;    % ENCODING
                
           clear filedata
       end

   
   %%% Second, for decoding epoch, collect amplitude mark data %%%%%%%%%%%%%%%%%%%%%%%%   
   
   ep;  % the decoding epoch #
   
   % Basic epoch data
   postimevec;
   epstart;
   epend;

   clear linpos
                  
   % Collect decoding spike data
   spktimes_d = cell(1,maxtet);
   posind_spike_d = cell(1,maxtet);
   markAll = cell(1,maxtet);
   minamp = inf;
   maxamp = -inf;
   for tet2 = selected_tets
       filedata = loadparamsfile(daydir,d,tet2);  % spike data from day directory       
       inds_ep =  (  filedata.params(:,1)/10000  >=  epstart  )  &  ( filedata.params(:,1)/10000 <= epend );
       inds_thresh = any(filedata.params(:,2:5) > spikethresh, 2) ;
       spktimes_d{tet2} = filedata.params(inds_ep & inds_thresh,1)/10000;
            centerpostimes2 = [postimevec ; postimevec(end)+posbin] - posbin/2;  % centers
       [ ~ , posind_spike_d{tet2} ]  =  histc(spktimes_d{tet2},centerpostimes2)   ;         
       markAll{tet2}(:,1) = posind_spike_d{tet2} ;         % Col 1: position time bin # , Col 2-5: spike amplitudes
       markAll{tet2}(:,2:5) = filedata.params(inds_ep & inds_thresh,2:5);
       % also, track of max and min spike amplitudes in decoding spike set
       % + in the remote-W encoding spike set
       minamp_tet_local = min(min(markAll{tet2}(:,2:5)));
       maxamp_tet_local = max(max(markAll{tet2}(:,2:5)));
       minamp_tet =  minamp_tet_local;
       maxamp_tet =  maxamp_tet_local ;
       if minamp_tet < minamp
           minamp = minamp_tet;
       end
       if maxamp_tet > maxamp
           maxamp = maxamp_tet;
       end
       clear filedata
   end

   % Set up the Amplitude mark space
   if isempty(spikethresh)
        mbins =  minamp : mdel : maxamp;
   else
        mbins = spikethresh : mdel : maxamp;
   end
   
   %%% Third, construct encoding model (Gaussian KDE) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   disp('Third, construct encoding model')
   
   Xnum = cell(1,maxtet);
   Lint = cell(1,maxtet);
   
   % need these for accurate occupancy spatial maps
   posinds_enc = [];
   numpossamps_enc = nan;  
   
   % actual animal positions in the encoding epoch
   % but to do correctly, need to only look at encoding times
   posinds_enc = logical(isExcluded(postimevec,encodeperiods));
   numpossamps_enc = sum(posinds_enc);
   X1  = repmat(xbins',1,numpossamps_enc);                % position bin vectors repeated as columns
   MU1 = ones(length(xbins),1)  *  io2(posinds_enc);  % actual linearized (+ concatenated) positions of the animal
   
   % 1.  occ: occupancy (linear position) map
   normpdf1 = normpdf(X1,MU1,sxker);
   occ = repmat(sum(normpdf1,2),1,length(mbins));  %normpdf1 * ones(numpossamps,length(mbins))   ;
   occ = occ/sum(occ(:,1));
   
   for tet3 = selected_tets
       
       % encoding spike positions
       X2  = repmat(xbins',1,numspikes_e(tet3));                                % position vectors repeated as columns
       MU2 = ones(length(xbins),1)  *  io2(posind_spike_e{tet3});    % actual (encoding) spike positions
       
       % 2.  Xnum: Gaussian kernel estimators for position, for each encoding spike
       Xnum{tet3} = normpdf(X2,MU2,sxker);
       
       % normalize for each spike
       %for sk = 1:size(Xnum{tet3},2)
       %    Xnum{tet3}(:,sk) = Xnum{tet3}(:,sk) / sum(Xnum{tet3}(:,sk));
       %end
       %keyboard
       
       if 1
           % kill nan spikes (can happen at end if dir pos not assigned)
           deadspikes = isnan(Xnum{tet3}(1,:));
           % delete spikes
           Xnum{tet3} (:,deadspikes) = [];
           posind_spike_e{tet3}(deadspikes) = [];
           mark0{tet3}(deadspikes,:) = [];
           % subtract spike count
           num_deadspikes = sum(deadspikes);
           numspikes_e(tet3) = numspikes_e(tet3) - num_deadspikes;
           if num_deadspikes > 0
               disp('killing deadspikes')
               if num_deadspikes > 10
                   keyboard
               end
           end
       end
       
       % 3.  Lint: firing occ-norm position map for tetrode
       Lint{tet3} = sum(Xnum{tet3},2) ./ occ(:,1) ./ dt;   %integral (?)
       % set non-occupied xbins indices to 0
       Lint{tet3} = Lint{tet3} ./ sum(Lint{tet3}) ;  % normalization
       
   end

  
   clear normpdf1 X1 MU1 X2 MU2
   
   
   %%% Fourth, collect decoding (all) spikes in each epoch
 
   disp('Fourth, collect decoding (all) spikes in each epoch')
   
        % variables to keep track of for both local and remote epochs
   spktimes_all         = [];       % all spike times in each epoch
   spktimes_all_sorted  = [];       % sorted spike times, sorted
   
   % Collect all spikes
   for tet = selected_tets
       numspikestet_all        =  length(spktimes_d{tet});
       %                     [  <spike time>     <parent tet #>                  <spike # on parent tet>    <pos inds of spikes>]
       spktimes_all        =  [spktimes_all     ;   spktimes_d{tet}   tet * ones(numspikestet_all,1)    (1:numspikestet_all)'    posind_spike_d{tet} ];
   end
   % for decoding epoch only, store this
   [spktimes_all_sorted , sinds] = sort(spktimes_all(:,1));
   num_d_spikes                      = length(spktimes_all_sorted);
       
     
   % Also construct a spike-tetrode matrix for the local (decoding) epoch
       % i.e. an "indicator matrix" : [  <spike #>  x  <which tetrode spikes> ]
      % initialize
  
   spike_d_parenttet = nan(num_d_spikes,1);
   spike_d_tetspikenum = nan(num_d_spikes,1);
       spike_d_parenttet = spktimes_all(sinds,2);
       spike_d_tetspikenum  = spktimes_all(sinds,3);
       
       clear spktimes_all
       
%    if decodemode > 1
%        disp(sprintf('calc tet_ind: %d spikes',length(spktimes_all{1})))
%        tic
%        for i = 1:length(spktimes_all{1})
%            parenttet = spktimes_all{1}(sinds(i),2) ;   % tet responsible for the spike
%            tet_ind(i,parenttet) = 1 ;
%        end
%        toc
%        tet_spikenum = tet_ind.*cumsum(tet_ind,1);  % [ <spike #>  x <index of spike for its parent tetrode>]
%    end
   
   
        
    
     
       
    %%% Fifth, obtain likelihood of not observing a spike  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('Fifth, collect all epikes in each epoch')    
    
    Lint_all = [];
    occ = [];    % save this for later
    
    Xnum_all = [];
    occ = nan(1,length(xbins));
    
    % occ: behavioral occupancy map in encoding epoch
    occ = normpdf( xbins'*ones(1,numpossamps_enc),...
        repmat(io2(posinds_enc),num_xbins,1),...
        sxker )   *    ones(numpossamps_enc,length(mbins) )    ;
    % remove invalid rows and normalize
    occ = occ/sum(occ(:,1));
    
    % Xnum: Gaussian kernel estimators for position
    if 1
        num_allspikes_e = 0;
        posind_allspikes_e = [];
        for vv = 1:maxtet
            if ~isnan(numspikes_e(vv))
                num_allspikes_e = num_allspikes_e + numspikes_e(vv);
                posind_allspikes_e = [posind_allspikes_e ; posind_spike_e{vv}(:)];
            end
        end
        Xnum_all = normpdf( xbins'*ones(1,num_allspikes_e),...           %  # of encoding spikes this epoch
            repmat(io2(posind_allspikes_e),num_xbins,1),...    %  # pos inds of all encoding spikes this epoch
            sxker );
    elseif 0
        Xnum_all = normpdf( xbins'*ones(1,length(spktimes_all)),...
            repmat(pospos2(posind_spike_all),num_xbins,1),...
            sxker );
    end
    
        if 1
            % kill nan spikes (can happen at end if dir pos not assigned)
            deadspikes = isnan(Xnum_all(1,:));
            Xnum_all(:,deadspikes) = [];
            num_deadspikes = sum(deadspikes);
            if num_deadspikes > 0
                disp('killing deadspikes')
                if num_deadspikes > 10
                    keyboard
                end
            end
        end
            
    
    % set non-occupied xbins indices to 0
    %Xnum_all(~validxbins,:) = 0;
    % normalize for each spike
    %Xnum_all = bsxfun(@rdivide,Xnum_all,sum(Xnum_all,1));
    %for skk = 1:size(Xnum_all,2)
    %    Xnum_all(:,skk) = Xnum_all(:,skk) / sum(Xnum_all(:,skk));
    %end
    % Lint: conditional intensity function for the unmarked case
    Lint_all = sum(Xnum_all,2)  ./  occ(:,1) ./  dt ; %integral
    % normalize
    Lint_all = Lint_all ./ sum(Lint_all)  ;


    clear Xnum_all
    
    %%% Sixth, decode (Bayes' Rule) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

    posteriors = [];    
    
    % Set up 1-ms epoch vector
    timevec_dec_ms = round(epstart*1000):nummsbin:round(epend*1000);   % 1-ms time vector (from position time)
        timevec_dec = timevec_dec_ms / 1000;
    numdecbin = length(timevec_dec_ms);
    eventperiods = [];  % initialize

    if decodemode == 1              % 1: full epoch
        disp('Decoding full epoch continuously')
        eventperiods = [0 timevec_dec(end)]; 
        numevents = 1;
        eventtimes = [];
    end

    posteriors{1} = [];
    posteriors{2} = [];
    
    % Initialize millisecond time vector

    % Iterate through 7 SD ripples (or events)
    for ww = 1:numevents
        
        % ms time of event (e.g. ripple) start
        evstart  = eventperiods(ww,1);
        
        postx = uniformprior;

        % Set up Prior
        Tparm = 3;
        if Tparm == 0
            small = 0.5;  % uniform
        else
            small = 10^(-Tparm);
        end
        Tmat = [(1 - small) small ; small (1-small)];
        
        % Initialize output matrices %%%%%%%%%%%%%%%%%%%%%%
            
            % initialize
            postxM_r = zeros(num_xbins,numdecbin) ;              % Matrix of posteriors
            SPK = spktimes_all_sorted;                   % Spike times in terms of 1-ms epoch vector
            
            tic
            % Iterate through each 1-ms time bin of the event
            ind_start = floor(dec_start * 1000 / nummsbin);     
            
            for t = ind_start:numdecbin
                    
                % Prior
                prior = Tmat * postx;  % one step prediction, Markovian predication
                               
                % Likelihood
                % initialized as flat
                L = ones( num_xbins, 1  );
                % obtain indices of spikes that occur closest to this 1-ms (epoch) time bin
                decbin = timevec_dec( t:(t + 1) ) - 0.5 * nummsbin / 1000;
                spind = find( isExcluded( SPK, decbin ) ) ;
                
                num_spikesinbin = length(spind);
                
                if num_spikesinbin == 0
                    
                    L = exp(-Lint_all .* dt) ;   %exp(-Lint_all.*dt) ;
                    
                elseif num_spikesinbin > 0
                    
                    % initialize
                    l = zeros(num_xbins , num_spikesinbin) ;
                    
                    % iterate through spikes (from all tetrodes) contained in this bin
                    for s = 1:num_spikesinbin
                        ind = spind(s);    % spike index in spktimes_all_sorted
                        if 0
                            % old indexing approach (matrix)
                            tet = find(tet_ind(ind,:)) ;
                            i = tet_spikenum(ind,tet);   % "running" spike # for its parent tetrode
                        else
                            % new indexing approach (vector)
                            tet = spike_d_parenttet(ind) ;
                            i = spike_d_tetspikenum(ind) ;
                        end
                        %l0: amplitude mark likelihoods for this single spike
                        l0 = normpdf(markAll{tet}(i,2)*ones(1,numspikes_e(tet)),...  % ch. 1 amp (scalar)  * ones vector length of all spikes
                                     mark0{tet}(:,1)',...                            % ch. 1 amplitude of all spikes
                                     smker).*...
                            normpdf(markAll{tet}(i,3)*ones(1,numspikes_e(tet)),...   % ch. 2 " "
                                     mark0{tet}(:,2)',...
                                     smker).*...
                            normpdf(markAll{tet}(i,4)*ones(1,numspikes_e(tet)),...
                                    mark0{tet}(:,3)',...
                                    smker).*...
                            normpdf(markAll{tet}(i,5)*ones(1,numspikes_e(tet)),...
                                    mark0{tet}(:,4)',...
                                    smker) ;    
                        l1 = Xnum{tet} * l0' ./ occ(:,1) ./ dt ;   % divide by occupancy map
                        l2 = l1 * dt .* exp(-Lint{tet} .* dt) ;  % exp factor
                            l2 = l2 + epsi;
                        l2 = l2 ./ sum(l2(isfinite(l2))) ;                       % normalize the posterior
                        l(:,s) = l2 ;
                    end
                    
                    L = prod(l,2) ; % take product of likelihoods across spikes
                    L = L./sum(L) ;
                    disp(sprintf('(outpos) %s d%d ep%d: %s',animalname(1:3),d,ep,num2str(t * nummsbin / 1000)))
                end
                
                postx = prior .*  L ;
                postx = postx / sum(postx) ;  % normalize
                
                postxM_r(:,t) = postx;
                clear onestep a l L;
                
                %if t == 306031
                %    keyboard
                %end
                
                % Print updates if doing full epoch / save halves
                if decodemode == 1
                    if t == round(1*numdecbin/100)
                        disp('1% done')
                    elseif t == round(5*numdecbin/100)
                        toc
                        disp('5% done')
                    elseif t == round(10*numdecbin/100)
                        disp('10% done')
                    elseif t == round(20*numdecbin/100)
                        disp('20% done')
                    elseif t == round(23*numdecbin/100)
                        disp('21% done')                
                    elseif t == round(50*numdecbin/100)
                        disp('50% done')
                    elseif t == round(90*numdecbin/100)
                        disp('90% done')
                    end
                    
                    % manual time stop
                    if ~isempty(manualstoptime)
                    
                    if t >= manualstoptime*1000 / nummsbin
                        tvec = (1:numdecbin)*nummsbin/1000;
                        decodeper = vec2list(sum(postxM_r,1) > 0,1:numdecbin);
                        plot_start = decodeper(1)*nummsbin/1000;
                        plot_end = decodeper(end)*nummsbin/1000;
                            a = lookup(plot_start,tvec);
                            b = lookup(plot_end,tvec);
                        aa = lookup(plot_start,postimevec - epstart);
                        bb = lookup(plot_end,postimevec - epstart);
                         figure('units','normalized','outerposition',[.5 .3 .7 .3])
                         dirtrace = postxM_r(1,a:b) - postxM_r(2,a:b); hold on
                         plot(tvec(a:b),dirtrace,'k-','linewidth',2); hold on
                         %colormap hot; 
                        colormap gray ; colormap(flipud(colormap));
                        set(gca,'ydir','normal')
                        plot(postimevec(aa:bb) - epstart,2*(io(aa:bb)-1.5),'color',[0 .8 0],'linewidth',2);
                        str1 = sprintf('sd: %d, dc: %d',sigma_randomwalk,DCoffset);
                        str2 = sprintf('Ha: %d (uv), Hx: %d (cm)',smker,sxker);
                        title({str1,str2},'fontsize',14,'fontweight','bold')
                        keyboard
                    end
                    end
                end
                
            end
            
            % (optional) manual plotting code (in function)
            if 0
                % POSITION
                figure;
                plot(postimevec_epenc - postimevec_epenc(1),outpos,'Color',[.8 .8 .8],'linewidth',2);
            end         
            if 0
                % DECODE
                tvec = (1:num_msbin_steps)*nummsbin/1000;
                decodeper = vec2list(sum(postxM_r,1) > 0,1:num_msbin_steps);
                plot_start = decodeper(1)*nummsbin/1000;
                plot_end = decodeper(end)*nummsbin/1000;
                a = lookup(plot_start,tvec);
                b = lookup(plot_end,tvec);
                aa = lookup(plot_start,postimevec - epstart);
                bb = lookup(plot_end,postimevec - epstart);
                figure; imagesc(tvec(a:b),xbins,postxM_r(:,a:b),[0 .025]); hold on
                colormap gray
                colormap(flipud(gray))
                set(gca,'ydir','normal')
                %plot(postimevec(aa:bb) - epstart,dirpos(aa:bb),'r-','linewidth',2);
                keyboard
            end            
                        
                        
            % if full epoch, save single continuous posterior to file
            % before proceeding to remote epoch
            keyboard
            cd(savedir)
            disp('done calculating full epoch posterior')
            savefilename = sprintf('%s_%d_%d_%s_%d_outpos',animalname(1:3),d,ep,'fullepoch');
            P = struct;
            P.nummsbin = nummsbin;
            P.animalname = animalname;
            P.dayep = [d ep];
            P.remoteep = remoteep;
            P.selected_tets = selected_tets;
            P.spikethresh = spikethresh;
            P.posvec = io;   % actual animal traj / positions
            P.stateM_gausnorm = stateM;
            P.posteriors = postxM_r;
            save(savefilename,'P','-v7.3')
            clear P postxM_r_cut
        
    end
    

    
end
        
        
      






end