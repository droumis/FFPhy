function kk_clusterlessdecode19(directoryname,fileprefix,dayeps,animalname,varargin)
  
% L vs. R (for Pro and Ret) clusterless decode w/ Uniform prior kk March 2018

animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander',...
                 'Iceland','Justice','Kapital','Laplace','Eight','Ten','Conley','Miles'} ;
animalinfo = animaldef(animalname);
   daydir = getdaydir(animalname);
   an = find(strcmp(animalname,animals_order));

% Default values   
epsilon = eps(0);   % infinitesimal numerical value
SKIP_SAVED_FILE = 0;
manual_time = [];
Excurbuffer = 0.5;  % in s     

% Options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'savedir'
            savedir = varargin{option+1};    
        case 'nummsbin'
            nummsbin = varargin{option+1};
        case 'nummsoverlap'
            nummsoverlap = varargin{option+1};
        case 'EXTYPE'
            EXTYPE = varargin{option+1};            
        case 'dt'
            dt = varargin{option+1};        
        case 'msig'
            msig = varargin{option+1};          
        case 'spikethresh'
            spikethresh = varargin{option+1};
        case 'SKIP_SAVED_FILE'
            SKIP_SAVED_FILE = varargin{option+1};                        
        case 'Priormodel'
            Priormodel = varargin{option+1};
        case 'manual_time'
            manual_time = varargin{option+1};
        case 'exclude_FS'
            exclude_FS = varargin{option+1};
        case 'Excurbuffer'
            Excurbuffer = varargin{option+1};
        case 'THETA_BINS'
            THETA_BINS = varargin{option+1};
    end
end

dt = nummsbin/1000; % postimevec step-size

% Identify day eps 
task = loaddatastruct(directoryname,fileprefix,'task');
    epochfilter = epochmaker('runW_rip');
        dayeps_all = evaluatefilter(task,epochfilter);
        if isempty(dayeps)
            dayeps = dayeps_all;
        end

% Iterate through day eps
for de = 1:size(dayeps,1)
    
    d = dayeps(de,1);
    ep = dayeps(de,2);

    % (optional) if saved, skip calculation
    if SKIP_SAVED_FILE
        cd(savedir)
        savefilename = sprintf('%s_%d_%d_LR_%d.mat',animalname(1:3),d,ep,EXTYPE);
        if ~isempty(dir(savefilename))
            disp(sprintf('>>>>> %s found, skipping',savefilename));
            continue
        end
    end
        
   % Identify the tetrodes to get spike data from %
        % In specific, identify tets w/ clustered non-N principal units
   load('/opt/data13/kkay/Unitgroups/Unitgroup_Principal_16_animals.mat','Principal_16_animals');
   adtc = Principal_16_animals;
   load('/opt/data13/kkay/Superlin_data/Classtable_09-Apr-2015.mat');
   adtc_ca2n = sortrows([ classtable6_adtc{2}{2} ; classtable6_adtc{2}{3} ],[1 2 3 4]);
   adtc_nonNunit = adtc(~ismember(adtc,adtc_ca2n,'rows'),:);
        anim_ind = adtc_nonNunit(:,1) == an;
        day_ind    = adtc_nonNunit(:,2) == d;
      selected_tets = unique(adtc_nonNunit(anim_ind & day_ind,3));
        selected_tets = selected_tets(:)';
      maxtet = max(selected_tets);
      if isempty(selected_tets)
          continue
      end

   % Linear position data %
   linpos = loaddatastruct(directoryname,fileprefix,'linpos',d);  
       % Time vector of position
       postimevec = linpos{d}{ep}.statematrix.time;
           epstart = postimevec(1);
           epend = postimevec(end);
       % Choice point position
       CP = choicepointer(linpos{d}{ep});   % CP's precise linear distance  
           maxC = ceil(CP);  % The 1-cm bin of the CP i.e. the last 1-cm bin of the Center arm
       % Animal's position
       lindist = linpos{d}{ep}.statematrix.lindist;
       segind = linpos{d}{ep}.statematrix.segmentIndex;
       % Identify W-track segment indices (in postimevec)
       seg1 = (segind == 1) ;
       seg2 = (segind == 2) ;
       seg3 = (segind == 3) ;
       seg4 = (segind == 4) ;
       seg5 = (segind == 5) ;
       
   % Moving periods (>4 cm/s)
   timefilterscript
   [dummy,~] = evaluatetimefilter(animalinfo{2},animalinfo{3},{ vel4 },[d ep]);
   moveperiods = dummy{d}{ep};
        movevec  = list2vec(moveperiods,postimevec);
        
   % Get LR data % %%%%%%%%%%%%%%%%%%%
    % Trajencode periods (1, 2, 3, 4) %
   trajencode = loaddatastruct(animalinfo{2},animalinfo{3},'trajencode',d);  % used to get center divider  
   right_outbound_vec = (trajencode{d}{ep}.trajstate == 1);  
   left_outbound_vec  = (trajencode{d}{ep}.trajstate == 3);       
   right_inbound_vec  = (trajencode{d}{ep}.trajstate == 2);
   left_inbound_vec   = (trajencode{d}{ep}.trajstate == 4);       
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
   
   % Initialize variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   
      % LR %
   LR = [];            % 1: LEFT, 2: RIGHT trajectory of the animal 
   LR2 = [];           % (same as above but row vector)
   xbins = [];          % position ("x") vector in 0,200,400 space (Dir mark space) 
   num_xbins = nan;     % total # of bins in Dir mark space

      % prior %
   T = [];           % Transition matrix (Markovian)

      % likelihood (encoding) %
   encodeperiods = [];              % [starttime endtime] of enc periods
   mark0 = cell(1,maxtet);          % {tet} Amp marks for enc spk
   Xnum = cell(1,maxtet);           % {tet} Dir Gauss kern for each enc spk
   Lint = cell(1,maxtet);           % {tet} Firing occnorm Dir map for each tet
   posind_spk_e = cell(1,maxtet);   % {tet} Dir inds of each enc spk
   numspk_e = nan(maxtet,1);        % {tet} # of enc spk on this tet
   Xnum_all = [];                   % (all enc spk) [ xbins x spk_e ] Gauss kern around the Dir of the rat during each enc spk (each column)
   Lint_all = [];                   % (all enc spk) [ xbins x 1 ] Adds spikes (each one probability 1) and divides by occupancy and pos time step
   
   
   %%% Preliminary A: Set-up Positional space (W-track linearization) %%%%%%%%%%%%%%%%%%%%%%%%      
       
   disp(['Encoding : Day ',num2str(d), ', Ep ',num2str(ep)])

   % Initialize Dir %%%%%%%%%% Dir: (W-track position linearized + concatenated (C, L, R)) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   LR  = nan(size(lindist)) ;
   LR2 = nan(size(lindist))';   

   % EXCURSION encoding periods, defined here
   if EXTYPE == 6           % Pro 
       Lind = left_outbound_vec;  % TRAJENCODE
       Rind = right_outbound_vec; % TRAJENCODE
   elseif EXTYPE == 7       % Ret
       Lind = left_inbound_vec;  % TRAJENCODE
       Rind = right_inbound_vec;  % TRAJENCODE
   else
       keyboard
   end

   LRvec = ( Lind | Rind ) ;  % Encoding bins (time bins)
   
   % Linearization concatenation to padded form 
   LR(Lind)           = -1;   % -1: LEFT  (left arm trajectory)
   LR(Rind)           = +1;   % +1: RIGHT (right arm trajectory) 
   LR(isnan(LR)) = nan;  % just a reminder that NaN means that spikes in this Dir bin was not used to encode!
   LR2 = LR(:)';

   % error checking for some early epochs
   if length(unique(LR(~isnan(LR)))) ~= 2
       disp('SKIPPING EPOCH SINCE LR NOT BOTH REPRESENTED')
        continue
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
   %%%% Construct Dir space  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xbins = [ -1 1 ];
        num_xbins = 2;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
            
  %%%%% Preliminary B: construct Uniform prior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  unifprior = [0.5 ; 0.5];                    % initialize
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      
  %%%% Preliminary C: construct Transition matrix (T) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  T = [nan nan ; nan nan];

  if Priormodel == 0   % uniform
      T = [0.5 0.5 ; 0.5 0.5 ];
  else
      keyboard
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
   %%%%%% Preliminary D: collect Encoding spikes' Marks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % Variables 
   mark0           ;    % {tet} Marks (4 spike amps) for enc spk
   numspk_e        ;    % {tet} Total # of enc spk on each tet
   posind_spk_e    ;    % {tet} Dir inds of enc spk

   % Iterate tetrodes
   for tet = selected_tets

       filedata = loadparamsfile(daydir,d,tet);  % spike data from day directory
       SPK = filedata.params;

       disp(filedata.paramnames{2})
       disp(filedata.paramnames{3})
       disp(filedata.paramnames{4})
       disp(filedata.paramnames{5})

       % Set-up encoding periods: moving + arm/dir selection 
       encodeperiods = vec2list(LRvec & movevec,postimevec);  
       %encodeperiods = vec2list(linmovevec & IOvec,postimevec);  

       % Filter for encoding spikes
       sinds_amp     = any( SPK(:,2:5) > spikethresh , 2) ;                  % spikes >thresh uV in at least one channel
       sinds_exc     = logical( isExcluded( SPK(:,1)/10000 ,encodeperiods));  % spikes occurring during excursion-specific encoding periods
       if THETA_BINS == 2       % if 1st half Theta cycle only, filter here
           sinds_1st_half     = logical( isExcluded( SPK(:,1)/10000 ,Theta1sthalfperiods));  % spikes occurring during 1st half of theta cycles
       else
           sinds_1st_half     = true(size(SPK(:,1)));
       end
           % (optional) Filter out FS spikes ("type") to get pyramidal neuron spikes
            if exclude_FS    
                   maxwidths = filedata.params(:,6)/10;
                   sinds_nonFS = maxwidths > .2;
               sinds_type = sinds_nonFS;    % filter
           else
               sinds_type = true(size(SPK(:,1)));
           end
       sinds_e      = sinds_type & sinds_amp & sinds_exc & sinds_1st_half;                   % final set of encoding spikes                

       % Get Mark values for Enc spikes
       disp(sprintf('%d encoding spk',sum(sinds_e)))
       spk_e{tet}           = SPK(sinds_e,1)/10000;                   % Enc spk times
            numspk_e(tet)   = length(spk_e{tet});                     % Total # of enc spk
       mark0{tet}           = SPK(sinds_e,2:5);                       % Enc spk amp mark  (4 channel amps) 
       posind_spk_e{tet}    = lookup(spk_e{tet},postimevec);          % Enc spk pos time indices

       clear filedata

   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%% Preliminary E: Collect Decoding spk (all spk) marks %%%%%%
        %%% also, construct Amplitude mark space %%%%%%%%%%%%%%%%%%%   
          
   % Initialize
   spktimes_d       = cell(1,maxtet);   % {tet} dec spk times
   posind_spike_d   = cell(1,maxtet);   % {tet} dec spk Dir inds
   markAll          = cell(1,maxtet);   % {tet} dec spk 5D mark (Dir + 4 amp)
   
   for tet2 = selected_tets
       
       filedata2 = loadparamsfile(daydir,d,tet2);  % spike data from day directory       
       SPK2 = filedata2.params;         % matrix of spike parameters
       
       % Obtain dec spk info
           sinds_ep  = ( SPK2(:,1)/10000  >=  epstart )  &  ( SPK2(:,1)/10000 <= epend );    % spikes in epoch
           sinds_amp = any( SPK2(:,2:5) > spikethresh, 2) ;                                    % spikes of sufficient amp
           if exclude_FS    
                   maxwidths = SPK2(:,6)/10;
                   sinds_nonFS = maxwidths > .2;
               sinds_type = sinds_nonFS;    % filter
           else
               sinds_type = true(size(SPK2(:,1)));
           end             
       spktimes_d{tet2}         =  SPK2(sinds_ep & sinds_amp & sinds_type,1)/10000    ;     % dec spk times  
       posind_spike_d{tet2}     =  lookup(spktimes_d{tet2},postimevec)   ;     % dec pos inds

       % Covariate - Mark of dec spk: [ Dir <amp1 amp2 amp3 amp4> ]
       markAll{tet2}(:,1)       = posind_spike_d{tet2} ;                             % Col 1:    Dir inds
       markAll{tet2}(:,2:5)     = SPK2(sinds_ep & sinds_amp & sinds_type,2:5);       % Col 2-5:  amps
       
       clear filedata
       
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%% Preliminary F: Identify Dec spikes parent tet + spike # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
   disp('Fourth, collect dec spike info')
  
   SPK_d                    = [];       % dec spk times
   SPK_d_sort               = [];       % dec spk, sorted
       spk_d_sort_parenttet     = [];       % dec spk parent tet, sorted
       spk_d_sort_tetspknum     = [];       % dec spk # on parent tet, sorted
   
   % Collect dec info across all tetrodes
   for tet = selected_tets
       numspktet        =  length(spktimes_d{tet});
       %                                    [ <time>           <parent_tet>              <spk # on parent tet>       <pos ind of spk> ]
       SPK_d        =  [SPK_d     ;      spktimes_d{tet}    tet * ones(numspktet,1)         (1:numspktet)'        posind_spike_d{tet} ];
   end
   [SPK_d_sort , sinds_d]    = sort(SPK_d(:,1));      % Sort in time   
       spk_d_sort_parenttet  = SPK_d(sinds_d,2);   % Parent tets
       spk_d_sort_tetspknum  = SPK_d(sinds_d,3);   % Spike # for parent tet
       
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
   
   %%% I. (Likelihood) Construct Dir Enc model w/ Gaussian KDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   disp('Third, construct Dir Encoding model')
   
   % Initialize
   Xnum;                % {tet} Dir Gauss kern for each enc spk
   Lint;                % {tet} Firing occnorm Dir map for each tet
   
   % Identify the rat's actual Dir (during enc periods)
   posinds_enc = [];
   numposbin_enc = nan;     
   posinds_enc = logical(isExcluded(postimevec,encodeperiods)); 
        numposbin_enc = sum(posinds_enc);
   
   % i. Calculate the occupancy map (rat's actual Dir)
   X1        = repmat(xbins',1,numposbin_enc);               % [ xbins x postimebins_enc ] 0,200,400 Dir vecs repeated as columns across enc pos time bins
   MU1       = ones(length(xbins),1)  *  LR2(posinds_enc);  % [ xbins x postimebins_enc ] Actual Dir of rat (for each enc time bin) in each column (every xbin is same value)
   deltapdf1 = ( X1 == MU1 ) ;                               % [ xbins x postimebins_enc ] Kronecker Delta kern at the Dir of the animal (each column)
   occ       = sum(deltapdf1,2);                              % [ xbins x 1         ] Rat's enc period Dir occupancy from Gauss kerns
        occ  = occ / sum(occ);                    
   
   % (tetrode - specific encoding)
   for tet3 = selected_tets
       
       % ii.  Xnum: Dir Gauss kern for each enc spk
           % (preliminary) collect each Enc spk's Dir
           X2  = repmat(xbins',1,numspk_e(tet3));                       % [ xbins x spk_e ]  0,200,400 Dir vecs, repeated for each enc spk (each column)
           MU2 = ones(length(xbins),1)  *  LR2(posind_spk_e{tet3});    % [ xbins x spk_e ]  Actual Dir of rat during each enc spk, repeated over (each column)
       Xnum{tet3} = ( X2 == MU2 );                              % [ xbins x spk_e ]  Gauss kern around the Dir of each enc spk (each column)
       
       % iii.  Lint: Firing occ-norm Dir map for this tet
       totalencdur = sum(encodeperiods(:,2) - encodeperiods(:,1));
       
       Lint{tet3} = sum(Xnum{tet3},2) ./ occ ./ totalencdur;            % [ xbins x 1 ]   Adds spikes (each one probability 1) and divides by occ (enc periods) and pos time step               
       
   end

   clear deltapdf1 X1 MU1 X2 MU2
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %%% II. Decode w/ Bayes' Rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    % i. (preliminary) Set-up decoding time bins %%%%
    
    Tvec = [];          % Centers of dec time bins (zeroed to ep start time)
    Win = [];          % Windows of dec time bins (zeroed to ep start time)
    Thetahalf = [];     % (optional) 1st or 2nd
    
    if THETA_BINS == 0          % Moving Window decode
        
        if nummsoverlap == 0       % Non-overlapping dec windows
            Tvec = (epstart:(nummsbin/1000):epend) - epstart;
            Win = [Tvec(1:(end-1))' Tvec(2:end)'] ;
            Tvec = mean(Win,2);  % Re-define time vector to bin centers
        elseif nummsoverlap > 0    % Overlapping dec windows
            Tvec  = ( (epstart + nummsbin/2000) : (nummsoverlap/1000) : (epend + nummsbin/2000) ) - epstart;
            numcenterbins = length(Tvec);
            Win =   nan(numcenterbins,2) ;
            for gg = 1:size(Win,1)
                Win(gg,:) = Tvec(gg) + [-nummsbin/2000  +nummsbin/2000] ;
            end
            Tvec = mean(Win,2);  % Re-define time vector to bin centers
        end
        
    elseif THETA_BINS == 1 || THETA_BINS == 2      % Theta-bin based decode
        
        % Tvec: [ center times of each theta half ] ;
        % Wins: [ win times of each theta half    ] ;
        
        % 1. Load 
        thetabins       = loaddatastruct(animalinfo{2},animalinfo{3},'thetabins',d);
            tph         = thetabins{d}.tph{ep};
            tvec_tph    = thetabins{d}.timevec_tph{ep};
            tbins       = thetabins{d}.thetabins{ep};      % tbins of full cycles
            
        % 2. Subdivide cycles into 1st & 2nd

            % i. Subdivide cycles into halves (pi/2 window) using Thetabinner
            Tbins               = cell(1,4);  % {1}: 1st, {2}: 2nd, {3}: full, {4}: 1st-2nd
                Tbins           = Thetabinner(tbins,tph,tvec_tph);
                Thalves_mid     = mean( Tbins{4} ,2);
                    Tvec        = Thalves_mid - epstart;   % mid time of each Theta half
                    if length(unique(Tvec)) ~= length(Tvec)
                        keyboard
                    end
                    Win         = Tbins{4} - epstart;      % [start end] of each Theta half
    
            % ii. Construct Thetahalf: indicator vector for 1st vs. 2nd
                Numcychalves = length(Tvec);        
            Thetahalf = nan(Numcychalves,1);
                Thetahalf( mod(1:Numcychalves,2) == 1 ) = 1;  % 1st half 
                Thetahalf( mod(1:Numcychalves,2) == 0 ) = 2;  % 2nd half
          
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ii. (preliminary) Identify Excursions (i.e. contiguous periods to decode) %%%%%%%%%%%%%%%%
    Excurlist  = [];  % [ starttime  endtime  excur_flag  percentage_traj ]
    Numexcur   = [];  
    excurlist_descript = '';
    if ~isempty(manual_time)  
        disp('Decoding manual time')
        Excurlist = [manual_time(1) manual_time(2) nan] + epstart;
        Numexcur = 1;
    else  
    disp('Decoding excursions')
        excursions          = loaddatastruct(animalinfo{2},animalinfo{3},'excursions',d);
        excurlist           = excursions{d}{ep}.excurlist_movefrags;
        excurlist_descript  = excursions{d}{ep}.excurlist_movefrags_cols;
        
        % Based on EXTYPE, select by excursion flag values
        if EXTYPE == 6
            exflags = [1 3 -11];              % Outbound trajs (1,3) + Home well trackbacks (-11)           
        elseif EXTYPE == 7
            exflags = [2 4 -11];              % Inbound trajs (2,40) + Home well trackbacks (-11) 
        end
        excurlist_raw   = excurlist(ismember(excurlist(:,3),exflags),:);  % select subset of excursions
        
        % For each EXTYPE, low-speed time buffer before and after
        %     and after doing so, "glue" excursions of the same type together 
        Excurlist = [];  % final output
        for XT = exflags
            % Identify specific excursion type
            excurlist_xt = excurlist_raw( excurlist_raw(:,3) == XT ,[1 2]);
                % Buffer at beginning and end
                excurlist_xt(:,1) = excurlist_xt(:,1) - Excurbuffer;
                excurlist_xt(:,2) = excurlist_xt(:,2) + Excurbuffer;
            % Glue buffered excursions of the same type together!
                mstimevec = epstart:0.001:epend;
            excurvec_xt = false(size(mstimevec));  % 1-ms indicator vector
            for exx = 1:size(excurlist_xt,1)
                ind_a = lookup(excurlist_xt(exx,1),mstimevec);
                ind_b = lookup(excurlist_xt(exx,2),mstimevec);
                excurvec_xt(ind_a:ind_b) = 1;
            end
            % Collect
            excurlist_xt = vec2list(excurvec_xt,mstimevec);
                numexcur_xt = size(excurlist_xt,1);
            Excurlist = [  Excurlist ;   excurlist_xt    XT * ones(numexcur_xt,1)  ]  ;
        end
        Excurlist = sortrows(Excurlist,[1]);   % sort in time
            Numexcur    = size(Excurlist,1);
        %keyboard    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % iii. Decode each excursion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Tvecs       = cell(Numexcur,1); % Tvecs for each excursion
    Wins        = cell(Numexcur,1); % Wins for each excursion
    Thetahalves = cell(Numexcur,1); % (optional) Thetahalf indicator vector for each excursion
    Posts       = cell(Numexcur,1); % Decodes for each excursion
    
    for ex = 1:Numexcur
        
        % Excursion info
        exstart     = Excurlist(ex,1);
        exend       = Excurlist(ex,2);
            a_dec = lookup(exstart - epstart,Tvec);  % First dec bin
            b_dec = lookup(exend - epstart,Tvec);    % Last dec bin
                numdecbin = b_dec - a_dec + 1;             % # of decoding bins in this excursion
        
        % Basic info %
        Tvecs{ex}       = Tvec(a_dec:b_dec);
        Wins{ex}        = Win(a_dec:b_dec,:);
        if THETA_BINS == 1 || THETA_BINS == 2
            Thetahalves{ex} = Thetahalf(a_dec:b_dec);
        end    
        % Initialize output %
        Posts{ex}       = zeros(num_xbins,numdecbin) ;              % Matrix of posteriors                 
        
        %%% Initialize posterior with uniform prior
        %postx = unifprior;
        
        tic
        % Iterate each dec time bin
        for bbb = 1:numdecbin
            
            t = bbb - 1 + a_dec;                % Tvec index
            decbin = Win(t,:) + epstart;   % [starttime endtime]
            midtime = mean(Win(t,:),2);    % [ midtime ]
            
                disp(sprintf('(LR %d) %s d%d ep%d: %s',EXTYPE,animalname(1:3),d,ep,num2str(midtime)))
            
            %%% Prior %%%%
            prior = unifprior;
            
            %%% Likelihood (across all tetrodes) %%%
                
            % Spikes that occurred w/in this dec bin
            sinds = find( isExcluded( SPK_d_sort, decbin ) ) ;
                parenttets = spk_d_sort_parenttet(sinds);
            
            % Initialize matrix of likelihoods from each tet ( [xbins x tetnum] )
            L_tets = zeros(num_xbins,length(selected_tets)) ;

            % Iterate to get likelihood from each tet
            for ttt = 1:length(selected_tets)
                
                tet = selected_tets(ttt);
              
                % Spikes from this tet in this bin
                sinds_tet       = sinds(parenttets == tet) ;
                    numspk_tet  = length(sinds_tet);
                
                if numspk_tet == 0
                    
                    L_tets(:,ttt) = exp(-Lint{tet} .* dt) ;   % % Eq'n 2.5 (Deng et. al. 2015)
                    
                elseif numspk_tet > 0
                    
                    % Initialize matrix of single-spike likelihoods
                    l = zeros(num_xbins,numspk_tet) ;   % [ xbins x spk ]
                    
                    % Iterate each dec spk
                    for s = 1:numspk_tet
                        
                        % Spike info
                        ind = sinds_tet(s);                % Sorted spike index
                        i = spk_d_sort_tetspknum(ind) ;    % Spike # on parent tet
                        
                        % l0: Amp mark likelihoods for this spike
                        amp      = markAll{tet}(i,2:5);
                        amp_enc  = mark0{tet};
                        enc_1vec = ones(1,numspk_e(tet));
                        
                        l0 = normpdf(   amp(1) * enc_1vec,...             % x:      amp  *   [ 1 x numspk_enc ]    % dec spk amplitude
                                        amp_enc(:,1)',...                 % mu:     amp_e    [ 1 x numspk_enc ]    % enc spk amplitudes
                                        msig ) .* ...                     % smker:  amp bandwidth                  % Gaus kern SD
                                        normpdf( amp(2) * enc_1vec,...
                                        amp_enc(:,2)',...
                                        msig ) .* ...
                                        normpdf( amp(3) * enc_1vec,...
                                        amp_enc(:,3)',...
                                        msig ) .* ...
                                        normpdf( amp(4) * enc_1vec,...
                                        amp_enc(:,4)',...
                                        msig ) ;
                        
                        % Eq'n 2.5 (Deng et. al. 2015)
                        l1 = Xnum{tet} * l0' ./ occ ./ totalencdur ;      % small lambda      ( Xnum is enc spk Dir Gaus kern )
                        l2 = ( l1 * dt ) .* exp( -Lint{tet} .* dt ) ;     % full likelihood
                        l2 = l2 + epsilon;                                % add infinitesimal value (prevents normalization error)
                        l2 = l2 ./ sum(l2(isfinite(l2))) ;                % normalize
                        l(:,s) = l2 ;
                        
                    end
                    
                    % Tetrode likelihood: product of single-spike likelihoods
                    L_tets(:,ttt) = prod(l,2) ;                             
                    
                    L_tets(:,ttt) = L_tets(:,ttt) / sum(L_tets(:,ttt)) ;      % normalize
                    
                    
                    
                end
                
            end
                
            L = prod(L_tets,2);    % Full likelihood: product of tetrode likelihoods
            
            if any(isnan(L))
                keyboard
            end
            
            Post        = prior .*  L ;         % Bayes' rule
                Post    = Post / sum(Post) ;  % normalize
            
            Posts{ex}(:,bbb) = Post;
            clear postx l L;
        
            % Plot manual decode
            if ~isempty(manual_time)
                if midtime >= manual_time(2)
                    decodeper = vec2list(sum(Posts{ex},1) > 0,1:numdecbin);
                    a = decodeper(1);
                    b = decodeper(end);
                    plot_start  = Tvecs{ex}(a);
                    plot_end    = Tvecs{ex}(b);
                    aa = lookup(plot_start,postimevec - epstart);
                    bb = lookup(plot_end,postimevec - epstart);
                    F = figure('units','normalized','outerposition',[.5 .3 .7 .3]);
                    
                    % LR
                    plot(postimevec(aa:bb) - epstart,1.02 * LR(aa:bb),'color',[.1 .8 .1],'linewidth',3); hold on
                    xlim([plot_start plot_end]);
                    
                    % Decode
                    plot(Tvecs{ex}(a:b),diff(Posts{ex}(:,a:b),[],1),'k-','linewidth',2); hold on
                    imagesc(Tvecs{ex}(a:b),2,diff(Posts{ex}(:,a:b),[],1),[-1 1]); colormap(redblue)                    

                    set(gca,'ydir','normal')
                    xlim([plot_start plot_end]);                    
                    ylim([-1.1 +3])
                    str1 = sprintf('%s %d %d',animalname(1:3),d,ep);
                    str2 = sprintf('Ha: %d uv, bin: %d ms, overlap: %d ms',msig,nummsbin,nummsoverlap);
                    title({str1,str2},'fontsize',14,'fontweight','bold');
                    keyboard
                end
            end

        end
        
       
    end
    
    % Save output in a file
    disp('Done calculating epoch')
    clear P
    P.date              = date;
    P.animalname        = animalname;
    P.dayep             = [d ep];
    P.EXTYPE            = EXTYPE;
    P.THETA_BINS        = THETA_BINS;
    P.divider           = '%%%%%%%%%%%%%%%%';    
    P.spikethresh       = spikethresh;
    P.nummsbin          = nummsbin;
    P.nummsoverlap      = nummsoverlap;
    P.selected_tets     = selected_tets;
    P.msig              = msig;
    P.divider1           = '%%%%%%%%%%%%%%%%';
    P.postimevec = postimevec;
    P.posvec = LR;                      % Actual animal position state
    P.yvec_pos = [-1 1];             
    P.divider2           = '%%%%%%%%%%%%%%%%';   
    P.excurlist          = Excurlist;
    P.excurlist_descript  = excurlist_descript;
    P.divider3           = '%%%%%%%%%%%%%%%%';    
    P.Priormodel        = Priormodel;   % Flag for which prior used
    P.Posts             = Posts;        % {excur} Posteriors 
    P.Tvecs             = Tvecs;        % {excur} Time vectors 
    P.Wins              = Wins;         % {excur} Time windows 
    P.Thetahalves       = Thetahalves;  % (optional) 1st vs. 2nd half theta cycle {excur}
    cd(savedir);
    savefilename = sprintf('%s_%d_%d_LR_%d',animalname(1:3),d,ep,EXTYPE);
    save(savefilename,'P','-v7.3')
    clear P Posts Tvecs
    
end

end








