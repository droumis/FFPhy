function kk_clusterlessdecode17(directoryname,fileprefix,dayeps,animalname,varargin)
  
% Pos (positional) clusterless decode w/ Uniform prior kk February 2018

animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander',...
                 'Iceland','Justice','Kapital','Laplace','Eight','Ten','Conley','Miles'} ;
animalinfo = animaldef(animalname);
   daydir = getdaydir(animalname);
   an = find(strcmp(animalname,animals_order));

% Default values   
epsilon = eps(0);   % infinitesimal numerical value
SKIP_SAVED_FILE = 0;
manual_time = [];
xdel = 1;
Excurbuffer = 0.5;  % in s
EXCISE = 1;        

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
        case 'xdel'
            xdel = varargin{option+1}; 
        case 'xsig'
            xsig = varargin{option+1}; 
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
        case 'EXCISE'
            EXCISE = varargin{option+1};
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
for de = size(dayeps,1):-1:1
    
    d = dayeps(de,1);
    ep = dayeps(de,2);

    % (optional) if saved, skip calculation
    if SKIP_SAVED_FILE
        cd(savedir)
        savefilename = sprintf('%s_%d_%d_Pos_%d.mat',animalname(1:3),d,ep,EXTYPE);
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
       
   %  Moving periods (>4 cm/s)
   timefilterscript
   [dummy,~] = evaluatetimefilter(animalinfo{2},animalinfo{3},{ vel4 },[d ep]);
   moveperiods = dummy{d}{ep};
        movevec  = list2vec(moveperiods,postimevec);
        
   %  Moving periods (>4 cm/s)
   timefilterscript
   [dummy,~] = evaluatetimefilter(animalinfo{2},animalinfo{3},{ linvel20 },[d ep]);
   linmoveperiods = dummy{d}{ep};
        linmovevec  = list2vec(linmoveperiods,postimevec);
        

   % Head direction data %
   headdir = linpos{d}{ep}.statematrix.segmentHeadDirection(:,1);   % head direction relative to center well -- values > 0 are outbound
        postimevec_nonnan = postimevec(~isnan(headdir));
        headdir_nonnan = headdir(~isnan(headdir));
        headdir2 = interp1(postimevec_nonnan,headdir_nonnan,postimevec,'linear');                
        outvec = (headdir2 >= 0);
        invec  = (headdir2 < 0);       
 
   % Initialize variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   
      % pos %
   Pos = [];         % Linearized concatenated (0,200,400) position of the animal 
   Pos2 = [];        % (same as above but row vector)
   xbins = [];       % position ("x") vector in 0,200,400 space (Pos mark space) 
   num_xbins = nan;   % total # of bins in Pos mark space
      C_length = nan; % length of Center arm
      L_length = nan; % length of Left arm
      R_length = nan; % length of Right arm

      % prior %
   T = [];           % Transition matrix (Markovian)

      % likelihood (encoding) %
   encodeperiods = [];              % [starttime endtime] of enc periods
   mark0 = cell(1,maxtet);          % {tet} Amp marks for enc spk
   Xnum = cell(1,maxtet);           % {tet} Pos Gauss kern for each enc spk
   Lint = cell(1,maxtet);           % {tet} Firing occnorm Pos map for each tet
   posind_spk_e = cell(1,maxtet);   % {tet} Pos inds of each enc spk
   numspk_e = nan(maxtet,1);        % {tet} # of enc spk on this tet
   Xnum_all = [];                   % (all enc spk) [ xbins x spk_e ] Gauss kern around the Pos of the rat during each enc spk (each column)
   Lint_all = [];                   % (all enc spk) [ xbins x 1 ] Adds spikes (each one probability 1) and divides by occupancy and pos time step
   
   
   %%% Preliminary A: Set-up Positional space (W-track linearization) %%%%%%%%%%%%%%%%%%%%%%%%      
       
   disp(['Encoding : Day ',num2str(d), ', Ep ',num2str(ep)])

   % Initialize Pos %%%%%%%%%% Pos: (W-track position linearized + concatenated (C, L, R)) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Pos  = nan(size(lindist)) ;
   Pos2 = nan(size(lindist))';   
        Carm = (lindist <= CP)         ;       % Before CP
        Larm = ~Carm & (seg2 | seg3)   ;       % L arm      
        Rarm = ~Carm & (seg4 | seg5)   ;       % R arm  

   % EXCUR - based encoding periods
   if EXTYPE == 1                   % LR Pro
       Cind = Carm & outvec;
       Lind = Larm & outvec;
       Rind = Rarm & outvec;
   elseif EXTYPE == 2               % LR Ret
       Cind = Carm & invec;
       Lind = Larm & invec;
       Rind = Rarm & invec;
   elseif EXTYPE == 3               % L inbound
       Cind = Carm & invec;
       Lind = Larm & invec;
       Rind = Rarm & outvec;
   elseif EXTYPE == 4               % R inbound
       Cind = Carm & invec;
       Lind = Larm & outvec;
       Rind = Rarm & invec;
   end

   CLRvec = ( Cind | Lind | Rind ) ;  % Encoding Pos bins (time bins) across the 3 arms

   % Linearization concatenation to padded form (0,200,400: C,L,R)
   Pos(Cind)           = lindist(Cind)          +   0;     % Center 
   Pos(Lind)           = lindist(Lind)  - CP    +   200;   % Left
   Pos(Rind)           = lindist(Rind)  - CP    +   400;   % Right
        maxL = ceil(max(Pos(Lind)));   % L arm maximum 1-cm bin in 0,200,400 space
        maxR = ceil(max(Pos(Rind)));   % R arm maximum 1-cm bin in 0,200,400 space
   Pos(isnan(Pos)) = nan;  % just a reminder that NaN means that spikes in this Pos bin was not used to encode!
   Pos2 = Pos(:)';

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
   %%%% Construct Pos space  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   x_edges = 0:xdel:maxR ;   % edge vector  (used if doing empirical transition matrix)
        xbins = 1:maxR;
        num_xbins = length(xbins);
        % indices of each Pos "island" of positions that were actually occupied (not the padding)
        aa = lookup(1,xbins);    bb = lookup(maxC,xbins);  % Center
        cc = lookup(200,xbins);  dd = lookup(maxL,xbins);  % Left
        ee = lookup(400,xbins);  ff = lookup(maxR,xbins);  % Right
       C_length = length(aa:bb) ;  % arm length in cm (1-cm bins)
       L_length = length(cc:dd) ;
       R_length = length(ee:ff) ;    
            Numlinbins = C_length + L_length + R_length;
   % Indicator vector of Pos ('x') bins that are actually occupied ('valid')        
   validxbins = zeros(size(xbins));
        validxbins(aa:bb) = 1 ;   % Center
        validxbins(cc:dd) = 1 ;   % Left
        validxbins(ee:ff) = 1 ;   % Right
        validxbins = logical(validxbins);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
            
 
       
   %%%%%% Preliminary D: collect Encoding spikes' Marks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % Variables 
   mark0           ;    % {tet} Marks (4 spike amps) for enc spk
   numspk_e        ;    % {tet} Total # of enc spk on each tet
   posind_spk_e    ;    % {tet} Pos inds of enc spk

   % Iterate tetrodes
   for tet = selected_tets

       filedata = loadparamsfile(daydir,d,tet);  % spike data from day directory
       SPK = filedata.params;

       disp(filedata.paramnames{2})
       disp(filedata.paramnames{3})
       disp(filedata.paramnames{4})
       disp(filedata.paramnames{5})

       % Set-up encoding periods: moving + arm/dir selection 
       encodeperiods = vec2list(movevec & CLRvec,postimevec);  
       %encodeperiods = vec2list(linmovevec & CLRvec,postimevec);  

       % Filter for encoding spikes
       sinds_amp     = any( SPK(:,2:5) > spikethresh , 2) ;                  % spikes >thresh uV in at least one channel
       sinds_exc     = logical( isExcluded( SPK(:,1)/10000 ,encodeperiods));  % spikes occurring during excursion-specific encoding periods
           % (optional) Filter out FS spikes ("type") to get pyramidal neuron spikes
            if exclude_FS > 0
                   maxwidths = filedata.params(:,6)/10;
                   sinds_nonFS = maxwidths > exclude_FS;
               sinds_type = sinds_nonFS;    % filter
           else
               sinds_type = true(size(SPK(:,1)));
           end
       sinds_e      = sinds_type & sinds_amp & sinds_exc;                   % final set of encoding spikes                

       % Get Mark values for Enc spikes
       disp(sprintf('%d encoding spk',sum(sinds_e)))
       spk_e{tet}           = SPK(sinds_e,1)/10000;                   % Enc spk times
            numspk_e(tet)   = length(spk_e{tet});                     % Total # of enc spk
       mark0{tet}           = SPK(sinds_e,2:5);                       % Enc spk amp mark  (4 channel amps) 
       posind_spk_e{tet}    = lookup(spk_e{tet},postimevec);          % Enc spk pos indices

       clear filedata

   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   %%% Preliminary E: Collect Decoding spk (all spk) marks %%%%%%
        %%% also, construct Amplitude mark space %%%%%%%%%%%%%%%%%%%   
          
   % Initialize
   spktimes_d       = cell(1,maxtet);   % {tet} dec spk times
   posind_spike_d   = cell(1,maxtet);   % {tet} dec spk Pos inds
   markAll          = cell(1,maxtet);   % {tet} dec spk 5D mark (Pos + 4 amp)
   
   for tet2 = selected_tets
       
       filedata2 = loadparamsfile(daydir,d,tet2);  % spike data from day directory       
       SPK2 = filedata2.params;         % matrix of spike parameters
       
       % Obtain dec spk info
           sinds_ep  = ( SPK2(:,1)/10000  >=  epstart )  &  ( SPK2(:,1)/10000 <= epend );    % spikes in epoch
           sinds_amp = any( SPK2(:,2:5) > spikethresh, 2) ;                                    % spikes of sufficient amp
           if exclude_FS > 0   
                   maxwidths = SPK2(:,6)/10;
                   sinds_nonFS = maxwidths > exclude_FS;
               sinds_type = sinds_nonFS;    % filter
           else
               sinds_type = true(size(SPK2(:,1)));
           end             
       spktimes_d{tet2}         =  SPK2(sinds_ep & sinds_amp & sinds_type,1)/10000    ;     % dec spk times  
       posind_spike_d{tet2}     =  lookup(spktimes_d{tet2},postimevec)   ;     % dec pos inds

       % Covariate - Mark of dec spk: [ Pos <amp1 amp2 amp3 amp4> ]
       markAll{tet2}(:,1)       = posind_spike_d{tet2} ;                % Col 1:    Pos inds
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
   if isempty(SPK_d)
       disp('Chapati?')
       continue
   end
   [SPK_d_sort , sinds_d] = sort(SPK_d(:,1));      % Sort in time   
       spk_d_sort_parenttet  = SPK_d(sinds_d,2);   % Parent tets
       spk_d_sort_tetspknum  = SPK_d(sinds_d,3);   % Spike # for parent tet
       
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
 %%%%% Preliminary B: construct Uniform prior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  unifprior = zeros(num_xbins,1);           % initialize
    unifprior(validxbins) = 1;              % fill only where 0,200,400 matrix is occupied
    unifprior = unifprior / sum(unifprior); % normalize
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      
  %%%% Preliminary C: construct Transition matrix (T) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  T = nan(num_xbins,num_xbins);

  if Priormodel == 0               % Uniform
      
      unidensity = max(unifprior);
      T = unidensity * ones(num_xbins,num_xbins);
      T(:,~validxbins) = 0;
      T(~validxbins,:) = 0;
      if ~all(unifprior == T(:,1)) % Sanity check
          keyboard
      end
      
  elseif Priormodel == 1            % 2D Gaussian (50 cm/s-based sigma) w/ Compression factor 10
      
      % (i)   Determine sigma of the 2D Gaussian
        dt;               % 0.020 s (20 ms)        
        vel = 50 * 5;         % 50 cm /s   -- 5 is initial compression factor
        sig = (nummsoverlap/1000) * vel;     % 1 cm 
      
      % (ii)  Transition matrix  (behavioral time-scale)
      T = Randwalk3_pad(C_length,L_length,R_length,sig);      
      
      % (iii) Compression factor (behavioral >> theta sequence)
      T = T^10;         % (Additional) compression factor  (Skaggs 1996, etc.)
                        
      % % % 50x speedup % % % % % % % % % % % % % % % % % % %
      
  elseif Priormodel == 2            % Empirical + Compression Factor
      
      disp('Empirical priormodel >> interpolating Pos to decoding bin size')
      postimevec_old = postimevec;
      postimevec_new = postimevec(1):(nummsbin/1000):postimevec(end);
      % interpolate data columns
      Posnew = interp1(postimevec(:),Pos(:),postimevec_new(:),'nearest');
      if 0
          % (optional) sanity check for interpolation procedure
          figure;
          plot(postimevec,Pos,'ro','markersize',10); hold on
          plot(postimevec_new,Posnew,'.','color','k')
          keyboard
      end
      
      posinds_enc = logical(list2vec(encodeperiods,postimevec_new));
      [ ~ , ind_x ]  =  histc(Posnew(posinds_enc),xbins);
                  % [ <current bin #>  <next bin #>  ]
      x_state  =  [ ind_x(1:end-1)   ind_x(2:end) ];
      % Iterate through each linear position (x)
      for x = 1:num_xbins
          next_x = x_state( x_state(:,1) == x , 2  ); % Collect next x from the animals data
          if ~isempty(next_x)
              T(:,x) = histc( next_x , 1:num_xbins ) / length(next_x)  ; % histogram and normalize the dist
          elseif isempty(next_x)
              T(:,x) = zeros(1,num_xbins);
          end
      end
      if 1
          T(:,~validxbins) = 0;
          for ccc = 1:num_xbins
              T(:,ccc) = T(:,ccc) + 1e-10;
              T(:,ccc) = T(:,ccc)/sum(T(:,ccc)); 
          end
          T(isnan(T)) = 0;
      end
      T = T^10;  % exponentiate by compression bin size
      
  else
      keyboard
  end
  
  if 0
      figure;
      imagesc(T); colormap gray; colormap(flipud(colormap));
      keyboard
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
   
   %%% I. (Likelihood) Construct Pos Enc model w/ Gaussian KDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   disp('Third, construct Pos Encoding model')
   
   % Initialize
   Xnum;                % {tet} Pos Gauss kern for each enc spk
   Lint;                % {tet} Firing occnorm Pos map for each tet
   
   % Identify the rat's actual Pos (during enc periods)
   posinds_enc = [];
   numposbin_enc = nan;     
   posinds_enc = logical(isExcluded(postimevec,encodeperiods)); 
        numposbin_enc = sum(posinds_enc);
   
   % i. Calculate the occupancy map (rat's actual Pos)
   X1        = repmat(xbins',1,numposbin_enc);               % [ xbins x postimebins_enc ] 0,200,400 Pos vecs repeated as columns across enc pos time bins
   MU1       = ones(length(xbins),1)  *  Pos2(posinds_enc);  % [ xbins x postimebins_enc ] Actual Pos of rat (for each enc time bin) in each column (every xbin is same value)
   normpdf1  = normpdf(X1,MU1,xsig);                        % [ xbins x postimebins_enc ] Gauss kern around the positions of the animal (each column)
   occ       = sum(normpdf1,2);                              % [ xbins x 1         ] Rat's enc period Pos occupancy from Gauss kerns
        if EXCISE
           occ(~validxbins) = 0; 
        end
        occ  = occ / sum(occ);                    
   
   % (tetrode - specific encoding)
   for tet3 = selected_tets
       
       % ii.  Xnum: Pos Gauss kern for each enc spk
           % (preliminary) collect each Enc spk's Pos
           X2  = repmat(xbins',1,numspk_e(tet3));                       % [ xbins x spk_e ]  0,200,400 Pos vecs, repeated for each enc spk (each column)
           MU2 = ones(length(xbins),1)  *  Pos2(posind_spk_e{tet3});    % [ xbins x spk_e ]  Actual Pos of rat during each enc spk, repeated over (each column)
       Xnum{tet3} = normpdf(X2,MU2,xsig);                               % [ xbins x spk_e ]  Gauss kern around the Pos of each enc spk (each column)
       if EXCISE 
          Xnum{tet3}(~validxbins,:) = 0; 
          for kk = 1:numspk_e(tet3)
              Xnum{tet3}(:,kk) = Xnum{tet3}(:,kk) / sum( Xnum{tet3}(:,kk) ); 
          end
       end
       
       % iii.  Lint: Firing occ-norm Pos map for this tet
       totalencdur = sum(encodeperiods(:,2) - encodeperiods(:,1));
       
       Lint{tet3} = sum(Xnum{tet3},2) ./ occ ./ totalencdur;            % [ xbins x 1 ]   Adds spikes (each one probability 1) and divides by occ (enc periods) and pos time step          
       if EXCISE
           Lint{tet3}(~isfinite(Lint{tet3})) = 0;
       end       
       
   end

   clear normpdf1 X1 MU1 X2 MU2
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %%% II. Decode w/ Bayes' Rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    % i. (preliminary) Set-up decoding time bins %%%%
    Tvec = [];          % Centers of dec time bins (zeroed to ep start time)
    Wins = [];          % Windows of dec time bins (zeroed to ep start time)
    if nummsoverlap == 0       % Non-overlapping dec windows
        Tvec = ( epstart:(nummsbin/1000):epend ) - epstart;
        Wins = [Tvec(1:(end-1))' Tvec(2:end)'] ;
        Tvec = mean(Wins,2);  % Re-define time vector to bin centers
    elseif nummsoverlap > 0    % Overlapping dec windows
        Tvec  = ( (epstart + nummsbin/2000) : (nummsoverlap/1000) : (epend + nummsbin/2000) ) - epstart;     
             numcenterbins = length(Tvec);
        Wins =   nan(numcenterbins,2) ;
        for gg = 1:size(Wins,1)
            Wins(gg,:) = Tvec(gg) + [-nummsbin/2000  +nummsbin/2000] ;
        end
        Tvec = mean(Wins,2);  % Re-define time vector to bin centers
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
        if EXTYPE == 1       % 1: Pro LR
            exflags = [1 3 -11 23 32];              % Outbound trajs (1,3) + Home well trackbacks (-11) + Out to Outs (23, 32)
        elseif EXTYPE == 2   % 2: Ret LR
            exflags = [2 4 -11 23 32];              % Inbound trajs (2,4) + Home well trackbacks (-11) + Out to Outs (23, 32)
        elseif EXTYPE == 3   % 3: In L
            exflags = [-33 4 32];                   % Inbound traj from L (4), L to R (32), Left well trackbacks (-33)    
        elseif EXTYPE == 4   % 4: In R
            exflags = [-22 2 23];                   % Inbound traj from R (2), R to L (23), Right well trackbacks (-22)            
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

    Tvecs = cell(Numexcur,1);
    Posts = cell(Numexcur,1);  
    
    for ex = 1:Numexcur
        
        % Excursion info
        exstart     = Excurlist(ex,1);
        exend       = Excurlist(ex,2);
        extype      = Excurlist(ex,3);
            a_dec = lookup(exstart - epstart,Tvec);  % First dec bin
            b_dec = lookup(exend - epstart,Tvec);    % Last dec bin
                numdecbin = b_dec - a_dec + 1;             % # of decoding bins in this excursion
        
        % Initialize outputs %
        Tvecs{ex} = Tvec(a_dec:b_dec);
        Posts{ex} = zeros(num_xbins,numdecbin) ;              % Matrix of posteriors                 
        
        % (if Prior not uniform) Initialize excursion posterior with uniform prior
        if Priormodel > 0
            Post = unifprior;
        end
        
        tic
        % Iterate each dec time bin
        for bbb = 1:numdecbin
            
            t = bbb - 1 + a_dec;                % Tvec index
            decbin = Wins(t,:) + epstart;   % [starttime endtime]
            midtime = mean(Wins(t,:),2);    % [ midtime ]

            %%% Prior %%%%
            if Priormodel == 0      % uniform
                prior = unifprior;  
            elseif Priormodel > 0  % Empirical
                prior = T * Post;
            end
            
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
                        l1 = Xnum{tet} * l0' ./ occ ./ totalencdur ;      % small lambda      ( Xnum is enc spk Pos Gaus kern )
                        l2 = ( l1 * dt ) .* exp( -Lint{tet} .* dt ) ;     % full likelihood
                        l2 = l2 + epsilon;                                % add infinitesimal value (prevents normalization error)
                        l2 = l2 ./ sum(l2(isfinite(l2))) ;                % normalize
                        l(:,s) = l2 ;
                        if any(isnan(l2(validxbins)))
                            keyboard
                        end
                        
                    end
                    
                    % Tetrode likelihood: product of single-spike likelihoods
                    if ~EXCISE
                        L_tets(:,ttt) = prod(l,2) ;                             
                    else
                        L_tets(validxbins,ttt) = prod(l(validxbins,:),2) ;                            
                    end

                    if sum(L_tets(:,ttt)) == 0
                        % FIX HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        keyboard
                    end
                    
                    L_tets(:,ttt) = L_tets(:,ttt) / sum(L_tets(:,ttt)) ;      % normalize
                    
                    disp(sprintf('(Pos %d) %s d%d ep%d: %s',EXTYPE,animalname(1:3),d,ep,num2str(midtime)))
                    
                end
                
            end
                
            L = prod(L_tets,2);    % Full likelihood: product of tetrode likelihoods
            
            if any(isnan(L))
                keyboard
            end
            
            Post = prior .*  L ;         % Bayes' rule
            %Post = L;
                Post = Post / sum(Post) ;  % normalize
            
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
                    imagesc(Tvecs{ex}(a:b),1:Numlinbins,Posts{ex}(validxbins,a:b),[0 .05]); hold on
                    % arm guide lines
                    plot([plot_start plot_end],[C_length C_length]+0.5,'k--');
                    plot([plot_start plot_end],[C_length C_length] + L_length +0.5,'k--');
                    %colormap hot;
                    colormap gray ; colormap(flipud(colormap));
                    set(gca,'ydir','normal')
                    plot(postimevec(aa:bb) - epstart,Pos(aa:bb),'color',[0 .8 0],'linewidth',2);
                    str1 = sprintf('%s %d %d',animalname(1:3),d,ep);
                    str2 = sprintf('Ha: %d uv, Hx: %d cm, nummsbin: %d, nummsoverlap: %d',msig,xsig,nummsbin,nummsoverlap);
                    title({str1,str2},'fontsize',14,'fontweight','bold')
                    keyboard
                end
            end

        end
        
        % Excise padding Pos bins before saving
        Posts{ex} = Posts{ex}(validxbins,:);
       
    end
    
    % Save output in a file
    disp('Done calculating epoch')
    clear P
    P.date              = date;
    P.animalname        = animalname;
    P.dayep             = [d ep];
    P.EXTYPE            = EXTYPE;
    P.EXCISE            = EXCISE;
    P.divider           = '%%%%%%%%%%%%%%%%';    
    P.spikethresh       = spikethresh;
    P.exclude_FS        = exclude_FS;
    P.selected_tets     = selected_tets;
    P.msig              = msig;
    P.xdel              = xdel;
    P.xsig              = xsig;
    P.Priormodel        = Priormodel;
    P.Excurbuffer       = Excurbuffer;
    P.nummsbin          = nummsbin;
    P.nummsoverlap      = nummsoverlap;
    P.divider1           = '%%%%%%%%%%%%%%%%';
    P.postimevec = postimevec;
    P.posvec = Pos;          % Actual animal position state (only encoding times -- otherwise nans)
    P.yvec_pos = 1:Numlinbins;             
    P.CLR_lengths        = [C_length L_length R_length];
    P.divider2           = '%%%%%%%%%%%%%%%%';   
    P.excurlist          = Excurlist;
    P.excurlist_descript  = excurlist_descript;
    P.divider3           = '%%%%%%%%%%%%%%%%';    
    P.Priormodel        = Priormodel;   % Flag for which prior used
    P.Posts             = Posts;        % Posteriors {excur}
    P.Tvecs             = Tvecs;        % Time vectors {excur}
    cd(savedir);
    savefilename = sprintf('%s_%d_%d_Pos_%d',animalname(1:3),d,ep,EXTYPE);
    save(savefilename,'P','-v7.3')
    clear P Posts Tvecs
    
end

end








