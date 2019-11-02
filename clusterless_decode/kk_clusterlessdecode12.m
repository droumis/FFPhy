function kk_clusterlessdecode12(directoryname,fileprefix,dayeps,animalname,varargin)
  
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
        savefilename = sprintf('%s_%d_%d_Pos_%d',animalname(1:3),d,ep,EXTYPE);
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
   [dummy,~] = evaluatetimefilter(animalinfo{2},animalinfo{3},{ linvel4 },[d ep]);
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
            
  %%%%% Preliminary B: construct Uniform prior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  unifprior = zeros(num_xbins,1);           % initialize
    unifprior(validxbins) = 1;              % fill only where 0,200,400 matrix is occupied
    unifprior = unifprior / sum(unifprior); % normalize
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      
  %%%% Preliminary C: construct Transition matrix (T) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  T = nan(num_xbins,num_xbins);

  if Priormodel == 0   % uniform
      unidensity = max(unifprior);
      T = unidensity * ones(num_xbins,num_xbins);
      T(:,~validxbins) = 0;
      T(~validxbins,:) = 0;
      if ~all(unifprior == T(:,1)) % Sanity check
          keyboard
      end
  else
      keyboard
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
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
            if exclude_FS    
                   maxwidths = filedata.params(:,6)/10;
                   sinds_nonFS = maxwidths > .2;
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
           if exclude_FS    
                   maxwidths = SPK2(:,6)/10;
                   sinds_nonFS = maxwidths > .2;
               sinds_type = sinds_nonFS;    % filter
           else
               sinds_type = true(size(SPK2(:,1)));
           end             
       spktimes_d{tet2}         =  SPK2(sinds_ep & sinds_amp & sinds_type,1)/10000    ;     % dec spk times  
       posind_spike_d{tet2}     =  lookup(spktimes_d{tet2},postimevec)   ;     % dec pos inds

       % Covariate - Mark of dec spk: [ Pos <amp1 amp2 amp3 amp4> ]
       markAll{tet2}(:,1)       = posind_spike_d{tet2} ;                % Col 1:    Pos inds
       markAll{tet2}(:,2:5)     = SPK2(sinds_ep & sinds_amp,2:5);       % Col 2-5:  amps
       
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
   [SPK_d_sort , sinds_d] = sort(SPK_d(:,1));      % Sort in time   
       spk_d_sort_parenttet  = SPK_d(sinds_d,2);   % Parent tets
       spk_d_sort_tetspknum  = SPK_d(sinds_d,3);   % Spike # for parent tet
       
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
   
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
        occ  = occ / sum(occ);                    
   
   % (tetrode - specific encoding)
   for tet3 = selected_tets
       
       % ii.  Xnum: Pos Gauss kern for each enc spk
           % (preliminary) collect each Enc spk's Pos
           X2  = repmat(xbins',1,numspk_e(tet3));                       % [ xbins x spk_e ]  0,200,400 Pos vecs, repeated for each enc spk (each column)
           MU2 = ones(length(xbins),1)  *  Pos2(posind_spk_e{tet3});    % [ xbins x spk_e ]  Actual Pos of rat during each enc spk, repeated over (each column)
       Xnum{tet3} = normpdf(X2,MU2,xsig);                              % [ xbins x spk_e ]  Gauss kern around the Pos of each enc spk (each column)
       
       % iii.  Lint: Firing occ-norm Pos map for this tet
       totalencdur = sum(encodeperiods(:,2) - encodeperiods(:,1));
       Lint{tet3} = sum(Xnum{tet3},2) ./ occ(:,1) ./ totalencdur;            % [ xbins x 1 ]   Adds spikes (each one probability 1) and divides by occ (enc periods) and pos time step
       
       % ******* Should we set non-occupied Pos (xbins) indices to 0? *******
   end

   clear normpdf1 X1 MU1 X2 MU2
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %%% II. (Likelihood) Calculation of no-spike likelihood  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('Fifth, calculate no-spike likelihood')    

    Xnum_all = [];      % [ xbins x spk_e ]     Gauss kern around the Pos of the rat for each enc spk (each column)
    Lint_all = [];      % [ xbins x 1 ]         Adds spikes (each one probability 1) and divides by occ (enc periods) and pos time step
    occ;
        
    % Xnum: Gauss kernel for Pos, from enc spks
    num_allspk_e = 0;        % Total # of enc spks
    posind_allspk_e = [];    % Pos inds of enc spks
    for vv = 1:maxtet   % Collect Pos of every enc spk -over all tets-
        if numspk_e(vv) > 0
            num_allspk_e = num_allspk_e + numspk_e(vv);
            posind_allspk_e = [posind_allspk_e ; posind_spk_e{vv}(:)];
        end
    end
    Xnum_all = normpdf( xbins'*ones(1,num_allspk_e),...     % [ xbins x spk_e ] Gauss kern around the Pos of the rat during each enc spk (each column)
                        repmat(Pos2(posind_allspk_e),num_xbins,1),...   
                        xsig );

    % Set non-occupied xbins indices to 0
    %Xnum_all(~validxbins,:) = 0;
    % normalize for each spike
    %Xnum_all = bsxfun(@rdivide,Xnum_all,sum(Xnum_all,1));
    %for skk = 1:size(Xnum_all,2)
    %    Xnum_all(:,skk) = Xnum_all(:,skk) / sum(Xnum_all(:,skk));
    %end
    
    % Lint: conditional intensity function for the no-spike case
    Lint_all = sum(Xnum_all,2)  ./  occ(:,1) ./  totalencdur ;       % [ xbins x 1 ] Adds spikes (each one probability 1) and divides by occupancy and pos time step
    %Lint_all = Lint_all ./ sum(Lint_all)  ;                     % Normalize
        clear Xnum_all
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% III. Decode w/ Bayes' Rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    % i. (preliminary) Set-up decoding time bins %%%%
    Tvec = [];          % Centers of dec time bins (zeroed to ep start time)
    Wins = [];          % Windows of dec time bins (zeroed to ep start time)
    if nummsoverlap == 0       % Non-overlapping dec windows
        Tvec = (epstart:(nummsbin/1000):epend) - epstart;
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
        Excurlist           = excursions{d}{ep}.excurlist;
        excurlist_descript  = excursions{d}{ep}.excurlist_cols;
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
        Excurlist   = Excurlist(ismember(Excurlist(:,3),exflags),:);  % select subset of excursions
        % Add extra time buffer before and after excursion
        Excurlist(:,1) = Excurlist(:,1) - Excurbuffer;
        Excurlist(:,2) = Excurlist(:,2) + Excurbuffer;
        Numexcur    = size(Excurlist,1);
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
        
        %%% Initialize posterior with uniform prior
        %postx = unifprior;
        
        tic
        % Iterate each dec time bin
        for bbb = 1:numdecbin
            
            t = bbb - 1 + a_dec;                % Tvec index
                decbin = Wins(t,:) + epstart;   % [starttime endtime]     
                midtime = mean(Wins(t,:),2);    % [ midtime ]
                
            %%% Prior %%%%
            prior = unifprior;  
            
            %%% Likelihood %%%
            
            L = [];  % initialize
            
            % Spikes that occurred w/in this dec bin    
            sinds = find( isExcluded( SPK_d_sort, decbin ) ) ;         
                numspk = length(sinds);
            
            if numspk == 0
                
                L = exp(-Lint_all .* dt) ;   % % Eq'n 2.5 (Deng et. al. 2015)
                
            elseif numspk > 0
                
                % Initialize matrix of single-spike likelihoods
                l = zeros(num_xbins,numspk) ;   % [ xbins x numspk_decbin ] 
                
                % Iterate each dec spk
                for s = 1:numspk
                    
                    % Spike info
                    ind = sinds(s);                    % Sorted spike index
                    tet = spk_d_sort_parenttet(ind) ;  % Parent tet
                    i = spk_d_sort_tetspknum(ind) ;    % Spike # on parent tet
                    
                    % l0: Amp mark likelihoods for this spike
                    amp      = markAll{tet}(i,2:5);
                    amp_enc  = mark0{tet};
                    enc_1vec = ones(1,numspk_e(tet));
                        
                    l0 = normpdf( amp(1) * enc_1vec,...   % x:      amp  *   [ 1 x numspk_enc ]    % dec spk amplitude
                                  amp_enc(:,1)',...       % mu:     amp_e    [ 1 x numspk_enc ]    % enc spk amplitudes
                                  msig ) .* ...           % smker:  amp bandwidth                  % Gaus kern SD
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
                    l1 = Xnum{tet} * l0' ./ occ(:,1) ./ totalencdur ;      % small lambda      ( Xnum is enc spk Pos Gaus kern )
                    l2 = ( l1 * dt ) .* exp( -Lint{tet} .* dt ) ;           % full likelihood
                        l2 = l2 + epsilon;                              % add infinitesimal value (prevents normalization error)
                            l2 = l2 ./ sum(l2(isfinite(l2))) ;           % normalize
                    l(:,s) = l2 ;  
                    
                end
                
                L = prod(l,2) ;     % independence: product of single-spike likelihoods
                L = L/sum(L) ;     % normalize
                disp(sprintf('(Pos %d) %s d%d ep%d: %s',EXTYPE,animalname(1:3),d,ep,num2str(midtime)))
            
            end
            
            postx = prior .*  L ;
            %postx = L;
            postx = postx / sum(postx) ;  % normalize
            
            Posts{ex}(:,bbb) = postx;
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
                    imagesc(Tvecs{ex}(a:b),1:Numlinbins,Posts{ex}(validxbins,a:b),[0 .04]); hold on
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
    cd(savedir); disp('Done calculating epoch')
    clear P
    P.date              = date;
    P.animalname        = animalname;
    P.dayep             = [d ep];
    P.EXTYPE            = EXTYPE;
    P.divider           = '%%%%%%%%%%%%%%%%';    
    P.spikethresh       = spikethresh;
    P.nummsbin          = nummsbin;
    P.nummsoverlap      = nummsoverlap;
    P.selected_tets     = selected_tets;
    P.mdel              = mdel;
    P.smker             = msig;
    P.xdel              = xdel;
    P.sxker             = xsig;
    P.divider1           = '%%%%%%%%%%%%%%%%';
    P.postimevec = postimevec;
    P.posvec = Pos;                      % Actual animal position state
    P.yvec_pos = 1:Numlinbins;             
    P.CLR_lengths   = [C_length L_length R_length];
    P.divider2           = '%%%%%%%%%%%%%%%%';   
    P.excurlist          = Excurlist;
    P.excurlist_descript  = excurlist_descript;
    P.divider3           = '%%%%%%%%%%%%%%%%';    
    P.Priormodel        = Priormodel;   % Flag for which prior used
    P.Posts             = Posts;        % Posteriors {excur}
    P.Tvecs             = Tvecs;        % Time vectors {excur}
    savefilename = sprintf('%s_%d_%d_Pos_%d',animalname(1:3),d,ep,EXTYPE);
    save(savefilename,'P','-v7.3')
    clear P Posts Tvecs
    
end

end








