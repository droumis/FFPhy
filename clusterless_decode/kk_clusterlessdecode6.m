function kk_clusterlessdecode6(directoryname,fileprefix,dayeps,animalname,varargin)

% dirpos 

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
sigma_transmat = 1;
plot_infunction = 0;
TETSET = 5;
nummsbin = 1; 5;
dec_start = 0; 72;
epsi = eps(0);
SKIPSAVED = 0;
CELLMAX = 0;

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
        case 'sigma_transmat'
            sigma_transmat = varargin{option+1};            
        case 'TETSET'
            TETSET = varargin{option+1};
        case 'spikethresh'
            spikethresh = varargin{option+1};
             
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
for de = 1:size(dayeps,1)
    
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
          selected_tets = adtc_nonN(animalind & dayind,3);
            selected_tets = selected_tets(:)';
          maxtet = max(selected_tets);
    end           
        
   % (if specified) Identify the remote W-track epoch # 
       epochs_enc = ep;
       remoteep = [];
   num_encodeeps = length(epochs_enc);   % [local epoch #, remote epoch #]
   
   % Initialize outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
      % positional variables
   xbins = cell(1,2);
   armdists = cell(1,2);        % two elements since up to two different W-tracks
   validxbins = cell(1,2);   
   centerarmmax = [nan nan];
   rightarmmax = [nan nan];
      % prior variables %
   stateM_gausnorm = cell(1,2);
      % likelihood variables %
   mark0 = cell(1,2);
   Xnum = cell(1,2);
   Lint = cell(1,2);
   numspikes_e = cell(1,2);
   posind_spike_e = cell(1,2);
      % internal use (not to save in output file later)
   num_xbins = [nan nan];
   armdists2 = cell(1,2);
   armdists_cat = cell(1,2);
      a_cut = cell(1,2);
      b_cut = cell(1,2);
      c_cut = cell(1,2);
      d_cut = cell(1,2);
      e_cut = cell(1,2);
      f_cut = cell(1,2);
  encodeperiods = cell(1,2);    

   %%% First, iterate through each encoding epoch to collect basic data %%%%%%%%%%%%%%%%%%%%%%%%      
   for nn = 1:num_encodeeps
       
       epenc = epochs_enc(nn);
       
       disp(['Encoding : Day ',num2str(d), ', Epoch ',num2str(epenc)])
       
       % Basic epoch data
       postimevec_epenc = linpos{d}{epenc}.statematrix.time;
            numpossamps_encep = length(postimevec_epenc);
            starttime_epenc = postimevec_epenc(1);
            endtime_epenc = postimevec_epenc(end);
       
       % Identify W-track segment indices
       lindist = linpos{d}{epenc}.statematrix.lindist;
       seg1 = linpos{d}{epenc}.statematrix.segmentIndex==1;
       seg2 = linpos{d}{epenc}.statematrix.segmentIndex==2;
       seg3 = linpos{d}{epenc}.statematrix.segmentIndex==3;
       seg4 = linpos{d}{epenc}.statematrix.segmentIndex==4;
       seg5 = linpos{d}{epenc}.statematrix.segmentIndex==5;
       CPbuff = choicepointer(linpos{d}{epenc});
       % Headdir
       headdir = linpos{d}{epenc}.statematrix.segmentHeadDirection(:,1);   % head direction relative to center well -- values > 0 are outbound
            postimevec_nonnan = postimevec_epenc(~isnan(headdir));
            headdir_nonnan = headdir(~isnan(headdir));
            headdir2 = interp1(postimevec_nonnan,headdir_nonnan,postimevec_epenc,'linear');          

       % Initialize dirpos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       dirpos{nn} = nan(size(lindist)) ;   % (horizontal vector)  0+ (OUT) +500 if R, 1000+ (IN), +500 if R  lin dist from center well + L vs. R
       dirpos2{nn} = nan(size(lindist));
            centerarmmax(nn) =  CPbuff;   % maximum linear distance of center arm 
                CAinds = lindist < CPbuff;
                Rinds = ~CAinds & (seg4 | seg5);
                Linds = ~CAinds & (seg2 | seg3);
                times1 = CAinds & (headdir2 >= 0);   % t1: center arm,  headed TO center junct
                times2 = CAinds & (headdir2 <  0);   % t2: center arm,  headed AWAY FROM center junct
                times3 = Rinds & (headdir2 >= 0);    % t3: right arm, headed AWAY FROM center junct
                times4 = Rinds & (headdir2 <  0);    % t4: right arm, headed TO center junct   
                times5 = Linds & (headdir2 >= 0);    % t5: left arm, headed AWAY FROM center junct   
                times6 = Linds & (headdir2 <  0);    % t6: left arm, headed TO center junct                   

                % (sanity check) Plot headdir segments    
                if 1
                    sttime = pos{d}{epenc}.data(1,1);
                    figure;       hold on
                        plot(pos{d}{epenc}.data(times1,1)-sttime,headdir2(times1),'.','markersize',10,'Color',[1 0 0]);   % 1: red
                        plot(pos{d}{epenc}.data(times2,1)-sttime,headdir2(times2),'.','markersize',10,'Color',[0 0 1]);   % 2: blue
                        plot(pos{d}{epenc}.data(times3,1)-sttime,headdir2(times3),'.','markersize',10,'Color',[1 .9 0]);  % 3: yellow
                        plot(pos{d}{epenc}.data(times4,1)-sttime,headdir2(times4),'.','markersize',10,'Color',[0 0 0]);   % 4: black
                        plot(pos{d}{epenc}.data(times5,1)-sttime,headdir2(times5),'.','markersize',10,'Color',[0 .8 0]);   % 5: green
                        plot(pos{d}{epenc}.data(times6,1)-sttime,headdir2(times6),'.','markersize',10,'Color',[1 .4 1]);  % 6: magenta                      
                    keyboard
                end  
                
        % Now convert headdir to Wu-Foster-2014's "inbound" and "outbound", which are in reference to the center junction
            % now, "INBOUND" is defined as TOWARD center junction (200, 400, 600)
            % now, "OUTBOUND" is defined as AWAY FROM center junction (1200, 1400, 1600)
        dirpos{nn}( times1 )   =   lindist( times1 )                      +0;  %+   200;       % >200: center arm + toward center junction
        dirpos{nn}( times2 )   =   lindist( times2 )                      +200;  %+   1200;          % >1200: center arm + away from center junction

        dirpos{nn}( times3 )   =   lindist(times3)    - centerarmmax(nn)  +1000   ;%1400;  % >1400: right arm + away from center junction
        dirpos{nn}( times4 )   =   lindist(times4)    - centerarmmax(nn)  +800   ;%400;  % >400: right arm + toward center junction

        dirpos{nn}( times5 )   =   lindist(times5)    - centerarmmax(nn)  +600   ;%1600;  % >1600: left arm + away from center junction
        dirpos{nn}( times6 )   =   lindist(times6)    - centerarmmax(nn)  +400   ;%600;  % >600: left arm + toward center junction

                dirpos{nn} = dirpos{nn}(:)   ;
                dirpos2{nn} = dirpos{nn}(:)' ;            
                                
                
              
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       %%%% Construct positional mark space (xbins)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       xbins{nn} = 0 : xdel : (max(dirpos{nn}) + xdel);  %  min(trajpos{nn}):xdel:(max(trajpos{nn}) + xdel);        %   
             
                num_xbins(nn) = length(xbins{nn});
             
                %indices for each positional "island" a to b, c to d, e to f
                a_cut{nn} = 1;
                b_cut{nn} = lookup(CPbuff,xbins{nn});
                
                c_cut{nn} = lookup(200,xbins{nn});
                d_cut{nn} = lookup(200 + CPbuff,xbins{nn});
                
                e_cut{nn} = lookup(400,xbins{nn});
                f_cut{nn} = lookup(400 + max(dirpos{nn}(times5 | times6)) - CPbuff + 1,xbins{nn});
                
                g_cut{nn} = lookup(600,xbins{nn});
                h_cut{nn} = lookup(600 + max(dirpos{nn}(times5 | times6)) - CPbuff + 1,xbins{nn});

                i_cut{nn} = lookup(800,xbins{nn});
                j_cut{nn} = lookup(800 + max(dirpos{nn}(times3 | times4)) - CPbuff + 1,xbins{nn});
                
                k_cut{nn} = lookup(1000,xbins{nn});
                l_cut{nn} = lookup(1000 + max(dirpos{nn}(times3 | times4)) - CPbuff + 1,xbins{nn});                
                
       validxbins{nn} = zeros(size(xbins{nn}));
       validxbins{nn}([a_cut{nn}:b_cut{nn}   c_cut{nn}:d_cut{nn}   e_cut{nn}:f_cut{nn} ...
                       g_cut{nn}:h_cut{nn}   i_cut{nn}:j_cut{nn}   k_cut{nn}:l_cut{nn}]) = 1 ;   % indicates which xbins indices were actually occupied by animal
            validxbins{nn} = logical(validxbins{nn});
                     
            xbins_cut{nn} = xbins{nn}(validxbins{nn});
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      %  Construct the transition matrix from animal's behavior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       stateM = zeros(num_xbins(nn),num_xbins(nn));
            [ ~ , ind_xbin ]  =  histc(dirpos{nn},xbins{nn});
                            % [ <current bin #>  <next bin #>  ]
                xbinstate  =  [ ind_xbin(1:end-1)   ind_xbin(2:end) ];  
       % Iterate through each linear position
       for lp = 1:num_xbins(nn)
           % nextstates : for a given bin, collection of next bin # of the animal's linear position from actual data
           next_xbin = xbinstate( xbinstate(:,1)==lp , 2  );
           if ~isempty(next_xbin)
               % histogram and normalizes the distribution
               stateM(:,lp) = histc( next_xbin , 1:num_xbins(nn) ) / length(next_xbin)  ;
           elseif isempty(next_xbin)
               stateM(:,lp) = zeros(1,num_xbins(nn));
           end
       end
       
       %    Smooth the transition matrix 
       if 0
           K = inline('exp(-(x.^2+y.^2)/2/sig^2)');                      % Gaussian
           [dx,dy] = meshgrid( [-1:1] );
           sig = sigma_transmat;
           weight = K(sig,dx,dy)/sum(sum(K(sig,dx,dy)));                 % Normalizing weights
           stateM_gaus = conv2(stateM, weight, 'same');                  % Gaussian smoothed
           stateM_gausnorm = stateM_gaus * diag(1./sum(stateM_gaus,1));  % Normalized to confine each bin (columns) probability to 1
           %stateM_gausnorm(isnan(stateM_gausnorm)) = 0;
       elseif 1
           K = inline('exp(-(x.^2+y.^2)/2/sig^2)');                      % Gaussian
           [dx,dy] = meshgrid( [-2:2] );
           sig = sigma_transmat;
           weight = K(sig,dx,dy)/sum(sum(K(sig,dx,dy)));                 % Normalizing weights
           stateM_gaus = conv2(stateM, weight, 'same');                  % Gaussian smoothed
           %stateM_gausnorm{nn} = stateM_gaus * diag(1./sum(stateM_gaus,1));  % Normalized to confine each bin (columns) probability to 1
           %stateM_gausnorm{nn}(  isnan(stateM_gausnorm{nn}) ) = 0; %0.001;     
           %stateM_gausnorm{nn} = stateM_gaus;  % don't normalize, cut off later 1.19.16 kk dl
           M = stateM_gaus;
           %stateM_gausnorm{nn}(~validxbins{nn},:) = 0;
           %stateM_gausnorm{nn}(:,~validxbins{nn}) = 0;
          
           for cc = 1:size(M,2)
               M(:,cc) = M(:,cc) / sum(M(:,cc));
           end
            M(isnan(M)) = 0;
           stateM_gausnorm{nn} = M;
           
       elseif 0
           gkernel = gaussian2(1,8);
           stateM_gaus = conv2(stateM,gkernel,'same');
           stateM_gausnorm = stateM_gaus * diag(1./sum(stateM_gaus,1));
           stateM_gausnorm(isnan(stateM_gausnorm)) = 0.0001;
       else
           stateM_gausnorm = stateM;
       end
       
       if 0
           keyboard
           figure; imagesc(stateM_gausnorm{nn})
       end
       
       %stateM_gausnorm = stateM_gausnorm + .0001;
       
       %%%%%%  Create encoding model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       % initialize epoch cell, will fill with tetrode data
       
       mark0{nn} = cell(1,maxtet);
       numspikes_e{nn} = nan(1,maxtet);
       posind_spike_e{nn} = cell(1,maxtet);
       
       for tet = selected_tets
           
           filedata = loadparamsfile(daydir,d,tet);  % spike data from day directory
           
           % Identify spikes for encoding 
           
           %  Spikes in experimenter-transcribed (in notebook) epoch
           inds_epenc =  (  filedata.params(:,1)/10000  >=  starttime_epenc  )  &  ( filedata.params(:,1)/10000 <= endtime_epenc );
           %  Spikes that has spikethresh (uV) voltage amplitude in at least one channel
           inds_thresh = any(filedata.params(:,2:5) > spikethresh, 2) ;
           %  Spikes in trajencode (non-low speed and non-well)
           timefilterscript
           [dummy,~] = evaluatetimefilter(animalinfo{2},animalinfo{3},{vel4},[d epenc]);
                vel4_periods = dummy{d}{epenc};
           inds_moving = isExcluded(filedata.params(:,1)/10000,vel4_periods)  ;
           
           % Determine choice of encoding spikes :  "model number"
           if modelnum == 2
               % exclusion case               
               % inds4 = isExcluded(filedata.params(:,1)/10000,gfperiods);
               inds_e = inds_epenc & inds_thresh & inds_moving;
                    encodeperiods{nn} = vel4_periods;
               disp(sprintf('%d spikes',sum(inds_e)))         
           else
               keyboard
           end
           
           % Spike times (recording clock) for E and D 
           spktimes_e{tet} = filedata.params(inds_e,1)/10000;    % exclusive, to ENCODE
                numspikes_e{nn}(tet) = length(spktimes_e{tet});
                
           % Amplitude mark vector (4 channel amplitudes) 
           mark0{nn}{tet} = filedata.params(inds_e,2:5);
           
           % Obtain position indices for each spike 
           posbin = (postimevec_epenc(2)-postimevec_epenc(1));
           centerpostimes = [postimevec_epenc ; postimevec_epenc(end) + posbin]  - posbin/2;
           [        ~ ,        posind_spike_e{nn}{tet} ]  =  histc(spktimes_e{tet},centerpostimes)  ;    % ENCODING
                
           clear filedata
       end
   end
   
   %%% Second, for decoding epoch, collect amplitude mark data %%%%%%%%%%%%%%%%%%%%%%%%   
   
   ep;  % the decoding epoch #
   
   % Basic epoch data
   postimevec = {};
   numpossamps = [nan nan];
   starttime = [nan nan];
   endtime = [nan nan]; 
   for zz = 1:num_encodeeps
        ep2 = epochs_enc(zz);
        postimevec{zz} = linpos{d}{ep2}.statematrix.time;
        numpossamps(zz) = length(linpos{d}{ep2}.statematrix.time);
        starttime(zz) = postimevec{zz}(1);
        endtime(zz) = postimevec{zz}(end);
   end
   clear linpos
                  
   % Collect decoding spike data
   spktimes_d = cell(1,maxtet);
   posind_spike_d = cell(1,maxtet);
   markAll = cell(1,maxtet);
   minamp = inf;
   maxamp = -inf;
   for tet2 = selected_tets
       filedata = loadparamsfile(daydir,d,tet2);  % spike data from day directory       
       inds_ep =  (  filedata.params(:,1)/10000  >=  starttime(1)  )  &  ( filedata.params(:,1)/10000 <= endtime(1) );
       inds_thresh = any(filedata.params(:,2:5) > spikethresh, 2) ;
       spktimes_d{tet2} = filedata.params(inds_ep & inds_thresh,1)/10000;
            posbin2 = (postimevec{1}(2)-postimevec{1}(1));
            centerpostimes2 = [postimevec{1} ; postimevec{1}(end)+posbin2] - posbin2/2;  % centers
       [ ~ , posind_spike_d{tet2} ]  =  histc(spktimes_d{tet2},centerpostimes2)   ;         
       markAll{tet2}(:,1) = posind_spike_d{tet2} ;         % Col 1: position time bin # , Col 2-5: spike amplitudes
       markAll{tet2}(:,2:5) = filedata.params(inds_ep & inds_thresh,2:5);
       % also, track of max and min spike amplitudes in decoding spike set
       % + in the remote-W encoding spike set
       minamp_tet_local = min(min(markAll{tet2}(:,2:5)));
       maxamp_tet_local = max(max(markAll{tet2}(:,2:5)));
       if remoteW == 1
            minamp_tet_remote_enc = min(min(mark0{2}{tet2}));
            maxamp_tet_remote_enc = max(max(mark0{2}{tet2}));
       else
            minamp_tet_remote_enc = [];
            maxamp_tet_remote_enc = [];           
       end
       minamp_tet = min( [ minamp_tet_local  minamp_tet_remote_enc  ] );
       maxamp_tet = max( [ maxamp_tet_local  maxamp_tet_remote_enc  ] );
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
        mbins = spikethresh:mdel:maxamp;
   end
   
   %%% Third, construct encoding model (Gaussian KDE) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   disp('Third, construct encoding model')
   
   Xnum{nn} = cell(1,maxtet);
   Lint{nn} = cell(1,maxtet);
   
   % need these for accurate occupancy spatial maps
   posinds_enc = cell(1,2);
   numpossamps_enc = [nan nan];  
   
   for nnn = 1:num_encodeeps

           % actual animal positions in the encoding epoch
                % but to do correctly, need to only look at encoding times
           posinds_enc{nnn} = logical(isExcluded(postimevec{nnn},encodeperiods{nnn}));
           numpossamps_enc(nnn) = sum(posinds_enc{nnn});     
           X1  = repmat(xbins{nnn}',1,numpossamps_enc(nnn));         % position bin vectors repeated as columns
           MU1 = ones(length(xbins{nnn}),1)  *  dirpos2{nnn}(posinds_enc{nnn});  % actual positions of the animal       

           % 1.  occ: occupancy (linear position) map
           normpdf1 = normpdf(X1,MU1,sxker);
           occ = repmat(sum(normpdf1,2),1,length(mbins));  %normpdf1 * ones(numpossamps,length(mbins))   ;
           occ(logical(~validxbins{nnn}),:) = 0;
           occ = occ/sum(occ(:,1));
           
           for tet3 = selected_tets
               
               % encoding spike positions
               X2  = repmat(xbins{nnn}',1,numspikes_e{nnn}(tet3));                                % position vectors repeated as columns
               MU2 = ones(length(xbins{nnn}),1)  *  dirpos2{nnn}(posind_spike_e{nnn}{tet3});    % actual (encoding) spike positions
               
               % 2.  Xnum: Gaussian kernel estimators for position, for each encoding spike
               Xnum{nnn}{tet3} = normpdf(X2,MU2,sxker);
               % normalize for each spike
               %for sk = 1:size(Xnum{nnn}{tet3},2)
               %    Xnum{nnn}{tet3}(:,sk) = Xnum{nnn}{tet3}(:,sk) / sum(Xnum{nnn}{tet3}(:,sk));
               %end
               %keyboard
               % 3.  Lint: firing occ-norm position map for tetrode
               Lint{nnn}{tet3} = sum(Xnum{nnn}{tet3},2) ./ occ(:,1) ./ dt;   %integral (?)
               % set non-occupied xbins indices to 0
               Lint{nnn}{tet3}(~validxbins{nnn}) = 0;
               Lint{nnn}{tet3} = Lint{nnn}{tet3} ./ sum(Lint{nnn}{tet3}) ;  % normalization
               
           end
   end
  
   clear normpdf1 X1 MU1 X2 MU2
   
   
   %%% Fourth, collect decoding (all) spikes in each epoch
 
   disp('Fourth, collect decoding (all) spikes in each epoch')
   
        % variables to keep track of for both local and remote epochs
   spktimes_all         = cell(1,2);       % all spike times in each epoch
   spktimes_all_sorted  = cell(1,2);       % sorted spike times, sorted
   
   for nn = 1:num_encodeeps
       
       % Collect all spikes
       for tet = selected_tets
           numspikestet_all        =  length(spktimes_d{tet});
                                                        %     [  <spike time>     <parent tet #>                  <spike # on parent tet>    <pos inds of spikes>]
           spktimes_all{nn}        =  [spktimes_all{nn}     ;   spktimes_d{tet}   tet * ones(numspikestet_all,1)    (1:numspikestet_all)'    posind_spike_d{tet} ];   
       end
       % for decoding epoch only, store this
       if nn == 1
           [spktimes_all_sorted{nn} , sinds] = sort(spktimes_all{nn}(:,1));
           num_d_spikes                      = length(spktimes_all_sorted{nn});
       end 
       
   end
     
   % Also construct a spike-tetrode matrix for the local (decoding) epoch
       % i.e. an "indicator matrix" : [  <spike #>  x  <which tetrode spikes> ]
      % initialize
  
   spike_d_parenttet = nan(num_d_spikes,1);
   spike_d_tetspikenum = nan(num_d_spikes,1);
       spike_d_parenttet = spktimes_all{1}(sinds,2);
       spike_d_tetspikenum  = spktimes_all{1}(sinds,3);
       
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
    
    Lint_all = cell(1,2);
    occ = cell(1,2);    % save this for later
    
    for n = 1:num_encodeeps
    
        Xnum_all = [];
        occ{n} = nan(1,length(xbins{n}));
        
        % occ: behavioral occupancy map in encoding epoch
        occ{n} = normpdf( xbins{n}'*ones(1,numpossamps_enc(n)),...
                       repmat(dirpos2{n}(posinds_enc{n}),num_xbins(n),1),...
                       sxker )   *    ones(numpossamps_enc(n),length(mbins) )    ;
            % remove invalid rows and normalize
            occ{n}(logical(~validxbins{n}),:) = 0;
            occ{n} = occ{n}/sum(occ{n}(:,1));       
            
        % Xnum: Gaussian kernel estimators for position
        if 1
            num_allspikes_e = 0;
            posind_allspikes_e = [];
            for vv = 1:maxtet
                if ~isnan(numspikes_e{n}(vv))
                    num_allspikes_e = num_allspikes_e + numspikes_e{n}(vv);
                    posind_allspikes_e = [posind_allspikes_e ; posind_spike_e{n}{vv}(:)];
                end 
            end
            Xnum_all = normpdf( xbins{n}'*ones(1,num_allspikes_e),...           %  # of encoding spikes this epoch
                repmat(dirpos2{n}(posind_allspikes_e),num_xbins(n),1),...    %  # pos inds of all encoding spikes this epoch
                sxker );            
        elseif 0
            Xnum_all = normpdf( xbins{n}'*ones(1,length(spktimes_all{n})),...
                repmat(dirpos2{n}(posind_spike_all{n}),num_xbins(n),1),...
                sxker );
        end
            % set non-occupied xbins indices to 0
            %Xnum_all(~validxbins{n},:) = 0;
            % normalize for each spike
            %Xnum_all = bsxfun(@rdivide,Xnum_all,sum(Xnum_all,1));
            %for skk = 1:size(Xnum_all,2)
            %    Xnum_all(:,skk) = Xnum_all(:,skk) / sum(Xnum_all(:,skk));
            %end
        % Lint: conditional intensity function for the unmarked case
        Lint_all{n} = sum(Xnum_all,2)  ./  occ{n}(:,1) ./  dt ; %integral
        % set non-occupied xbins indices to 0
        Lint_all{n}(~validxbins{n}) = 0;
        % normalize
        Lint_all{n} = Lint_all{n} ./ sum(Lint_all{n})  ;

    end
    
    clear Xnum_all
    
    %%% Sixth, decode (Bayes' Rule) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

    posteriors = cell(1,2);    
    xbins_cut = cell(1,2);
    
    % Set up 1-ms epoch vector
    ms_bin_timevec_ep = round(postimevec{1}(1)*1000):nummsbin:round(postimevec{1}(end)*1000);   % 1-ms time vector (from position time)
    num_msbin_bins_ep = length(ms_bin_timevec_ep);
    eventperiods = [];  % initialize

    if decodemode == 1              % 1: full epoch

        disp('Decoding full epoch continuously')
        
        eventperiods = [1 num_msbin_bins_ep]; 
        numevents = 1;
        eventtimes = [];
    end

    posteriors{1} = [];
    posteriors{2} = [];
    
    % Initialize millisecond time vector

    % Initialize spike matrix (for plotting, for event-by-event decoding)
    num_msbin_steps = num_msbin_bins_ep;

    % Iterate through 7 SD ripples (or events)
    for ww = 1:numevents
        
        % ms time of event (e.g. ripple) start
        ev_startms  = eventperiods(ww,1);

        % Initialize output matrices %%%%%%%%%%%%%%%%%%%%%%
        for nnn = 1:num_encodeeps
            
            % initialize
            postxM_r = zeros(num_xbins(nnn),num_msbin_steps) ;              % Matrix of posteriors
            SPK = round(spktimes_all_sorted{1}*1000);                   % Spike times in terms of 1-ms epoch vector
            
            % Uniform prior %
            clear postx ;
            uniformprior = ones(num_xbins(nnn),1) / sum(validxbins{nnn}) ;     % uniform prior
            uniformprior(logical(~validxbins{nnn}),:) = 0;
                uniformprior = uniformprior / sum(uniformprior);
                postx = uniformprior;
                
            tic
            % Iterate through each 1-ms time bin of the event
            ind_start = dec_start * 1000 / nummsbin + 1;            
            for t = ind_start:num_msbin_steps
                    
                % Prior probability
                if 1
                    prior = stateM_gausnorm{nnn} * postx;  % one step prediction, Markovian predication
                elseif 1
                    prior = uniformprior ;  % one step prediction, Markovian predication
                end
                    prior(~validxbins{nnn}) = 0;
                
                if all(isnan(prior))
                    keyboard
                end
                
                % Likelihood
                % initialized as flat
                L = ones( num_xbins(nnn), 1  );
                % obtain indices of spikes that occur closest to this 1-ms (epoch) time bin
                spind = find( SPK == ms_bin_timevec_ep(t + ev_startms - 1) ) ;
                
                num_spikesinbin = length(spind);
                
                if num_spikesinbin == 0
                    
                    L = exp(-Lint_all{nnn} .* dt) ;   %exp(-Lint_all.*dt) ;
                    L(~validxbins{nnn},:) = 0;
                    
                elseif num_spikesinbin > 0
                    
                    % initialize
                    l = zeros(num_xbins(nnn) , num_spikesinbin) ;
                    
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
                        l0 = normpdf(markAll{tet}(i,2)*ones(1,numspikes_e{nnn}(tet)),...  % ch. 1 amp (scalar)  * ones vector length of all spikes
                                     mark0{nnn}{tet}(:,1)',...                            % ch. 1 amplitude of all spikes
                                     smker).*...
                            normpdf(markAll{tet}(i,3)*ones(1,numspikes_e{nnn}(tet)),...   % ch. 2 " "
                                     mark0{nnn}{tet}(:,2)',...
                                     smker).*...
                            normpdf(markAll{tet}(i,4)*ones(1,numspikes_e{nnn}(tet)),...
                                    mark0{nnn}{tet}(:,3)',...
                                    smker).*...
                            normpdf(markAll{tet}(i,5)*ones(1,numspikes_e{nnn}(tet)),...
                                    mark0{nnn}{tet}(:,4)',...
                                    smker) ;
                        l1 = Xnum{nnn}{tet} * l0' ./ occ{nnn}(:,1) ./ dt ;   % divide by occupancy map
                        l2 = l1 * dt .* exp(-Lint{nnn}{tet} .* dt) ;  % exp factor
                        l2 = l2 ./ sum(l2(isfinite(l2))) ;                       % normalize the posterior
                        l(:,s) = l2 ;
                    end
                    
                    L = prod(l,2) ; % take product of likelihoods across spikes
                    L(logical(~validxbins{nnn}),:) = 0;
                    L = L./sum(L) ;
                    disp(sprintf('(dirpos) %s d%d ep%d: %s',animalname(1:3),d,ep,num2str(t * nummsbin / 1000)))
                end
                
                postx = prior .*  L ;
                postx = postx / sum(postx) ;  % normalize
                
                postxM_r(:,t) = postx;
                clear onestep a l L;
                
                % Print updates if doing full epoch / save halves
                if decodemode == 1
                    if t == round(1*num_msbin_steps/100)
                        disp('1% done')
                    elseif t == round(5*num_msbin_steps/100)
                        toc
                        disp('5% done')
                    elseif t == round(10*num_msbin_steps/100)
                        disp('10% done')
                    elseif t == round(20*num_msbin_steps/100)
                        disp('20% done')
                    elseif t == round(23*num_msbin_steps/100)
                        disp('21% done')
                        if 0
                            tvec = (1:num_msbin_steps)*nummsbin/1000;
                            decodeper = vec2list(sum(postxM_r,1) > 0,1:num_msbin_steps);
                            plot_start = decodeper(1)*nummsbin/1000;
                            plot_end = decodeper(end)*nummsbin/1000;
                                a = lookup(plot_start,tvec);
                                b = lookup(plot_end,tvec);
                                aa = lookup(plot_start,postimevec_epenc - postimevec_epenc(1));
                                bb = lookup(plot_end,postimevec_epenc - postimevec_epenc(1));
                            figure; imagesc(tvec(a:b),xbins{1},postxM_r(:,a:b),[0 .1]); hold on
                            colormap hot
                            %colormap(flipud(hot))
                            set(gca,'ydir','normal')
                            plot(postimevec_epenc(aa:bb) - postimevec_epenc(1),dirpos{1}(aa:bb),'b-','linewidth',2);
                            keyboard
                        end                        
                    elseif t == round(50*num_msbin_steps/100)
                        disp('50% done')
                        cd(savedir)
                        disp('done calculating half of posterior')
                        savefilename = sprintf('%s_%d_%d_%s_%d_dirpos_firsthalf',animalname(1:3),d,ep,'fullepoch',nnn);
                        P = struct;
                        P.nummsbin = nummsbin;
                        P.animalname = animalname;
                        P.dayep = [d ep];
                        P.remoteep = remoteep;
                        P.selected_tets = selected_tets;
                        P.spikethresh = spikethresh;
                        P.posvec = dirpos{1};   % actual animal traj / positions
                        P.stateM_gausnorm = stateM_gausnorm;
                        P.linposbins{1} = (1:length(xbins_cut{1})) * xdel ;
                        P.validxbins = validxbins;
                        P.posteriors{nnn} = postxM_r;
                        save(savefilename,'P','-v7.3')
                        if 0
                            figure; 
                             plot(postimevec_epenc - postimevec_epenc(1),dirpos{1},'Color',[.8 .8 .8],'linewidth',2);
                        end

                    elseif t == round(90*num_msbin_steps/100)
                        disp('90% done')
                    end
                end
                
            end
            
            % if full epoch, save single continuous posterior to file
            % before proceeding to remote epoch

            cd(savedir)
            disp('done calculating full epoch posterior')
            savefilename = sprintf('%s_%d_%d_%s_%d_dirpos',animalname(1:3),d,ep,'fullepoch',nnn);
            P = struct;
            P.nummsbin = nummsbin;
            P.animalname = animalname;
            P.dayep = [d ep];
            P.remoteep = remoteep;
            P.selected_tets = selected_tets;
            P.spikethresh = spikethresh;
            P.posvec = dirpos{1};   % actual animal traj / positions
            P.stateM_gausnorm = stateM_gausnorm;
            P.linposbins{1} = (1:length(xbins_cut{1})) * xdel ;
            P.validxbins = validxbins;
            P.posteriors{nnn} = postxM_r;
            save(savefilename,'P','-v7.3')
            clear P postxM_r_cut

        end
        
    end
    

    
end
        
        
      






end