function kk_clusterlessdecode(directoryname,fileprefix,dayeps,animalname,varargin)
%
% kk_linear2dayprocess(directoryname,fileprefix,days, options)
%
%  W-track specific day process (C, L, R)
%       Relies on linpos from the original kk_lineardayprocess
%
%
%Runs linearizeposition for all run epochs in each day and saves the data in
%'linpos' in the directoryname folder.  See LINEARIZEPOSITION for the definitions 
%of the options.
         
animalinfo = animaldef(animalname);
blank = '';
tetfilter = '';
savedir = '/opt/data13/kkay/Superposteriors_data';
decodemode = 1;  % default is all epoch
modelnum = 1;
spikethresh = 0;
mdel = 10;  % in uV
xdel = 1;  % in cm
dt = .001;
sigma_transmat = 1;
old_version = 0;
extratime = 500;  % ms of extra time before and after ripple to decode
plot_powertraces = 0;
calctraj = 0;

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'savedir'
            savedir = varargin{option+1};           
        case 'decodemode'
            decodemode = varargin{option+1};            
        case 'modelnum'
            modelnum = varargin{option+1};           
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
        case 'tetfilter'
            tetfilter = varargin{option+1};
        case 'spikethresh'
            spikethresh = varargin{option+1};
        case 'extratime'
            extratime = varargin{option+1};        
        case 'old_version'
            old_version = varargin{option+1};        
        case 'plot_powertraces'
            plot_powertraces = varargin{option+1};  
        case 'calctraj'
            calctraj = varargin{option+1};                
    end
end

%days = days(:)';
days = dayeps(1,:);

for de = 1:size(dayeps,1)
    
    day = dayeps(de,1);
    ep = dayeps(de,2);
    
   % Load data
   task = loaddatastruct(directoryname,fileprefix,'task',day);
   linpos = loaddatastruct(directoryname,fileprefix,'linpos',day);
   pos = loaddatastruct(directoryname,fileprefix,'pos',day);
   tetinfo = loaddatastruct(directoryname,fileprefix,'tetinfo',day);  
   trajencode = loaddatastruct(directoryname,fileprefix,'trajencode',day);
   if plot_powertraces
       out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'rippletrace', day, ep);
       riptrace = zscorer(out{day}{ep}{1}.powertrace);
       riptrace_timevec = out{day}{ep}{1}.eegtimesvec_ref;
       out = loadtracestruct(animalinfo{2}, animalinfo{3}, 'fastgammatrace', day, ep);
       wgtrace = zscorer(out{day}{ep}{2}.powertrace);
       wgtrace_timevec = out{day}{ep}{2}.eegtimesvec_ref;
   end
   
   % Verify epoch
    if (ep > length(task{day})) || isempty(task{day}{ep}) || ~isfield(task{day}{ep},'type') || ~strcmp(task{day}{ep}.type,'run')
        continue
    end   
   
   % Identify tetrodes
   dummy = evaluatefilter(tetinfo,tetfilter);
        selected_tets = unique(dummy((dummy(:,1) == day),3))';    
    
   % Identify 2 SD SWRs (here, to use to exclude from encoding model) 
   ripout = kk_getconstimes(animalinfo{2},animalinfo{3},[day ep],'ripplescons',1,...
       'consensus_numtets',3,'minthresh',2,...
       'exclusion_dur',0,'minvelocity',0,'maxvelocity',4);
   consvec_rip2 = ripout{day}{ep}.cons;
   consvectimes_rip2 = ripout{day}{ep}.time;
   periodtimes_rip2 = vec2list(consvec_rip2,consvectimes_rip2);
       
   % Process Epoch
   if (~isempty(task{day}{ep})) && ~isfield(task{day}{ep},'type')
       
       disp(sprintf('this epoch (d%d e%d) has no .type field!',day,ep))
       
   elseif ( (~isempty(task{day}{ep})) && (strcmp(task{day}{ep}.type,'run')) && ...
           (~strcmp(task{day}{ep}.environment,'OpenFieldA')) && ...
           (~strcmp(task{day}{ep}.environment,'OpenFieldB')) && ...
           (~strcmp(task{day}{ep}.environment,'LinearA'))   )
       
       disp(['Day ',num2str(day), ', Epoch ',num2str(ep)])
       
       % Basic epoch data
       postimevec = linpos{day}{ep}.statematrix.time;
            numpossamps = length(postimevec);
       starttime_ep = postimevec(1);
       endtime_ep = postimevec(end);
       
       % Identify W-track segment indices
       lindist = linpos{day}{ep}.statematrix.lindist;
       seg1 = linpos{day}{ep}.statematrix.segmentIndex==1;
       seg2 = linpos{day}{ep}.statematrix.segmentIndex==2;
       seg3 = linpos{day}{ep}.statematrix.segmentIndex==3;
       seg4 = linpos{day}{ep}.statematrix.segmentIndex==4;
       seg5 = linpos{day}{ep}.statematrix.segmentIndex==5;
       
       % Initialize position variables : arm dists %%%%%%%%%%%%%%%%%%%
       armdists = nan(1,numpossamps) ;   % (horizontal vector)   +200 / +400 / +600 linear distance from center well
       armdists2 = nan(numpossamps,1) ;  % (vertical vector)         
            centerarmmax =  max( lindist(seg1) );   % maximum linear distance of center arm
            rightarmmax  =  max( lindist(seg4 | seg5)    - centerarmmax );
            armdists(seg1)        = lindist(seg1)                           +   200;  % >200: center arm
            armdists(seg4 | seg5) = lindist(seg4 | seg5)    - centerarmmax  +   400;  % >400: right arm
            armdists(seg2 | seg3) = lindist(seg2 | seg3)    - centerarmmax  +   600;  % >600: left arm
                armdists = armdists(:);
                armdists2 = armdists(:)';
       armdists_cat = nan(size(armdists));  % position vectors with no spacing
            armdists_cat(seg1) = lindist(seg1)  ;
            armdists_cat(seg4 | seg5) = lindist(seg4 | seg5);
            armdists_cat(seg2 | seg3) = lindist(seg2 | seg3) + rightarmmax;
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       %%%% Construct positional mark space (xbins) + transition matrix (stateM) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       xbins = min(armdists):xdel:(max(armdists) + xdel);        %    % -3 : 0.1 : 3
             num_xbins = length(xbins);
                %indices for each positional "island" a to b, c to d, e to f
                a_cut = 1;
                b_cut = lookup(max(armdists(seg1)),xbins);
                c_cut = lookup(min(armdists(seg4 | seg5)),xbins);
                d_cut = lookup(max(armdists(seg4 | seg5)),xbins);
                e_cut = lookup(min(armdists(seg2 | seg3)),xbins);
                f_cut = lookup(max(armdists(seg2 | seg3)),xbins);
       validxbins = zeros(size(xbins));
       validxbins([a_cut:b_cut c_cut:d_cut e_cut:f_cut]) = 1 ;   % indicates which xbins indices were actually occupied by animal
            validxbins = logical(validxbins);
            
       %  Construct the transition matrix from animal's behavior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       stateM = nan(num_xbins,num_xbins);
            [ ~ , ind_xbin ]  =  histc(armdists,xbins);
                xbinstate  =  [ ind_xbin(1:end-1)   ind_xbin(2:end) ];  % [ <current bin #>  <next bin #>  ]
       % Iterate through each linear position
       for lp = 1:num_xbins
           % nextstates : for a given bin, collection of next bin # of the animal's linear position from actual data
           next_xbin = xbinstate( xbinstate(:,1)==lp , 2  );
           if ~isempty(next_xbin)
               % histogram and normalizes the distribution
               stateM(:,lp) = histc( next_xbin , 1:num_xbins ) / length(next_xbin)  ;
           elseif isempty(next_xbin)
               stateM(:,lp) = zeros(1,num_xbins);
           end
       end
       
       %    Smooth the transition matrix 
       if 0
           K = inline('exp(-(x.^2+y.^2)/2/sig^2)');                      % Gaussian
           [dx,dy] = meshgrid( [-2:2] );
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
           stateM_gausnorm = stateM_gaus * diag(1./sum(stateM_gaus,1));  % Normalized to confine each bin (columns) probability to 1
           stateM_gausnorm(isnan(stateM_gausnorm)) = 0; %0.001;         
       elseif 1
           gkernel = gaussian2(1,8);
           stateM_gaus = conv2(stateM,gkernel,'same');
           stateM_gausnorm = stateM_gaus * diag(1./sum(stateM_gaus,1));
           stateM_gausnorm(isnan(stateM_gausnorm)) = 0.0001;
       else
           stateM_gausnorm = stateM;
       end
       
       if 0
           %keyboard
           figure; imagesc(stateM_gausnorm)
       end
       
       %stateM_gausnorm = stateM_gausnorm + .0001;
       
       %%%%%%  Collect spike data + Create encoding model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       mark0 = {};       % mark0_t1
       spktimes_d = {};    % time_t1
       spktimes_e = {};   % time2_t1
       markAll = {};
       Xnum = {};
       Lint = {};
       
       for tet = selected_tets
           
           daydir = getdaydir(animalname);
           filedata = loadparamsfile(daydir,day,tet);
           
           % Identify spikes for encoding
           
           %  Spikes in experimenter-transcribed (in notebook) epoch
           inds_ep =  (  filedata.params(:,1)/10000  >=  starttime_ep  )  &  ( filedata.params(:,1)/10000 <= endtime_ep );
           %  Spikes that has spikethresh (uV) voltage amplitude in at least one channel
           inds_thresh = any(filedata.params(:,2:5) > spikethresh, 2) ;
           %  Spikes in trajencode (non-low speed and non-well)
                trajencode_periods = vec2list(trajencode{day}{ep}.trajstate > 0,trajencode{day}{ep}.timevec);
           inds_trajencode = isExcluded(filedata.params(:,1)/10000,trajencode_periods)  ;
%            % inds_vel :      Spikes @ > 6 cm/s
%                 spikeposinds = lookup(filedata.params(:,1)/10000,postimevec);
%                 spikevels = pos{day}{ep}.data(spikeposinds,5);
%            inds_vel = spikevels > 6;
           % inds_swr :      Spikes that occur within SWR (2 SD) periods
           inds_swr = isExcluded(filedata.params(:,1)/10000, periodtimes_rip2)  ;
           
           % Determine final set of inds
           if modelnum == 1
               % basic case
               inds_e = inds_ep & inds_thresh;
               inds_d = inds_ep & inds_thresh;
           elseif modelnum == 2
               % exclusion case               
               % inds4 = isExcluded(filedata.params(:,1)/10000,gfperiods);
               inds_e = inds_ep & inds_thresh & inds_trajencode;
               inds_d = inds_ep & inds_thresh;
               disp(sprintf('%d of %d spikes excluded from encoding',sum(inds_d)-sum(inds_e),sum(inds_d)))
           elseif modelnum == 3
               inds_e = inds_ep & inds_thresh & ~inds_swr;
               inds_d = inds_ep & inds_thresh;
               disp(sprintf('%d of %d spikes excluded from encoding',sum(inds_d)-sum(inds_e),sum(inds_d)))               
           else
               keyboard
           end
           
           % Spike times for E vs. D
           spktimes_e{tet} = filedata.params(inds_e,1)/10000;    % exclusive, to ENCODE
           spktimes_d{tet} = filedata.params(inds_d,1)/10000;           % inclusive, to DECODE           
                numspikes_e = length(spktimes_e{tet});
                
           % Amplitude mark vector (4 channel amplitudes) %%%%%%%%%%%%
           mark0{tet} = filedata.params(inds_e,2:5);
           
           % Obtain position indices for each spike
           [        ~ ,        posind_spike_d{tet} ]  =  histc(spktimes_d{tet},postimevec)   ;
           [        ~ ,        posind_spike_e{tet} ]  =  histc(spktimes_e{tet},postimevec)  ;    % ENCODING
           
           % Mark vector for decoding spikes %%
           markAll{tet}(:,1) = posind_spike_d{tet} ;         % Col 1: position time bin # , Col 2-5: spike amplitudes
           markAll{tet}(:,2:5) = filedata.params(inds_d,2:5);
           
           % Set up Amplitude mark space %%%%%%%%%%%%%%%%%%%%%%%%
               minamp = min(min(markAll{tet}(:,2:5)));
               maxamp = max(max(markAll{tet}(:,2:5)));
           mbins =  minamp : mdel : maxamp;
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%% Encoding model  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
           
           % actual animal positions
           X1  = repmat(xbins',1,numpossamps);         % position bin vectors repeated as columns
           MU1 = ones(length(xbins),1)  *  armdists2;  % actual positions of animal
           
           % encoding spike positions
           X2  = repmat(xbins',1,numspikes_e);                                % position vectors repeated as columns
           MU2 = ones(length(xbins),1)  *  armdists2(posind_spike_e{tet});    % actual (encoding) spike positions

               % 1.  occ: occupancy (linear position) map
               normpdf1 = normpdf(X1,MU1,sxker);
               occ = repmat(sum(normpdf1,2),1,length(mbins));  %normpdf1 * ones(numpossamps,length(mbins))   ;
                    occ(logical(~validxbins),:) = 0;
                    occ = occ/sum(occ(:,1));
               % 2.  Xnum: Gaussian kernel estimators for position, for each encoding spike
               Xnum{tet} = normpdf(X2,MU2,sxker);
                    % set non-occupied xbins indices to 0
                    Xnum{tet}(~validxbins,:) = 0;
                    % normalize for each spike
                    %Xnum{tet} = bsxfun(@rdivide,Xnum{tet},sum(Xnum{tet},1));
                    for sk = 1:size(Xnum{tet},2)
                        Xnum{tet}(:,sk) = Xnum{tet}(:,sk) / sum(Xnum{tet}(:,sk));
                    end
                    %keyboard
               % 3.  Lint: firing occ-norm position map for tetrode
               Lint{tet} = sum(Xnum{tet},2) ./ occ(:,1) ./ dt;   %integral (?)
                    % set non-occupied xbins indices to 0
                    Lint{tet}(~validxbins) = 0;
                    Lint{tet} = Lint{tet} ./ sum(Lint{tet}) ;  % normalization
               
           
       end
       

       % Decoding: collect all spikes and their marks %%%%%%%%%%%%%%%
       
            spktimes_all = [];          % all spike times
            spktimes_all_sorted = [];   % sorted spike times, sorted
            posind_spike_all = [];
            sinds = [];                 
            % iterate through tetrodes 
            for tet = selected_tets
                numspikestet_e        =  length(spktimes_d{tet});
                spktimes_all        =  [spktimes_all ;  spktimes_d{tet}  tet * ones(numspikestet_e,1) ];   % * tag with the tet number    
                posind_spike_all    =  [posind_spike_all ; posind_spike_d{tet}];
            end
            
            % Sort all decoding spikes in time
            [spktimes_all_sorted , sinds] = sort(spktimes_all(:,1));
            
            % Construct tet_ind "indicator matrix" : [  <spike #>  x  <which tetrode spikes> ]
            tet_ind = zeros(length(spktimes_all_sorted),max(selected_tets))  ;
            for i=1:length(spktimes_all)
                parenttet = spktimes_all(sinds(i),2) ;   % tet responsible for the spike
                tet_ind(i,parenttet) = 1 ;
            end
            tet_spikenum = tet_ind.*cumsum(tet_ind,1);  % [ <spike #>  x <index of spike for its parent tetrode>]
            
            numpostimes = size(postimevec,1);
            
            %%% Decoding spike marks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                % occ: columns are identical; occupancy based on position; denominator
                occ = normpdf( xbins'*ones(1,numpostimes),...
                               repmat(armdists2,num_xbins,1),...
                               sxker )   *    ones(numpostimes,length(mbins))    ;
                    % remove invalid rows and normalize
                    occ(logical(~validxbins),:) = 0;
                    occ = occ/sum(occ(:,1));                
                % Xnum: Gaussian kernel estimators for position
                Xnum_all = normpdf( xbins'*ones(1,length(spktimes_all)),...
                                    repmat(armdists2(posind_spike_all),num_xbins,1),...
                                    sxker);
                    % set non-occupied xbins indices to 0
                    Xnum_all(~validxbins,:) = 0;
                    % normalize for each spike
                    %Xnum_all = bsxfun(@rdivide,Xnum_all,sum(Xnum_all,1));
                    for skk = 1:size(Xnum_all,2)
                        Xnum_all(:,skk) = Xnum_all(:,skk) / sum(Xnum_all(:,skk));
                    end
                % Lint: conditional intensity function for the unmarked case
                Lint_all = sum(Xnum_all,2)  ./  occ(:,1) ./  dt ; %integral
                    % set non-occupied xbins indices to 0
                    Lint_all(~validxbins) = 0;
                    % normalize
                    Lint_all = Lint_all ./ sum(Lint_all)  ;

            
            %%%% DECODING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Set up 1-ms epoch vector
            mstimevec_ep = round(postimevec(1)*1000):1:round(postimevec(end)*1000);   % 1-ms time vector (from position time)
            nummsbins_ep = length(mstimevec_ep);
            
            eventperiods = [];  % initialize
            
            if decodemode == 1              % 1: full epoch
                
                eventperiods = [1 nummsbins_ep];
                numevents = 1;
                
            elseif decodemode == 2          % 2: SWRs
                
                % SWR parameters
                cons_name1 = 'ripplescons';
                tetfilter1 = 1;     % 1 for validripples CA1, 2 for CA3-DG
                consensus_numtets_ripc = 3;
                minthresh_ripc = 7;
                exclusion_dur_ripc = 0;
                minvelocity_ripc = 0;
                maxvelocity_ripc = 4;
               
                % retrieve ripples
                output1 = kk_getconstimes(animalinfo{2},animalinfo{3},[day ep],cons_name1,tetfilter1,...
                    'consensus_numtets',consensus_numtets_ripc,'minthresh',minthresh_ripc,...
                    'exclusion_dur',exclusion_dur_ripc,'minvelocity',minvelocity_ripc,'maxvelocity',maxvelocity_ripc);
                consvec = output1{day}{ep}.cons;
                consvectimes = output1{day}{ep}.time;
                riplist = vec2list(consvec,consvectimes);
                numevents = size(riplist,1);
                
                disp(sprintf('%d %s detected',size(riplist,1),cons_name1))
    
                ripplevec = list2vec(riplist,mstimevec_ep/1000);  % <time bin> x <consvec> 
                eventperiods = vec2list(ripplevec,1:nummsbins_ep); % ripple periods in 1-ms indices
               
            end


            
            
            
            % iterate through ripples (or events)
            for ww = 1:numevents
                
                % Event 1-ms times
                eventstart = eventperiods(ww,1);
                eventend = eventperiods(ww,2);
                
                % ms time bins to decode
                if decodemode == 1
                    mstimevec = eventstart:eventend;                 
                else
                    mstimevec = (eventstart - extratime ) : (eventstart + extratime ); 
                        
                end
                
                num_mssteps = length(mstimevec) ;
                
                % ms time of ripple start
                ev_startms = lookup(eventperiods(ww,1),mstimevec);
                ev_endms = lookup(eventperiods(ww,2),mstimevec);
                
                % recording file clock times for the event
                ev_starttime = starttime_ep + eventstart/1000;
                ev_endtime = starttime_ep + eventend/1000;
                
                % position of the animal during the event (should barely change)
                ev_startind = lookup(ev_starttime - extratime/1000, postimevec);
                ev_endind = lookup(ev_starttime + extratime/1000, postimevec);
                if ev_startind <= 0 || ev_endind > length(postimevec)
                    disp('cant plot event at edge of epoch')
                    continue
                end
                winpos = armdists_cat( ev_startind:ev_endind );  % armdists_cut (linear) positions over the course of the ripple
                winpos_timevec = 1000 * ( postimevec(ev_startind:ev_endind) - ev_starttime );  % 1-ms indices over the course of the ripple
                
                
                % Uniform prior %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                clear postx postxM_r ;
                postx = ones(num_xbins,1) / sum(validxbins) ;     % uniform prior
                    postx(logical(~validxbins),:) = 0;
                   
                % Initialize output matrices %%%%%%%%%%%%%%%%%%%%%%
                postxM_r = zeros(num_xbins,num_mssteps) ;              % Matrix of posteriors
                spike_mat = zeros(max(selected_tets),num_mssteps);     % SPIKE MATRIX : <tetrode #> x <time>  
                SPK = round(spktimes_all_sorted*1000);                   % Spike times in terms of 1-ms epoch vector
                
                tic
                % Iterate through each 1-ms time bin of the event
                for t = 1:num_mssteps
                    
                    tt = mstimevec(t);    % tt: recording clock time / epoch ms index of this decoding bin
                    
                    % Print updates if doing full epoch
                    if decodemode == 1
                       if t == round(1*num_mssteps/100)
                           disp('1% done')
                       elseif t == round(5*num_mssteps/100)
                           disp('5% done')
                       elseif t == round(10*num_mssteps/100)
                           disp('10% done')
                       elseif t == round(20*num_mssteps/100)
                           disp('20% done')
                       elseif t == round(50*num_mssteps/100)
                           disp('50% done')
                       elseif t == round(90*num_mssteps/100)
                           disp('90% done')
                       end
                    end
                    
                    % Prior probability
                    prior = stateM_gausnorm * postx;  % one step prediction, Markovian predication
                    
                    if all(isnan(prior))
                        keyboard
                    end
                    
                    % Likelihood
                        % initialized as flat
                    L = ones( num_xbins, 1  );       
                        % obtain indices of spikes that occur closest to this 1-ms (epoch) time bin
                    spind = find( SPK == mstimevec_ep(tt) ) ;  
                    
                    num_spikesinbin = length(spind);
                    if num_spikesinbin == 0
                        L = exp(-Lint_all .* dt) ;   %exp(-Lint_all.*dt) ;
                          L(~validxbins,:) = 0;
                    elseif num_spikesinbin > 0
                        
                        % initialize
                        l = zeros(num_xbins , num_spikesinbin) ;   
                        
                        % iterate through spikes (from all tetrodes) contained in this bin
                        for s = 1:num_spikesinbin
                            ind = spind(s);    % spike index in spktimes_all_sorted
                            tet = find(tet_ind(ind,:)) ;
                            numspikestet_e = length(spktimes_e{tet});
                            % register in spike matrix (for raster plot later)
                            spike_mat(tet,t) = 1 ;    
                            % 
                            i = tet_spikenum(ind,tet);   % "running" spike # for its parent tetrode
                            
                            %l0: amplitude mark likelihoods for this single spike
                            l0 = normpdf(markAll{tet}(i,2)*ones(1,numspikestet_e),...  % ch. 1 amp (scalar)  * ones vector length of all spikes
                                             mark0{tet}(:,1)',...                            % ch. 1 amplitude of all spikes
                                             smker).*...         
                                     normpdf(markAll{tet}(i,3)*ones(1,numspikestet_e),...   % ch. 2 " "
                                             mark0{tet}(:,2)',...
                                             smker).*...
                                     normpdf(markAll{tet}(i,4)*ones(1,numspikestet_e),...
                                             mark0{tet}(:,3)',...
                                             smker).*...
                                     normpdf(markAll{tet}(i,5)*ones(1,numspikestet_e),...
                                             mark0{tet}(:,4)',...
                                             smker) ;
                            l1 = Xnum{tet} * l0' ./ occ(:,1) ./ dt ;   % divide by occupancy map
                            l2 = l1 * dt .* exp(-Lint{tet} .* dt) ;  % exp factor
                            l2 = l2 ./ nansum(l2) ;                       % normalize the posterior
                            l(:,s) = l2 ;
                        end

                        L = prod(l,2) ; % take product of likelihoods across spikes
                            L(logical(~validxbins),:) = 0;
                                L = L./sum(L) ;
                    end        
                    
                    postx = prior .*  L ;
                        postx = postx / sum(postx) ;  % normalize

                        
                    postxM_r(:,t) = postx;
                    clear onestep a l L;
                end
                toc
               
                xvecms = (1:num_mssteps)-ev_startms;
                ytetvec = 1:max(selected_tets);
                
                % some additional processing of the posteriors %%%%%%%%%
                posteriors = postxM_r;
                xbins_cut = xbins - 200;
                xbins_cut(d_cut:e_cut) = [];
                xbins_cut(b_cut:c_cut) = [];
                    linposbins = (1:length(xbins_cut)) * xdel ;
                posteriors(d_cut:e_cut,:) = [];
                posteriors(b_cut:c_cut,:) = [];
                
                
                if decodemode == 1
                    cd(savedir)
                    savefilename = sprintf('%s_%d_%d',animalname(1:3),day,ep);
                    disp('done calculating full epoch posterior')
                    save(savefilename,'posteriors','xvecms','linposbins','spike_mat','selected_tets')
                    return
                end
                
                
                
                %%% Figure plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 
                figure; hold on
                
                % Spike raster plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                subplot(3,2,[1 3 5]);
                
                    % plot 2 SD SWRs that occur in the widnow %%%%%%%%%
                maxtet = max(selected_tets);
                consvec = output1{day}{ep}.cons;
                c = lookup(ev_starttime - 0.5,consvectimes);      
                d = lookup(ev_starttime + 0.5,consvectimes); 
                if 0
                    plot(1000 * (consvectimes_rip2(c:d) - rip_starttime), consvec_rip2(c:d) + 0.5,'-','Color',[1 .7 .7],'linewidth',8);      hold on          
                else
                    % identify which ripples to plot
                    plotwinvec = logical(list2vec([ev_starttime-0.5   ev_starttime+0.5],ripout{day}{ep}.time))';
                    consvec_rip2_toplot = consvec_rip2 & plotwinvec;
                    rip_toplot = 1000 * (vec2list(consvec_rip2_toplot,ripout{day}{ep}.time) - ev_starttime);
                    numrip_toplot = size(rip_toplot,1);
                    % plot all 2 SD ripples in the raster window
                    for rp = 1:numrip_toplot
                        patch([rip_toplot(rp,1) rip_toplot(rp,2) rip_toplot(rp,2) rip_toplot(rp,1)],...
                              [1 1 (maxtet + 4) (maxtet + 4)],...
                              [1 .9 .9],'edgecolor','none'); hold on
                    end
                end
                    % plot EVENT as a thick red line (i.e. ripple or wave gamma) %%%%%%%%%%%
                hold on
                plot([ev_startms ev_endms] - ev_startms,[.5 .5],'Color','r','linewidth',4)
                if 0
                    set(gca,'xtick',[-300:100:300]);
                    set(gca,'xticklabel',{'-300','','','0','','','+300'});
                    xlim([-300 300]);
                else
                    set(gca,'xtick',[-500:100:500]);
                    set(gca,'xticklabel',{'-500','','','','','0','','','','','+500'});
                    xlim([-500 500]);
                end                
               
                % plot SPIKES raster %% 
                for tet = 1:max(selected_tets)
                    spikebins = find(spike_mat(tet,:) > 0);
                    numlines = length(spikebins);
                    for ll = 1:numlines
                        plot([spikebins(ll) spikebins(ll)] - ev_startms,[tet-0.5 tet+0.5],...
                            'linestyle','-','Color',[0 0 0],'LineWidth',1)
                    end
                end
                set(gca,'ytick',1:max(selected_tets))
                %set(gca,'YTick',1:5:30,'YTickLabel',1:5:30,'FontSize',14);
                xlabel('Time (ms)') ;   ylabel('Tetrode');
                title([cons_name1 ' # : ',num2str(ww)],'fontweight','bold','fontsize',14);
                set(gca,'tickdir','out');
                
                box off
                
                % plot WG and SWR power traces at bottom%%%%%%%%%%%%%%%%
                a = lookup(ev_starttime - 0.5,riptrace_timevec) ;
                b = lookup(ev_starttime + 0.5,riptrace_timevec) ;
                plot(1000 * (riptrace_timevec(a:b) - ev_starttime), -riptrace(a:b)/2 + maxtet + 4,'r-','linewidth',2); hold on
                a = lookup(ev_starttime - 0.5,wgtrace_timevec) ;
                b = lookup(ev_starttime + 0.5,wgtrace_timevec) ;                
                plot(1000 * (wgtrace_timevec(a:b) - ev_starttime), 2 * -wgtrace(a:b)/2 + maxtet + 4,'-','Color',[0 .5 1],'linewidth',2);
                
                set(gca,'ydir','reverse')
                set(gca,'fontsize',12)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                % Decoded posterior image plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                subplot(3,2,[2 4 6]);
                imagesc(xvecms,linposbins,posteriors);

                %title('postx, clusterless','FontSize',12);
                ylabel('Linearized position','FontSize',12);
                %set(gca,'FontSize',14,'YTick',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3],'YTickLabel',[-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3]);
                %colormap(flipud(hot(256)));
                colormap(hot);
                caxis([0 0.1]);
                ylim([-.01 max(linposbins)])
                set(gca,'tickdir','out');

                    set(gca,'xtick',[-500:100:500]);
                    set(gca,'xticklabel',{'-500','','','','','0','','','','','+500'});
                    xlim([-500 500])
                set(gca,'fontsize',12)
                xlabel('Time (ms)')
                
                hold on
                % plot current position (dark grey line)
                plot(winpos_timevec,winpos,'-','linewidth',6,'Color',[.8 .8 .8])
                
                % plot line indicating start and end of ripple
                %plot([ripms_start ripms_end] - ripms_start,[50 50],'Color','r','linewidth',4)
                plot([0 0],[0 max(linposbins)],'Color','r','linewidth',1)
                plot([ev_endms - ev_startms ev_endms - ev_startms],...
                     [0 max(linposbins)],'Color','r','linewidth',1)
                
                % plot lines demarcating positions of arms
                    % between Center and Right
                plot([-500 500],[centerarmmax centerarmmax],'--','Color','w','linewidth',2)
                    % between Right and Left
                plot([-500 500],[centerarmmax+rightarmmax ...
                                 centerarmmax+rightarmmax],'--','Color','w','linewidth',2)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %keyboard
                %close all
            end
            
            
      
   end
   
end
      
end










