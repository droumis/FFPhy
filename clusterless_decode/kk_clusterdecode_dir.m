function kk_clusterdecode_dir(directoryname,fileprefix,dayeps,animalname,varargin)
%
% kk_clusterdecode_dir(directoryname,fileprefix,days, options)
%
%  W-track specific day process (C, L, R)
%       Relies on linpos from the original kk_lineardayprocess

         
animalinfo = animaldef(animalname);
   daydir = getdaydir(animalname);
blank = '';
tetfilter = '';
savedir = '/opt/data13/kkay/Superposteriors_data';
decodemode = 1;  % default is all epoch
modelnum = 1;
remoteW = 0;
extratime = 500;  % ms of extra time before and after ripple to decode
plot_powertraces = 0;
plot_infunction = 0;
xdel = 1;
xkern = 2;  % # of xdel bins to smooth occupancy / firing rate maps
winsize = 0.5;
min_activecells = 3;
manual_period = [];
minplaceunits = 10;

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'savedir'
            savedir = varargin{option+1};    
        case 'plot_infunction'
            plot_infunction = varargin{option+1};
        case 'decodemode'
            decodemode = varargin{option+1};            
        case 'modelnum'
            modelnum = varargin{option+1};  
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
        case 'minplaceunits'
            minplaceunits = varargin{option+1};
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
           if numunits(1) < minplaceunits && numunits(2) < minplaceunits
               disp('too few units')
               continue
           end
   end
   
   % Initialize outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
      % positional variables
   xbins = cell(1,2);
   armdists = cell(1,2);        % two elements since up to two different W-tracks
   armvect = cell(1,2);  
   validxbins = cell(1,2);   
   occ = cell(1,2);
   occ_sm = cell(1,2);
   centerarmmax = [nan nan];
   rightarmmax = [nan nan];
   maxarmdist = [nan nan];
        lastcenterbin = [nan nan];
        lastrightbin = [nan nan];
   numbinsingle = [nan nan];     
   
      % spike data
   spk_e = cell(1,2);
    
   % internal use (not to save in output file later)
   armdists2 = cell(1,2);
   armdists_cat = cell(1,2);
      a_cut = cell(1,2);
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
       armdists{nn} = nan(1,numpossamps_encep) ;   % (horizontal vector) +200 / +400 / +600 linear distance from center well         
       armvect{nn} = nan(1,numpossamps_encep);     %   dists + direction +200 / +400 / +600  << head dir = AWAY FROM center junction
                                                   %                +1200 / +1400 / +1600 << head dir = TOWARD center junction
       
        % Armdists (pseudo 1D)                                 
            centerarmmax(nn) =  max( lindist(seg1) );                                   % maximum linear distance of center arm
            rightarmmax(nn)  =  max( lindist(seg4 | seg5) ) - centerarmmax(nn) ;         % right arm, distance from center junction to well
            maxarmdist(nn) =  max ( lindist(seg2 | seg3) ) + rightarmmax(nn)  ;
                lastcenterbin(nn) = ceil(centerarmmax(nn));
                lastrightbin(nn) = ceil(rightarmmax(nn))+ceil(centerarmmax(nn));            
            armdists{nn}(seg1)        = lindist(seg1)                               +   200;  % >200: center arm
            armdists{nn}(seg4 | seg5) = lindist(seg4 | seg5)    - centerarmmax(nn)  +   400;  % >400: right arm
            armdists{nn}(seg2 | seg3) = lindist(seg2 | seg3)    - centerarmmax(nn)  +   600;  % >600: left arm
                 armdists{nn} = armdists{nn}(:);                
                                                   
                                                   
       % Retrieve head direction that is in reference to well 1 (the center well)                                            
       headdir = linpos{day}{epenc}.statematrix.segmentHeadDirection(:,1);   % head direction relative to center well -- values > 0 are outbound
                                                                                                 %                       values < 0 are inbound
            postimevec_nonnan = postimevec{nn}(~isnan(headdir));
            headdir_nonnan = headdir(~isnan(headdir));
            headdir2 = interp1(postimevec_nonnan,headdir_nonnan,postimevec{nn},'linear');   
       
                % (optional) plot head direction to validate the sign of headdir
                       % what we expect is that > 0 is W-track "outbound" , and < 0 is W-track "inbound"
                if 0
                    % head angle
                    headang = linpos{day}{epenc}.statematrix.headdir;
                    xlen = 5 * cos(headang);
                    ylen = 5 * sin(headang);
                    time_a = 415;
                    time_b = 418;
                    sttime = pos{day}{epenc}.data(1,1);
                    figure;
                    a = lookup(time_a,postimevec{nn}-postimevec{nn}(1));
                    b = lookup(time_b,postimevec{nn}-postimevec{nn}(1));
                    subplot(4,1,1:3)
                    plot(pos{day}{epenc}.data(:,2),pos{day}{epenc}.data(:,3),'.','Color',[.7 .7 .7]); hold on
                    plot(pos{day}{epenc}.data(a:b,2),pos{day}{epenc}.data(a:b,3),'-','linewidth',3,'Color',[0 0 0]);
                        % headang
                        for ind = a:b
                            % arrow
                            plot([pos{day}{epenc}.data(ind,2)  pos{day}{epenc}.data(ind,2) + xlen(ind) ],...
                                 [pos{day}{epenc}.data(ind,3)  pos{day}{epenc}.data(ind,3) + ylen(ind)],'-','linewidth',1,'Color',[0 0 0]);
                            % arrowhead
                            plot(pos{day}{epenc}.data(ind,2) + xlen(ind) ,...
                                 pos{day}{epenc}.data(ind,3) + ylen(ind),'x','linewidth',2,'Color',[1 0 0]);
                        end
                    subplot(4,1,4)
                    plot(pos{day}{epenc}.data(:,1)-sttime,headdir2,'.','Color',[.7 .7 .7]); hold on
                    plot(pos{day}{epenc}.data(a:b,1)-sttime,headdir2(a:b),'-','linewidth',3,'Color',[0 0 0]);
                end            
            
                times1 = seg1 & (headdir2 >= 0);   % center arm,  headed TO center junct
                times2 = seg1 & (headdir2 <  0);   % center arm,  headed AWAY FROM center junct
                times3 = (seg4 | seg5) & (headdir2 >= 0);  % right arm, headed AWAY FROM center junct
                times4 = (seg4 | seg5) & (headdir2 <  0);  % right arm, headed TO center junct   
                times5 = (seg2 | seg3) & (headdir2 >= 0);  % left arm, headed AWAY FROM center junct   
                times6 = (seg2 | seg3) & (headdir2 <  0);  % left arm, headed TO center junct   
                
                
                % check segments
                if 0
                    sttime = pos{day}{epenc}.data(1,1);
                    figure;       
                    plot(pos{day}{epenc}.data(:,1)-sttime,headdir2,'.','Color',[.7 .7 .7]); hold on
                    if 1
                        plot(pos{day}{epenc}.data(times1,1)-sttime,headdir2(times1),'.','markersize',10,'Color',[1 0 0]);  
                        plot(pos{day}{epenc}.data(times2,1)-sttime,headdir2(times2),'.','markersize',10,'Color',[0 0 1]); 
                        plot(pos{day}{epenc}.data(times3,1)-sttime,headdir2(times3),'.','markersize',10,'Color',[1 .5 0]);  
                        plot(pos{day}{epenc}.data(times4,1)-sttime,headdir2(times4),'.','markersize',10,'Color',[0 .9 .4]); 

                         plot(pos{day}{epenc}.data(times5,1)-sttime,headdir2(times5),'.','markersize',10,'Color',[0 1 0]);  
                        plot(pos{day}{epenc}.data(times6,1)-sttime,headdir2(times6),'.','markersize',10,'Color',[1 .4 1]);                        
                    end
                end
                
            
            % Now convert headdir to Wu-Foster-2014's "inbound" and "outbound", which are in reference to the center junction
                % now, "INBOUND" is defined as TOWARD center junction (200, 400, 600)
                % now, "OUTBOUND" is defined as AWAY FROM center junction (1200, 1400, 1600)
            armvect{nn}( times1 )   =   lindist( times1 )            +   200;       % >200: center arm + toward center junction
            armvect{nn}( times2 )   =   lindist( times2 )            +   1200;          % >1200: center arm + away from center junction
            
            armvect{nn}( times3 )   =   lindist(times3)    - centerarmmax(nn)  +   1400;  % >1400: right arm + away from center junction
            armvect{nn}( times4 )   =   lindist(times4)    - centerarmmax(nn)  +   400;  % >400: right arm + toward center junction

            armvect{nn}( times5 )   =   lindist(times5)    - centerarmmax(nn)  +   1600;  % >1600: left arm + away from center junction
            armvect{nn}( times6 )   =   lindist(times6)    - centerarmmax(nn)  +   600;  % >600: left arm + toward center junction
             
                

                
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       %%%% Construct positional space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % armvect %%%%%%%%%%%%%%
       xbins{nn} = 0:xdel:2000;      
       %xbins{nn} = min(armvect{nn}):xdel:(max(armvect{nn}) + xdel);      
       
       %indices for each positional "island" a to b, c to d, e to f
       a_cut{nn} = lookup(min(armdists{nn}(seg1)),xbins{nn},-1);
       b_cut{nn} = lookup(max(armdists{nn}(seg1)),xbins{nn},-1);
       c_cut{nn} = lookup(min(armdists{nn}(seg4 | seg5)),xbins{nn},-1);
       d_cut{nn} = lookup(max(armdists{nn}(seg4 | seg5)),xbins{nn},-1);
       e_cut{nn} = lookup(min(armdists{nn}(seg2 | seg3)),xbins{nn},-1);
       f_cut{nn} = lookup(max(armdists{nn}(seg2 | seg3)),xbins{nn},-1);
       
       validxbins{nn}([a_cut{nn}:b_cut{nn} ...
                       c_cut{nn}:d_cut{nn} ...
                       e_cut{nn}:f_cut{nn} ...
                       (a_cut{nn}:b_cut{nn}) + 1000 ...
                       (c_cut{nn}:d_cut{nn}) + 1000 ...
                       (e_cut{nn}:f_cut{nn}) + 1000 ...
                       ]) = 1 ;   % indicates which xbins indices were actually occupied by animal
            disp(sprintf('%d bins valid',sum(validxbins{nn})))
            validxbins{nn} = logical(validxbins{nn});  
                numbinsingle(nn) = sum(validxbins{nn})/2;
       

       % armdists_cat (used for plotting and output later -- same as armdists but with no 200-400-600 buffer) 
       armdists_cat{nn}              = nan(1,numpossamps_encep) ;
       armdists_cat{nn}(seg1)        = lindist(seg1)  ;  % 
       armdists_cat{nn}(seg4 | seg5) = lindist(seg4 | seg5)    ;  %
       armdists_cat{nn}(seg2 | seg3) = lindist(seg2 | seg3)   +   rightarmmax(nn) ;  % 
       
       % armvect_cat (used for plotting and output later -- same as armvect but with no 200-400-600 buffer) 
       armvect_cat{nn}                = nan(1,numpossamps_encep) ;
            % headed toward center junct
       armvect_cat{nn}(times1)        = lindist(times1)  ;                                          % 
       armvect_cat{nn}(times4)        = lindist(times4)    ;                                        %
       armvect_cat{nn}(times6)        = lindist(times6)  + rightarmmax(nn) ;   %
            % headed away from center junct
       armvect_cat{nn}(times2)        = lindist(times2)                    + numbinsingle(nn) ;                                          % 
       armvect_cat{nn}(times3)        = lindist(times3)                    + numbinsingle(nn)  ;                                        %
       armvect_cat{nn}(times5)        = lindist(times5)  + rightarmmax(nn) + numbinsingle(nn)  ;   %       
       
  
       

        % Calculate occupancy (position + direction, thus using armvect)
        if modelnum == 4
            validposinds = logical(isExcluded(postimevec{nn},hispeedperiods));
           armvect_filt = armvect{nn}( validposinds );
           occ{nn} = histc(armvect_filt,xbins{nn}) * 1/29.97;
        else
            occ{nn} = histc(armvect{nn},xbins{nn}) * 1/29.97;
        end   
            if 0
               figure;
                    % island plot
               subplot(1,2,1)
               plot(1:2001,occ{nn});
                    % stitch plot
               subplot(1,2,2)
               cutocc = occ{nn}(validxbins{nn});
               plot(1:length(cutocc),cutocc); hold on
                        % stitch dividers, inbounds
                        plot([lastcenterbin(nn) lastcenterbin(nn)],[0 max(cutocc)],'k--')
                        plot([lastrightbin(nn)  lastrightbin(nn)],[0 max(cutocc)],'k--')
                        plot([numbinsingle(nn)       numbinsingle(nn)],[0 max(cutocc)],'k--')
                        % stitch dividers, outbounds
                        OFFSET = numbinsingle(nn);
                        plot([lastcenterbin(nn) lastcenterbin(nn)] + OFFSET,[0 max(cutocc)],'g--')
                        plot([lastrightbin(nn)  lastrightbin(nn)] + OFFSET,[0 max(cutocc)],'g--')
                        plot([numbinsingle(nn)       numbinsingle(nn)] + OFFSET,[0 max(cutocc)],'g--')
            end       
       
       occ_sm{nn} = smoothvect(occ{nn},gaussian(xkern,8*xkern))';
            occ_sm{nn} = occ_sm{nn}(:);
            % ** Later, exclude bins with very low occupancy and replace with NaN
            lowocc = occ_sm{nn} < xdel*.1 ;  % Maggie used threshold of 0.1 s / cm               
            
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

       %%%%%%  Create encoding model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       spk_e{nn} = cell(1,numunits(nn));
       unitmap{nn} = cell(1,numunits(nn));
       
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
               encodevec = []; % would need to code this 
               disp(sprintf('%d of %d spikes excluded from encoding',sum(inds_d)-sum(inds_e),sum(inds_d)))
           elseif modelnum == 3
               % non-SWR
               inds_pre = inds_epenc;
               inds_e = inds_epenc & ~inds_swr;
                    % also, formulate encoding periods in a period format (necessary
                    % to create accurate occupancy map below)
                    encodeperiods{nn} = vec2list(~consvec_rip2,consvectimes_rip2);
                    %disp(sprintf('%d of %d spikes excluded from encoding',sum(inds_pre)-sum(inds_e),sum(inds_pre)))
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
           
           % Spike times (recording clock) for Encoding + Decoding
           spk_e{nn}{cc} = spiketimes(inds_e);  
           
           % Calculate firing rate map for each unit
           posinds = lookup(spiketimes(inds_e),postimevec{nn});
                spk_armvects = armvect{nn}(posinds); 
                    spk_bincounts = histc(spk_armvects,xbins{nn});
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
   
   
   %%% Second, for decoding data, collect amplitude mark data %%%%%%%%%%%%%%%%%%%%%%%%   
   
   ep;  % the decoding epoch #
   
   clear linpos
                  
   % Decode for each encoding epoch  (since unit common lists might be different,
                            %         doing two separate encoding epochs)
 
   binvec = [];
   binvec_c = [];

   S = cell(1,2);      % [ <times> <bin #> <cellnumber> ]
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
   for nn = 1:num_encodeeps
   
       % initialize these outputs (for later concatenation)
       S{nn} = [];
       spikecounts{nn} = [];
       
       % Go through each cell and calculate the binned spike counts
       for cc2 = 1:numunits(nn)
           
           tet2 = detc{nn}(cc2,3);
           cellnum2 = detc{nn}(cc2,4);
           spk = spikes{day}{ep}{tet2}{cellnum2}.data;   % decode all clustered spikes
           if ~isempty(spk)
               spk_d = spk(:,1);
           else
               spk_d = [];
           end
                
           % Assign each spike into a bin
           [~,spikebins] = histc(spk_d, binvec);
                validinds = spikebins > 0;
                spikebins = spikebins(validinds);
                    numspk = length(spikebins);
                    
           % Concatenate to celldata ( [ <spike times> <bin #> <cellnumber> ] ) and
           if numspk > 0
               tmpcelldata = [ spk_d(validinds)   spikebins(:)   cc2*ones(numspk,1)  ];
               S{nn} = [S{nn} ; tmpcelldata];
           end
                      
           % Register spikecounts: # of spikes from each cell in each bin (cell # x bin #)
           spikecount_unit = zeros(1,numbins);
           for i = 1:length(spikebins)
               spikecount_unit(spikebins(i)) = spikecount_unit(spikebins(i)) + 1;
           end
           
           spikecounts{nn} = [spikecounts{nn} ; spikecount_unit];
           
       end
       
        % Sort unit data
        S{nn} = sortrows(S{nn},1);
      
        % For each time bin, count # of distinct cells contributing at least one spike
        cellcounts{nn} = sum((spikecounts{nn} > 0),1);       
        
        % Find all decoding bins with enough cells
        activebins{nn} = find(cellcounts{nn} >= min_activecells);   % by bin index
        
       
   end
     
   
    %%% Sixth, decode (Bayes' Rule) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

    posteriors = cell(1,2);
    numcutlinpos = [nan nan];
    
    
    for nnn = 1:num_encodeeps
    
        % define these to make code cleaner
        numlinpos = length(xbins{nnn});     % total # of spatial locations (traj-dist combination)
        
        % initialize output
        posteriors{nnn} = nan(numlinpos,numbins);
        
        for bb = 1:numbins
            
            % inds of spikes in this bin
            inds = ( S{nnn}(:,2) == bb ) ;
            
            % Identify active ("a") encoding cells in this bin
            cellnums_a = unique( S{nnn}(inds,3) )' ;
                numactive = length(cellnums_a);
                
            % If no clustered spiking in this bin, then posterior set to 0s
            if isempty(cellnums_a)
                posteriors{nnn}(:,bb) = zeros(numlinpos,1);   % just report a series of 0s if there are no active spikes this bin
                continue
            end
            
            % Iterate through each cell, calculating their individual likelihood functions
            L = nan(numactive,numlinpos); 
            for vv = 1:numactive
                ccc = cellnums_a(vv);  % cell number in detc list
                bin_numspk = spikecounts{nnn}(ccc,bb);
                tet = detc{nnn}(ccc,3);
                cellno = detc{nnn}(ccc,4);
                % iterate over spatial locations
                for pp = 1:numlinpos
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
        
        % elminate invalid position bins %%%%%%%%%
        posteriors{nnn} = posteriors{nnn}(validxbins{nnn},:);
            numcutlinpos(nnn) = size(posteriors{nnn},1);
    
    end
    
        % save output
        cd(savedir)
        disp('done calculating full epoch posteriors')
        savefilename = sprintf('%s_%d_%d',animalname(1:3),day,ep);
        D = struct;
        D.animalname = animalname;
        D.dayep = [day ep];
        D.remoteep = remoteep;
        D.detc = detc;
        D.binvec_c = binvec_c;
        D.armdists = armdists_cat{1};   % actual animal positions
        D.armvect = armvect_cat{1};  % actual animal position-dir, offset by 1000 to distinguish INBOUND (toward center junction, < 1000) vs. OUTBOUND (away, > 1000)
        D.linposbins{1} = (1:numcutlinpos(1)) * xdel ;
        D.lastcenterbin = lastcenterbin;
        D.lastrightbin = lastrightbin;
        D.lastbin = numbinsingle;
        D.activebins = activebins;
        D.posteriors = posteriors;
        
        if remoteW && ~isnan(numcutlinpos(2)) > 1
            D.linposbins{2} = (1:numcutlinpos(2)) * xdel ;
        end
        
        save(savefilename,'D')
        
        clear D
            
    
    
end


