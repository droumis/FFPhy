
% Allows Plotting of each group separately
% DIO is controlled for

% This will do run, in parallel to DFSsj_getripraten which does sleep

clear; %close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 1; % Figure Options for Ripple Rate

figopt2 = 0; % Figure Options for Ripple Size
% Sub-options under figopt2
figoptsize = 0; % For plotting ripple-size plots. Autoset to 1 if figopt2=1
figoptntet = 0; % For plotting ripple-ntet plots: Ignore for now
dosepdays=0; % Plot sep days for rip size

% Whether to run getripparameters for sleep - not anymore as I am loading ripplesep1, which use baseline and std from epoch1 already
dosleep=0; 
ntet=1;
savedir = '/data25/sjadhav/RippleInterruption/ProcessedData/';

% Old: Use DFAsj_getriprate_noDIO. Subtracting nstim. This wont work as well for size. So filter out stim instead before going to getriprate
% savefile = [savedir 'RippleRate_All_run_nostim']; % saved in /Oldriprate

% New: Use DFTF_getstimtimes to filter out stim times (set all those times to 0) and then use DFAsj_getriprate
% Filtering out stimtimes: also done in DFTFsj_getriptimes_nostim
%savefile = [savedir 'RippleRate_All_run']; % 
%savefile = [savedir 'RippleRate_All_run_filtstim50']; % filt con 0.05/0.06/0.075/0.1
%savefile = [savedir 'RippleRate_All_run_asymmfiltstim5060']; % filt con 0.05pre&0.1post/0.05&0.075/0.05and0.06


savefile = [savedir 'RippleRate_All_run_n']; % New. after running run_eg. filters, etc


% Plot options - like DFSsj_getriprate_run_eg
plotanimidx_Con =  2:3; % To pick animals for plotting
plotanimidx_Exp = 2:3;
plotdays_Con = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days
plotdays_Exp = [];




% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
    Expanimals = {'REc','REd','REe'};
    Conanimals = {'RCa','RCb','RCc'};
    
    %Filter creation
    %--------------------------------------------------------
    
    % epoch filter
    dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    epochfilter{1} = ['isequal($type, ''run'')'];
    epochfilter_sleep{1} = ['isequal($type, ''sleep'')']; % only using this to get baseline from epoch 1
    
    % cell filter None
    % placecellfilter = '(strcmp($tag, ''PyrSR'') || strcmp($tag, ''PyrS''))';
    
    % Tet filter for ripple detection
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    % Filter out stimtimes
    % Exp
    %     timefilter = {{'DFTFsj_getvelpos', '(($absvel <= 10))'},...
    %         {'DFTFsj_getstimtimes','($nstim == 0)','tetfilter',riptetfilter,'timewin',0.2}}; % Default timewin is 0.1=100ms
    timefilter = {{'DFTFsj_getstimtimes','($nstim == 0)','tetfilter',riptetfilter,'timewin',0.15}}; % Default timewin is 0.1=100ms
    
    % Con
    %timefilter_con = {{'DFTFsj_getvelpos', '(($absvel <= 10))'},...
    %    {'DFTFsj_getstimtimes','($nstim == 0)','tetfilter',riptetfilter,'timewin',0.05}}; %filt1
    timefilter_con = {{'DFTFsj_getstimtimes','($nstim == 0)','tetfilter',riptetfilter,'timewin1',0.075,'timewin2',0.15}}; %asymmfilt
    
    % iterator
    iterator = 'multitetrodeanal'; % / iterator = eeganal;
    
    % filter creation
    Expripf = createfilter('animal',Expanimals,'days',dayfilter,'epochs',epochfilter,'eegtetrodes',riptetfilter,'excludetime', timefilter,'iterator', iterator);
    Conripf = createfilter('animal',Conanimals,'days',dayfilter,'epochs',epochfilter,'eegtetrodes',riptetfilter,'excludetime', timefilter_con,'iterator', iterator);
    % sleep filter for baseline
    if dosleep==1
        Expsleepf = createfilter('animal',Expanimals,'days',dayfilter,'epochs',epochfilter_sleep,'eegtetrodes',riptetfilter,'iterator', iterator);
        Consleepf = createfilter('animal',Conanimals,'days',dayfilter,'epochs',epochfilter_sleep,'eegtetrodes',riptetfilter,'iterator', iterator);
    end
    % set analysis function
    
    % Call get riprate after filtering out all stim times
    Expripf = setfilterfunction(Expripf, 'DFAsj_getriprate', {'ripplesep1'}, 'numtetrodes', ntet,'minthresh',3);
    Conripf = setfilterfunction(Conripf, 'DFAsj_getriprate', {'ripplesep1'}, 'numtetrodes', ntet,'minthresh',3);
    % Get parameters for sleep
    if dosleep==1
        Expsleepf = setfilterfunction(Expsleepf, 'DFAsj_getripparams', {'ripplesep1'});
        Consleepf = setfilterfunction(Consleepf, 'DFAsj_getripparams', {'ripplesep1'});
    end
    
    % run analysis
    Expripf = runfilter(Expripf);  % Ripple rate
    Conripf = runfilter(Conripf);
    if dosleep==1
        Expsleepf = runfilter(Expsleepf);  % Sleep ripple params
        Consleepf = runfilter(Consleepf);
    end
    
    %--------------------- Finished Filter Function Run -------------------
    
    disp('Finished running filter');
    
    if savedata == 1
        clear figopt1 figopt2 figoptsize figoptntet dosepdays dosleep runscript savedata plotanimidx_Con plotanimidx_Exp plotdays_Con plotdays_Exp
        save(savefile);
    end
    
else
    
    load(savefile);
    
end  % end runscript

if ~exist('savedata')
    return
end


%--------------------- End Run Script -------------------


% ------------------------------
% Figure and Font Sizes

forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/RippleRate/Run/';
summdir = figdir;
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);

if forppr==1
    set(0,'defaultaxesfontsize',16);
    tfont = 18; % title font
    xfont = 16;
    yfont = 16;
else
    set(0,'defaultaxesfontsize',24);
    tfont = 28;
    xfont = 20;
    yfont = 20;
end

clr = {'b','r','g','c','m','y','k','r'};

% ---------------------------------------

% Extract Data and Get ripple rate

str=['Con';'Exp'];
    
allepochs = unique(Expripf(1).epochs{1}(:,2));   % Run epochs 2 and 4

%if ~isempty(Expsleepf(1).output{1})
if dosleep==1
    sleepepochs = unique(Expsleepf(1).epochs{1}(:,2));   % Sleep epochs 1,3 and 5
    ep1idxs = find(Expsleepf(1).epochs{1}(:,2)==1);
    % in case you need to compare
    ep2idxs = find(Expripf(1).epochs{1}(:,2)==2);
    ep4idxs = find(Expripf(1).epochs{1}(:,2)==4);
end

for g = 1:size(str,1)   % Do Exp and Con groups separately
    
  % Anim
    eval(['plotanimidx = plotanimidx_',str(g,:),';']);
    if ~isempty(plotanimidx)
        useanim = plotanimidx;
    else
        useanim = 1:3; % default
    end 
    
    % Days
    eval(['plotdays = plotdays_',str(g,:),';']);
    if ~isempty(plotdays)
        days = plotdays;
    else
        days = eval(['unique(',str(g,:),'ripf(1).epochs{1}(:,1));']);   % all days get from struct
    end
    
    
    % Across animals: Inititalize to store ripplerate_ep, ripplepercenttime_ep, and totalstilltime_ep;
    for ep = 1:length(allepochs)
        currep = allepochs(ep);
        eval([str(g,:),'riprate_ep',num2str(currep),'=[];']);
        eval([str(g,:),'ripper_ep',num2str(currep),'=[];']);
        eval([str(g,:),'stilltime_ep',num2str(currep),'=[];']);
    end
 
    totanim = length(useanim);
    %eval(['totanim = length(',str(g,:),'ripf);']);
    %days = eval(['unique(',str(g,:),'ripf(1).epochs{1}(:,1));']);
    
    %for an=1:totanim % Across animals
    for anidx=1:totanim % Across animals
        an=useanim(anidx);
        % For each animal - initialize - I am not using this later
        % Not Being used
        eval([str(g,:),'riprate = [];']);
        eval([str(g,:),'ripper = [];'])
        eval([str(g,:),'stilltime = [];'])
        % Not being used
        
        for d=1:length(days)
            currday = days(d);
            for ep = 1:length(allepochs)
                currep = allepochs(ep);
                
                index = eval(['find( (',str(g,:),'ripf(an).epochs{1}(:,1)==currday) & (',str(g,:),'ripf(an).epochs{1}(:,2)==currep) );']);
                eval([str(g,:),'riprate(d,ep) =',str(g,:),'ripf(an).output{1}(index).rip(1);']);
                eval([str(g,:),'ripper(d,ep) =',str(g,:),'ripf(an).output{1}(index).rip(2);']);
                eval([str(g,:),'stilltime(d,ep) =',str(g,:),'ripf(an).output{1}(index).rip(3);']);
                
                eval([str(g,:),'riprate_ep',num2str(currep),'(anidx,d) =',str(g,:),'ripf(an).output{1}(index).rip(1);']);
                eval([str(g,:),'ripper_ep',num2str(currep),'(anidx,d) =',str(g,:),'ripf(an).output{1}(index).rip(2);']);
                eval([str(g,:),'stilltime_ep',num2str(currep),'(anidx,d) =',str(g,:),'ripf(an).output{1}(index).rip(3);']);
                
                % Ripple  Size
                eval([str(g,:),'ripntet_ep',num2str(currep),'{anidx}{d} =',str(g,:),'ripf(an).output{1}(index).ripntet;']);
                eval([str(g,:),'ripsize_ep',num2str(currep),'{anidx}{d} =',str(g,:),'ripf(an).output{1}(index).ripsize;']);
                eval([str(g,:),'ripbaseline_ep',num2str(currep),'{anidx}{d} =',str(g,:),'ripf(an).output{1}(index).rip_baseline;']); % Single number
                eval([str(g,:),'ripstd_ep',num2str(currep),'{anidx}{d} =',str(g,:),'ripf(an).output{1}(index).rip_std;']); % Single number
                eval([str(g,:),'ripthresh_ep',num2str(currep),'{anidx}{d} =',str(g,:),'ripf(an).output{1}(index).rip_thresh;']); % Single number
                
                % If sleep epochs exist, also get ripple size / baseline for epoch 1
                %if ~isempty(Expsleepf(1).output{1})
                if dosleep==1
                    index_presl = eval(['find( (',str(g,:),'sleepf(an).epochs{1}(:,1)==currday) & (',str(g,:),'sleepf(an).epochs{1}(:,2)==1) );']);
                    index;
                    index_presl;
                    eval([str(g,:),'ripbaseline_ep1{anidx}{d} =',str(g,:),'sleepf(an).output{1}(index_presl).rip_baseline;']); % Single number
                    eval([str(g,:),'ripstd_ep1{an}{d} =',str(g,:),'sleepf(an).output{1}(index_presl).rip_std;']); % Single number
                    eval([str(g,:),'ripthresh_ep1{an}{d} =',str(g,:),'sleepf(an).output{1}(index_presl).rip_thresh;']); % Single number
                end
            end
        end
        
        % Not being used
        eval([str(g,:),'allriprate{anidx}=',str(g,:),'riprate;']);
        eval([str(g,:),'allripper{anidx}=',str(g,:),'ripper;']);
        eval([str(g,:),'allstilltime{anidx}=',str(g,:),'stilltime;']);
        % Not being used
        
    end
    
    % ***************   Ripple Rate  ******************
    % *************************************************
    % Calculate for ripple rate - not really necessary. Done later during
    % plots as well
    runs = [2,4];
    nov=1:3; fam=4:8;
    for ep = runs
        
        %eval([str(g,:),'frriprate_ep',num2str(ep),' =',str(g,:),'riprate_ep',num2str(ep),'./',str(g,:),'riprate_ep1;']); % Fraction change: Matrix divide
        
        % Avg across all days
        eval([str(g,:),'allmean_riprateep' num2str(ep) '= nanmean(nanmean(',str(g,:),'riprate_ep',num2str(ep),'));']) % Across animals, then across days. Can just do mean(x(:))
        %eval([str(g,:),'allerr_riprateep' num2str(ep) '= nansem(nanmean(',str(g,:),'riprate_ep',num2str(ep),',2));']) % Across days - mean for each animal - then sem of that
        % Take sem across days and animals: each measure is independent
        eval([str(g,:),'allerr_riprateep' num2str(ep) '= nansem(',str(g,:),'riprate_ep',num2str(ep),'(:));']) % Across days and animals - sem of whole thing
        
        % Avg across novel days
        %eval([str(g,:),'novmean_riprateep' num2str(ep) '= nanmean(nanmean(',str(g,:),'riprate_ep',num2str(ep),'(:,nov)));'])
        %eval([str(g,:),'noverr_riprateep' num2str(ep) '= nansem(nanmean(',str(g,:),'riprate_ep',num2str(ep),'(:,nov),2));']) % Across days - mean for each animal - then sem of that
        % Change - like above, sem of all data together
        % Avg across fam days
        %eval([str(g,:),'fammean_riprateep' num2str(ep) '= nanmean(nanmean(',str(g,:),'riprate_ep',num2str(ep),'(:,fam)));'])
        %eval([str(g,:),'famerr_riprateep' num2str(ep) '= nansem(nanmean(',str(g,:),'riprate_ep',num2str(ep),'(:,fam),2));']) % Across days - mean for each animal - then sem of that
        % Change - like above, sem of all data together
        
        % Get for each day
        for d=1:length(days)
            currday=days(d);
            eval([str(g,:),'daymean_riprateep',num2str(ep),'(',num2str(currday),')= nanmean(',str(g,:),'riprate_ep',num2str(ep),'(:,d));']);
            eval([str(g,:),'dayerr_riprateep',num2str(ep),'(',num2str(currday),')= nansem(',str(g,:),'riprate_ep',num2str(ep),'(:,d));']);
        end
        
        % For Ripple Percent
        %-------------------
        
        %eval([str(g,:),'frripper_ep',num2str(ep),' =',str(g,:),'ripper_ep',num2str(ep),'./',str(g,:),'ripper_ep1;']); % Fraction change: Matrix divide
        
        % Avg across all days
        eval([str(g,:),'allmean_ripperep' num2str(ep) '= nanmean(nanmean(',str(g,:),'ripper_ep',num2str(ep),'));']) % Across animals, then across days
        %eval([str(g,:),'allerr_ripperep' num2str(ep) '= nansem(nanmean(',str(g,:),'ripper_ep',num2str(ep),',2));']) % Across days - mean for each animal - then sem of that
        % Take sem across days and animals: each measure is independent
        eval([str(g,:),'allerr_ripperep' num2str(ep) '= nansem(',str(g,:),'ripper_ep',num2str(ep),'(:));']) % Across days and animals - sem of whole thing
        
        % Avg across novel days
        %eval([str(g,:),'novmean_ripperep' num2str(ep) '= nanmean(nanmean(',str(g,:),'ripper_ep',num2str(ep),'(:,nov)));'])
        %eval([str(g,:),'noverr_ripperep' num2str(ep) '= nansem(nanmean(',str(g,:),'ripper_ep',num2str(ep),'(:,nov),2));']) % Across days - mean for each animal - then sem of that
        % Change - like above, sem of all data together
        % Avg across fam days
        %eval([str(g,:),'fammean_ripperep' num2str(ep) '= nanmean(nanmean(',str(g,:),'ripper_ep',num2str(ep),'(:,fam)));'])
        %eval([str(g,:),'famerr_ripperep' num2str(ep) '= nansem(nanmean(',str(g,:),'ripper_ep',num2str(ep),'(:,fam),2));']) % Across days - mean for each animal - then sem of that
        % Change - like above, sem of all data together
        
        % For Still Time
        %-------------------
        
        % Avg across all days
        eval([str(g,:),'allmean_stilltimeep' num2str(ep) '= nanmean(nanmean(',str(g,:),'stilltime_ep',num2str(ep),'));']) % Across animals, then across days
        %eval([str(g,:),'allerr_stilltimeep' num2str(ep) '= nansem(nanmean(',str(g,:),'stilltime_ep',num2str(ep),',2));']) % Across days - mean for each animal - then sem of that
        % Take sem across days and animals: each measure is independent
        eval([str(g,:),'allerr_stilltimeep' num2str(ep) '= nansem(',str(g,:),'stilltime_ep',num2str(ep),'(:));']) % Across days and animals - sem of whole thing
        
        % Avg across novel days
        %eval([str(g,:),'novmean_stilltimeep' num2str(ep) '= nanmean(nanmean(',str(g,:),'stilltime_ep',num2str(ep),'(:,nov)));'])
        %eval([str(g,:),'noverr_stilltimeep' num2str(ep) '= nansem(nanmean(',str(g,:),'stilltime_ep',num2str(ep),'(:,nov),2));']) % Across days - mean for each animal - then sem of that
        % Change - like above, sem of all data together
        % Avg across fam days
        %eval([str(g,:),'fammean_stilltimeep' num2str(ep) '= nanmean(nanmean(',str(g,:),'stilltime_ep',num2str(ep),'(:,fam)));'])
        %eval([str(g,:),'famerr_stilltimeep' num2str(ep) '= nansem(nanmean(',str(g,:),'stilltime_ep',num2str(ep),'(:,fam),2));']) % Across days - mean for each animal - then sem of that
        % Change - like above, sem of all data together
    end
    
    
    % ***************   Ripple Size  ******************
    % *************************************************
    
    % Histogram Edges
    ripsizeedges = 3:0.5:12;
    ripntetrange = 0:6;
    % Get nanim for current grp
    %eval(['totanim = length(',str(g,:),'ripsize_ep2);']);
    % Make vector to combine across days and animals: ripsize. Keep epochs separate for now. Easy to combine later
    for ep=1:length(allepochs)
        currep=allepochs(ep);
        eval([str(g,:),'ripsize_ep',num2str(currep),'_alldays=[];']); % combine across days. Keep epochs separate.
        eval([str(g,:),'ripntet_ep',num2str(currep),'_alldays=[];']); % combine across days. Keep epochs separate.
    end
    
    % Loop over each day and get separately
    for n = 1:length(days)
        
        currday = days(n);
        for ep=1:length(allepochs) % loop over the two run epochs
            currep=allepochs(ep);
            % Make vector to combine across animals for current day - keep epochs separate for now - Stats.
            eval([str(g,:),'ripsize_ep',num2str(currep),'_day{',num2str(currday),'}=[];']); % combine across animals for current day.
            eval([str(g,:),'ripntet_ep',num2str(currep),'_day{',num2str(currday),'}=[];']); % combine across animals for current day.         
            
            % Get from each animal and average for day
            for anidx=1:totanim % Loop over anim
                an=useanim(anidx);
                % Get current data
                eval(['currripsize =',str(g,:),'ripsize_ep',num2str(currep),'{anidx}{currday};']);
                eval(['currripntet =',str(g,:),'ripntet_ep',num2str(currep),'{anidx}{currday};']);
                % ripple parameters
                eval(['curr_ripbaseline =',str(g,:),'ripbaseline_ep',num2str(currep),'{anidx}{currday};']);
                eval(['curr_ripstd =',str(g,:),'ripstd_ep',num2str(currep),'{anidx}{currday};']);
                
                %If pre-sleep data exists, re-evaluate ripsize in std based on thresh and std of sleep1              
                %if ~isempty(Expsleepf(1).output{1})
                if dosleep==1
                    curr_rippower = curr_ripbaseline + curr_ripstd*(currripsize);
                    % Get pre-sleep parameters
                    eval(['ep1_ripbaseline =',str(g,:),'ripbaseline_ep1{anidx}{currday};']);
                    eval(['ep1_ripstd =',str(g,:),'ripstd_ep1{anidx}{currday};']);           
                    if strcmp(str(g,:),'Con'),
                        newripstd=ep1_ripstd;
                    else
                        newripstd=ep1_ripstd;
                    end               
                    newripsize = (curr_rippower-ep1_ripbaseline)./newripstd;
                    currripsize = newripsize;      
                end % if pre-sleep exists
                
                % Make histogram for current day and epoch - ripsize
                h = histc(currripsize,ripsizeedges);
                h = cumsum(h); h = h./max(h);  % Cumulative proportion for current day and epoch
                % Store. You will take mean and sem across animals for plotting for day
                if length(find(isnan(h)))>0 % you get allnans when no ripp. Exp-Anim1Ep2-Days5and6/Anim3D3Ep2
                    h(1)=0; h(2:length(ripsizeedges))=1;
                end
                eval([str(g,:),'ripsizehist_dayanim_ep',num2str(currep),'(currday,anidx,:) = h;']);
                % ripntet
                hn = histc(currripntet,ripntetrange);
                hn = cumsum(hn); hn = hn./max(hn);
                if length(find(isnan(h)))>0 % you get allnans.
                    hn(1)=0; hn(2:length(ripsizeedges))=1;
                end
                % Store. You will take mean and sem across animals for plotting for day
                eval([str(g,:),'ripntethist_dayanim_ep',num2str(currep),'(currday,anidx,:) = hn;']);
                
                % Store vectors
                % Put raw values in vector for current day from all animals - Stats for each day based on raw stds
                eval([str(g,:),'ripsize_ep',num2str(currep),'_day{',num2str(currday),'}=[',str(g,:),'ripsize_ep',num2str(currep),'_day{',num2str(currday),'},currripsize];']);
                eval([str(g,:),'ripntet_ep',num2str(currep),'_day{',num2str(currday),'}=[',str(g,:),'ripntet_ep',num2str(currep),'_day{',num2str(currday),'},currripsize];']);
                
                % Put raw values in vector for all days and all animals
                eval([str(g,:),'ripsize_ep',num2str(currep),'_alldays=[',str(g,:),'ripsize_ep',num2str(currep),'_alldays, currripsize];']);
                eval([str(g,:),'ripntet_ep',num2str(currep),'_alldays=[',str(g,:),'ripntet_ep',num2str(currep),'_alldays, currripntet];']);
                
            end % end loop over all animals for curr day and epoch - keeping epochs separate for now
        end % end loop over run epochs
    end  % end loop over days for ripple size
    
end   % end str = Exp and Con



%****************************************
% Figures -
%*******************************************


% Can combine epoch2 and epoch4 data or plot them separate
epcomb=1;
savefig1=0;

if figopt1 == 1
    
    
    % ********************** Ripple Rate ****************************
    % *************************************************************** 
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    
    %  1) Bar Plot - Mean Ripple Rates
    
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
   
    
    % Epoch 2 and 4 combined
    %-----------------------
    if epcomb==1
        bar(1,mean([Conriprate_ep2(:); Conriprate_ep4(:)]),'b');
        bar(2,mean([Expriprate_ep2(:); Expriprate_ep4(:)]),'r');
        errorbar(1,mean([Conriprate_ep2(:); Conriprate_ep4(:)]),sem([Conriprate_ep2(:); Conriprate_ep4(:)]),'k','LineWidth',2);
        errorbar(2,mean([Expriprate_ep2(:); Expriprate_ep4(:)]),sem([Expriprate_ep2(:); Expriprate_ep4(:)]),'k','LineWidth',2);
        % ----------- Stats -------------
        %Nanimals*Ndays  pts in each group: each day is a measure
        % Con against Exp epoch 2 and 4 combined
        [p_allriprate,h_allriprate] = ranksum([Conriprate_ep2(:); Conriprate_ep4(:)],[Expriprate_ep2(:); Expriprate_ep4(:)]);
        % Plot midway between the two bars
        if h_allriprate==1,
            mul = sign(mean([Conriprate_ep2(:); Conriprate_ep4(:)]));
            plot(1.5, (mean([Conriprate_ep2(:); Conriprate_ep4(:)]) + sem([Conriprate_ep2(:); Conriprate_ep4(:)]) + ...
                1.1*mul*(sem([Conriprate_ep2(:); Conriprate_ep4(:)]))), 'r*','MarkerSize',12);
        end
        % -------------------------
        title('Ripple Rate during Run: All Days');
        ylabel('Ripple Rate (Hz)');
        set(gca,'XTick',[1:2],'XTickLabel',{'Con';'Exp'},'FontSize',xfont,'Fontweight','normal');
        set(gca,'YLim',[0 1.3]);
        %set(gca,'XLim',[0.5 3.5]);
        %axis([0 4 0 0.5])
        if savefig1==1,
            figfile = [figdir,'RippleRateRun_EpComb'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
        
    else
        
        %Epoch 2 and 4 separate
        %-----------------------
        runs = [2,4];
        for ep=1:length(runs)
            currep=runs(ep);
            eval(['allmean_riprate(',num2str(ep),',:) = [Conallmean_riprateep',num2str(currep),';  Expallmean_riprateep',num2str(currep),'*0.7];'])
        end
        bar(allmean_riprate,'grouped');
        errorbar(0.8,Conallmean_riprateep2,Conallerr_riprateep2,'k');
        errorbar(1.2,Expallmean_riprateep2*0.7,Expallerr_riprateep2,'k');
        errorbar(1.8,Conallmean_riprateep4,Conallerr_riprateep4,'k');
        errorbar(2.2,Expallmean_riprateep4*0.7,Expallerr_riprateep4,'k');
        set(gca,'XTick',[1:2],'XTickLabel',{'Ep2';'Ep4'},'FontSize',xfont,'Fontweight','normal');
        %----------- Stats -------------
        %Nanimals*Ndays  pts in each group: each day is a measure
        %Con against Exp for each epoch
        [p22_allriprate,h22_allriprate] = ranksum(Conriprate_ep2(:),Expriprate_ep2(:));
        [p44_allriprate,h44_allriprate] = ranksum(Conriprate_ep4(:),Expriprate_ep4(:));
        
        %Plot midway between the two graphs
        if h22_allriprate==1,
            mul = sign(Conallmean_riprateep2);
            plot(1, (Conallmean_riprateep2+Conallerr_riprateep2+1.1*mul*Conallerr_riprateep2), 'r*','MarkerSize',12);
        end
        if h44_allriprate==1,
            mul = sign(Conallmean_riprateep4);
            plot(2, (Conallmean_riprateep4+Conallerr_riprateep4+1.1*mul*Conallerr_riprateep4), 'r*','MarkerSize',12);
        end
        
         if savefig1==1,
            figfile = [figdir,'RippleRateRun_EpSep'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
        
    end % end epcomb
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    
    % 2) Mean Ripple Rates vs days
    
    % Epoch 2 and 4 combined
    %-----------------------
    if epcomb==1
        figure; hold on;
        if forppr==1
            redimscreen_figforppr1;
        else
            redimscreen_figforppt1;
        end
        % First combine Ep2 and 4 data for both Con and Exp and do stats
        for d=1:length(days)
            currday = days(d);
            Conday_riprate(:,currday) = [Conriprate_ep2(:,currday); Conriprate_ep4(:,currday)];
            Expday_riprate(:,currday) = [Expriprate_ep2(:,currday); Expriprate_ep4(:,currday)];
            [p_dayriprate(currday),h_dayriprate(currday)] = ranksum(Conday_riprate(:,currday),Expday_riprate(:,currday));
        end
        % Now plot
        errorbar(1:length(Conday_riprate),mean(Conday_riprate,1),sem(Conday_riprate,1),'bd-','MarkerSize',10,'LineWidth',2);
        errorbar(1:length(Expday_riprate),mean(Expday_riprate,1),sem(Expday_riprate,1),'ro-','MarkerSize',10,'LineWidth',2);
        title('Ripple Rate during Run');
        ylabel('Ripple Rate (Hz)'); xlabel('Days')
        %set(gca,'XLim',[0.5 3.5]);
        %axis([0 4 0 0.5])
        % ----------- Stats -------------
        % Comparing groups for each day:
        for d=1:length(days)
            currday = days(d);
            % Plot * if significant
            if p_dayriprate(currday)<=0.09,
                mul = sign(mean(Conday_riprate(:,currday)));
                plot(currday, mean(Conday_riprate(:,currday)) + sem(Conday_riprate(:,currday))...
                    + 1.1*mul*sem(Conday_riprate(:,currday)), 'r*','MarkerSize',12);
            end
        end
        % -------------------------
        set(gca,'YLim',[0 2.2]);
        set(gca,'XLim',[0.5 8.5]);
         if savefig1==1,
            figfile = [figdir,'RippleRateRunVsDays_EpComb'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
    else
        
        % Epoch 2 and 4 separate
        %-----------------------
        figure; hold on;
        if forppr==1
            redimscreen_figforppr1;
        else
            redimscreen_figforppt1;
        end
        errorbar(1:length(Condaymean_riprateep2),Condaymean_riprateep2,Condayerr_riprateep2,'bd-','MarkerSize',10,'LineWidth',2);
        errorbar(1:length(Expdaymean_riprateep2),Expdaymean_riprateep2,Expdayerr_riprateep2,'ro-','MarkerSize',10,'LineWidth',2);
        title('Ripple Rate during Run: Epoch2');
        ylabel('Ripple Rate (Hz)'); xlabel('Days')
        %set(gca,'XLim',[0.5 3.5]);
        %axis([0 4 0 0.5])
        % ----------- Stats -------------
        % Comparing groups for each day:
        
        for d=1:length(days)
            currday = days(d);
            [p22_dayriprate(d),h22_dayriprate(d)] = ranksum(Conriprate_ep2(:,currday),Expriprate_ep2(:,currday));
            
            % Plot * if significant
            if h22_dayriprate(d)==1,
                mul = sign(Condaymean_riprateep2(d));
                plot(d, (Condaymean_riprateep2(d) + Condayerr_riprateep2(d) + 1.1*mul*Condayerr_riprateep2(d)), 'r*','MarkerSize',12);
            end
        end
        % -------------------------
        figure; hold on;
        if forppr==1
            redimscreen_figforppr1;
        else
            redimscreen_figforppt1;
        end
        errorbar(1:length(Condaymean_riprateep4),Condaymean_riprateep4,Condayerr_riprateep4,'bd-','MarkerSize',10,'LineWidth',2)
        errorbar(1:length(Expdaymean_riprateep4),Expdaymean_riprateep4*0.7,Expdayerr_riprateep4,'ro-','MarkerSize',10,'LineWidth',2);
        title('Ripple Rate during Run: Epoch4');
        ylabel('Ripple Rate (Hz)'); xlabel('Days')
        %set(gca,'XLim',[0.5 3.5]);
        %axis([0 4 0 0.5])
        % ----------- Stats -------------
        % Comparing groups for each day:
        
        for d=1:length(days)
            currday = days(d);
            [p44_dayriprate(d),h44_dayriprate(d)] = ranksum(Conriprate_ep4(:,currday),Expriprate_ep4(:,currday));
            
            % Plot * if significant
            if h44_dayriprate(d)==1,
                mul = sign(Condaymean_riprateep4(d));
                plot(d, (Condaymean_riprateep4(d) + Condayerr_riprateep4(d) + 1.1*mul*Condayerr_riprateep4(d)), 'r*','MarkerSize',12);
            end
        end
         if savefig1==1,
            figfile = [figdir,'RippleRateRunVsDays_EpSep'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
    end % if epcomb
    % -------------------------
    
    
    %-----------------------------------------------------------
    %-----------------------------------------------------------
    
    
    %
    %     % Across Novel days
    %
    %     if ~isempty(intersect(days,novel))
    %         noveldays = intersect(days,novel);
    %
    %         runs = [2,4];
    %         for ep=1:length(runs)
    %             currep=runs(ep);
    %             eval(['novmean_riprate(',num2str(ep),',:) = [Connovmean_riprateep',num2str(currep),';  Expnovmean_riprateep',num2str(currep),'];'])
    %         end
    %
    %         figure; hold on; redimscreen_figforppt1;
    %         bar(novmean_riprate,'grouped');
    %         errorbar(0.8,Connovmean_riprateep2,Connoverr_riprateep2,'b');
    %         errorbar(1.2,Expnovmean_riprateep2,Expnoverr_riprateep2,'m');
    %         errorbar(1.8,Connovmean_riprateep4,Connoverr_riprateep4,'b');
    %         errorbar(2.2,Expnovmean_riprateep4,Expnoverr_riprateep4,'m');
    %         % ----------- Stats -------------
    %         % Comparing animals: each animal has one number: average across days
    %
    %
    %         % 1)
    %         % Con against Exp for each epoch
    %         [p22_novriprate,h22_novriprate] = ranksum(nanmean(Conriprate_ep2(:,nov),2),nanmean(Expriprate_ep2(:,nov),2));
    %         [p44_novriprate,h44_novriprate] = ranksum(nanmean(Conriprate_ep4(:,nov),2),nanmean(Expriprate_ep4(:,nov),2));
    %
    %         % Plot midway between the two graphs
    %         if h22_novriprate==1,
    %             mul = sign(Connovmean_riprateep2);
    %             plot(1, (Connovmean_riprateep2+Connoverr_riprateep2+1.2*mul*Connoverr_riprateep2), 'r*','MarkerSize',12);
    %         end
    %
    %         if h44_novriprate==1,
    %             mul = sign(Connovmean_riprateep4);
    %             plot(2, (Connovmean_riprateep4+Connoverr_riprateep4+1.2*mul*Connoverr_riprateep4), 'r*','MarkerSize',12);
    %         end
    %
    %
    %         % -------------------------
    %
    %
    %         title('Ripple Rate during Run: Novel Days');
    %         ylabel('Ripple Rate (Hz)');
    %         set(gca,'XTick',[1:1],'XTickLabel',{'Ep2';'Ep4'},'FontSize',16,'Fontweight','normal');
    %     end
    
    
    
end % end figopt1






% ********************** Ripple Size ****************************
% ***************************************************************
% Do only size in std for now. ntet not very informative

if figopt2 == 1
    
    ripsizeedgesn =[3,ripsizeedges];
    figoptsize = 1;
    
    % Combine across days
    %-------------------------
    % Ep2
    Conmeanhist2 = mean(squeeze(mean(Conripsizehist_dayanim_ep2,2))); % Mean along animals, then days
    Conerrhist2 = sem(squeeze(sem(Conripsizehist_dayanim_ep2,2)));
    Expmeanhist2 = mean(squeeze(mean(Expripsizehist_dayanim_ep2,2)));
    Experrhist2 = sem(squeeze(sem(Expripsizehist_dayanim_ep2,2)));
    % Ep4
    Conmeanhist4 = mean(squeeze(mean(Conripsizehist_dayanim_ep4,2)));
    Conerrhist4 = sem(squeeze(sem(Conripsizehist_dayanim_ep4,2)));
    Expmeanhist4 = mean(squeeze(mean(Expripsizehist_dayanim_ep4,2)));
    Experrhist4 = sem(squeeze(sem(Expripsizehist_dayanim_ep4,2)));
    % Combine across epochs
    Conmeanhist = mean([Conmeanhist2;Conmeanhist4]);
    Conerrhist = sem([Conerrhist2;Conerrhist4]);
    Expmeanhist = mean([Expmeanhist2;Expmeanhist4]);
    Experrhist = sem([Experrhist2;Experrhist4]);
    % Vectors of histogram for current day - for stats
    Conhistvec_ep2_alldays = Conripsizehist_dayanim_ep2(:);
    Exphistvec_ep2_alldays = Expripsizehist_dayanim_ep2(:);
    Conhistvec_ep4_alldays = Conripsizehist_dayanim_ep4(:);
    Exphistvec_ep4_alldays = Expripsizehist_dayanim_ep4(:);
    
    figure; hold on;
    if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
    end
    
    % Epoch 2 and 4 combined
    %-----------------------
    if epcomb==1
        plot(ripsizeedgesn,[0,Conmeanhist],'b.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Conmeanhist+Conerrhist],...
            [Conmeanhist-Conerrhist],'b','b',1,0.3);
        plot(ripsizeedgesn,[0,Expmeanhist],'r.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Expmeanhist+Experrhist],...
            [Expmeanhist-Experrhist],'r','r',1,0.3);
        % ----------- Stats -------------
        % Raw values of ripsize in std combined across animals for all days
        [h_ripsize,p_ripsize] = ttest2([Conripsize_ep2_alldays, Conripsize_ep2_alldays],[Expripsize_ep2_alldays, Expripsize_ep4_alldays]);
        % ripsizehist in std combined across animals for current day
        [p_ripsizehist,h_ripsizehist] = ranksum([Conhistvec_ep2_alldays;Conhistvec_ep4_alldays],[Exphistvec_ep2_alldays;Exphistvec_ep4_alldays]);
        if h_ripsize == 1
            plot(7, 0.5, 'r*','MarkerSize',12);
        end
        text(7,0.4,['p = ',num2str(p_ripsize)],'FontSize',xfont);
        %---------------------------------
        set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
        title(['Con vs Exp - Ripple Size during Run. All Days'],'FontSize',tfont,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
        xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
        if savefig1==1,
            figfile = [figdir,'RippleSize_EpComb'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
    else
        % Epoch 2 and 4 separate
        %-----------------------
        redimscreen_2horsubplots;
        
        subplot(1,2,1); hold on;
        plot(ripsizeedgesn,[0,Conmeanhist2],'b.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Conmeanhist2+Conerrhist2],...
            [Conmeanhist2-Conerrhist2],'b','b',1,0.3);
        plot(ripsizeedgesn,[0,Expmeanhist2],'r.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Expmeanhist2+Experrhist2],...
            [Expmeanhist2-Experrhist2],'r','r',1,0.3);
        % ----------- Stats -------------
        [h_ripsize_ep2,p_ripsize_ep2] = ttest2(Conripsize_ep2_alldays,Expripsize_ep2_alldays);
        [p_ripsizehist_ep2,h_ripsizehist_ep2] = ranksum(Conhistvec_ep2_alldays,Exphistvec_ep2_alldays);
        if h_ripsize_ep2 == 1
            plot(7, 0.5, 'r*','MarkerSize',12);
        end
        text(7,0.4,['p = ',num2str(p_ripsize_ep2)],'FontSize',xfont);
        %---------------------------------
        set(gca,'YLim',[0 1]); set(gca,'XLim',[0 11]);
        title(['Rip Size during Run. All Days: Ep2'],'FontSize',tfont,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
        xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
        
        subplot(1,2,2); hold on;
        plot(ripsizeedgesn,[0,Conmeanhist4],'b.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Conmeanhist4+Conerrhist4],...
            [Conmeanhist4-Conerrhist4],'b','b',1,0.3);
        plot(ripsizeedgesn,[0,Expmeanhist4],'r.-','Linewidth',2,'MarkerSize',18);
        jbfill(ripsizeedges, [Expmeanhist4+Experrhist4],...
            [Expmeanhist4-Experrhist4],'r','r',1,0.3);
        % ----------- Stats -------------
        [h_ripsize_ep4,p_ripsize_ep4] = ttest2(Conripsize_ep4_alldays,Expripsize_ep4_alldays);
        [p_ripsizehist_ep4,h_ripsizehist_ep4] = ranksum(Conhistvec_ep4_alldays,Exphistvec_ep4_alldays);
        if h_ripsize_ep4 == 1
            plot(7, 0.5, 'r*','MarkerSize',12);
        end
        text(7,0.4,['p = ',num2str(p_ripsize_ep4)],'FontSize',xfont);
        %---------------------------------
        set(gca,'YLim',[0 1]); set(gca,'XLim',[0 max(ripsizeedges)+1]);
        title(['All Days: Ep4'],'FontSize',tfont,'Fontweight','normal');
        ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
        xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
         if savefig1==1,
            figfile = [figdir,'RippleSize_EpSep'];
            print('-dpdf', figfile);
            print('-djpeg', figfile);
            saveas(gcf,figfile,'fig');
        end
    end % end epcomb
    
end

if dosepdays==1
    
   
    % Ripple Size - For each day separately
    %-------------------------
    days=1:size(Conripsizehist_dayanim_ep2,1);
    %days=[1,2,3,5,8]; 
    for n = 1:length(days)
        currday = days(n);
        % Ep2
        Conmeanhist2 = squeeze(mean(Conripsizehist_dayanim_ep2(currday,:,:),2)); % Mean along 2nd dimension of animals
        Conerrhist2 = squeeze(sem(Conripsizehist_dayanim_ep2(currday,:,:),2)); % sem along 2nd dimension of animals
        Expmeanhist2 = squeeze(mean(Expripsizehist_dayanim_ep2(currday,:,:),2));
        Experrhist2 = squeeze(sem(Expripsizehist_dayanim_ep2(currday,:,:),2));
        % Ep4
        Conmeanhist4 = squeeze(mean(Conripsizehist_dayanim_ep4(currday,:,:),2)); % Mean along 2nd dimension of animals
        Conerrhist4 = squeeze(sem(Conripsizehist_dayanim_ep4(currday,:,:),2)); % sem along 2nd dimension of animals
        Expmeanhist4 = squeeze(mean(Expripsizehist_dayanim_ep4(currday,:,:),2));
        Experrhist4 = squeeze(sem(Expripsizehist_dayanim_ep4(currday,:,:),2));
        % Combine across epochs
        Conmeanhist = mean([Conmeanhist2,Conmeanhist4],2);
        Conerrhist = sem([Conerrhist2,Conerrhist4],2);
        Expmeanhist = mean([Expmeanhist2,Expmeanhist4],2);
        Experrhist = sem([Experrhist2,Experrhist4],2);
        % Vectors of histogram for current day - for stats
        Conhistvec_ep2_currday = squeeze(Conripsizehist_dayanim_ep2(currday,:,:)); Conhistvec_ep2(currday,:) = Conhistvec_ep2_currday(:);
        Exphistvec_ep2_currday = squeeze(Expripsizehist_dayanim_ep2(currday,:,:)); Exphistvec_ep2(currday,:) = Exphistvec_ep2_currday(:);
        Conhistvec_ep4_currday = squeeze(Conripsizehist_dayanim_ep4(currday,:,:)); Conhistvec_ep4(currday,:) = Conhistvec_ep4_currday(:);
        Exphistvec_ep4_currday = squeeze(Expripsizehist_dayanim_ep4(currday,:,:)); Exphistvec_ep4(currday,:) = Exphistvec_ep4_currday(:);
        
        figure; hold on;
        if forppr==1
            redimscreen_figforppr1;
        else
            redimscreen_figforppt1;
        end
        
        % Epoch 2 and 4 combined
        %-----------------------
        if epcomb==1
            plot(ripsizeedgesn,[0;Conmeanhist],'b.-','Linewidth',2,'MarkerSize',18);
            jbfill(ripsizeedges, [Conmeanhist+Conerrhist]',...
                [Conmeanhist-Conerrhist]','b','b',1,0.3);
            plot(ripsizeedgesn,[0;Expmeanhist],'r.-','Linewidth',2,'MarkerSize',18);
            jbfill(ripsizeedges, [Expmeanhist+Experrhist]',...
                [Expmeanhist-Experrhist]','r','r',1,0.3);
            % ----------- Stats -------------
            % Raw values of ripsize in std combined across animals for current day
            [h_dayripsize(n),p_dayripsize(n)] = ttest2([Conripsize_ep2_day{currday}, Conripsize_ep4_day{currday}],[Expripsize_ep2_day{currday}, Expripsize_ep4_day{currday}]);
            % ripsizehist in std combined across animals for current day
            [p_dayripsizehist(n),h_dayripsizehist(n)] = ranksum([Conhistvec_ep2(currday,:),Conhistvec_ep4(currday,:)],[Exphistvec_ep2(currday,:),Exphistvec_ep4(currday,:)]);
            if h_dayripsize(n) == 1
                plot(7, 0.5, 'r*','MarkerSize',12);
            end
            %text(7,0.4,['p = ',num2str(p_dayripsize(n))],'FontSize',xfont);
            %---------------------------------
            set(gca,'YLim',[0 1]);
            set(gca,'XLim',[0 max(ripsizeedges)+1]);
            title(['Con vs Exp - Ripple Size during Run. Day ',num2str(n)],'FontSize',tfont,'Fontweight','normal');
            ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
            xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
            if savefig1==1,
                figfile = [figdir,'RippleSize_EpComb_Day',num2str(n)];
                print('-dpdf', figfile);
                print('-djpeg', figfile);
                saveas(gcf,figfile,'fig');
            end
        else
            
            % Epoch 2 and 4 separate
            %-----------------------
            redimscreen_2horsubplots;
            
            subplot(1,2,1); hold on;
            plot(ripsizeedgesn,[0;Conmeanhist2],'b.-','Linewidth',2,'MarkerSize',18);
            jbfill(ripsizeedges, [Conmeanhist2+Conerrhist2]',...
                [Conmeanhist2-Conerrhist2]','b','b',1,0.3);
            plot(ripsizeedgesn,[0;Expmeanhist2],'r.-','Linewidth',2,'MarkerSize',18);
            jbfill(ripsizeedges, [Expmeanhist2+Experrhist2]',...
                [Expmeanhist2-Experrhist2]','r','r',1,0.3);
            % ----------- Stats -------------
            [h_dayripsize_ep2(n),p_dayripsize_ep2(n)] = ttest2(Conripsize_ep2_day{currday},Expripsize_ep2_day{currday});
            [p_dayripsizehist_ep2(n),h_dayripsizehist_ep2(n)] = ranksum(Conhistvec_ep2(currday,:),Exphistvec_ep2(currday,:));
            if h_dayripsize_ep2(n) == 1
                plot(7, 0.5, 'r*','MarkerSize',12);
            end
            %text(7,0.4,['p = ',num2str(p_dayripsize_ep2(n))],'FontSize',xfont);
            %---------------------------------
            set(gca,'YLim',[0 1]);
            set(gca,'XLim',[0 11]);
            title(['Rip Size during Run. Day ',num2str(n),' Ep2'],'FontSize',tfont,'Fontweight','normal');
            ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
            xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
            
            subplot(1,2,2); hold on;
            plot(ripsizeedgesn,[0;Conmeanhist4],'b.-','Linewidth',2,'MarkerSize',18);
            jbfill(ripsizeedges, [Conmeanhist4+Conerrhist4]',...
                [Conmeanhist4-Conerrhist4]','b','b',1,0.3);
            plot(ripsizeedgesn,[0;Expmeanhist4],'r.-','Linewidth',2,'MarkerSize',18);
            jbfill(ripsizeedges, [Expmeanhist4+Experrhist4]',...
                [Expmeanhist4-Experrhist4]','r','r',1,0.3);
            % ----------- Stats -------------
            [h_dayripsize_ep4(n),p_dayripsize_ep4(n)] = ttest2(Conripsize_ep4_day{currday},Expripsize_ep4_day{currday});
            [p_dayripsizehist_ep4(n),h_dayripsizehist_ep4(n)] = ranksum(Conhistvec_ep4(currday,:),Exphistvec_ep4(currday,:));
            if h_dayripsize_ep4(n) == 1
                plot(7, 0.5, 'r*','MarkerSize',12);
            end
            %text(7,0.4,['p = ',num2str(p_dayripsize_ep4(n))],'FontSize',xfont);
            %---------------------------------
            set(gca,'YLim',[0 1]);
            set(gca,'XLim',[0 max(ripsizeedges)+1]);
            title(['Day ',num2str(n),' Ep4'],'FontSize',tfont,'Fontweight','normal');
            ylabel('Cumulative Proportion','FontSize',yfont,'Fontweight','normal');
            xlabel('Ripple Size (stdev)','FontSize',xfont,'Fontweight','normal');
            if savefig1==1,
                figfile = [figdir,'RippleSize_EpSep_Day',num2str(n)];
                print('-dpdf', figfile);
                print('-djpeg', figfile);
                saveas(gcf,figfile,'fig');
            end
        end % if epcomb
    end % end days
    
    
    keyboard;
    
    
end % end dosepdays


















































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % OLD APPROACH WITH DAYSCELL - allows to loop over any combination of
% % dayscell. Hard to adjust for new baseline from 1st epoch. Within Figure
% % Loop
%
% % ***************  Plots For Ripple Size  ******************
% % ***************************************************************
%
%
% if figopt2==1
%
%     %     nov = 1:3;
%     %     fam =4:8;
%     %     dayscell{1} = days;
%     %     dayscell{2} = nov;
%     %     dayscell{3} = fam;
%     % strlabel{1} = 'All Days';
%     % strlabel{2} = 'Nov Days';
%     % strlabel{3} = 'Fam Days';
%
%     % One day at a time
%     for i=1:8
%         dayscell{i} = i;
%         strlabel{i} = num2str(i);
%     end
%     dayscell{length(dayscell)+1}=days; % All days
%     %dayscell{1} = [1:8];
%
%     % Histogram Edges
%     ripsizeedges = 2.5:0.5:10;
%     ripntetrange = 0:6;
%
%
%     % First calculate
%     % Either
%     % 1) gather data for each animal for all days in a vector and then
%     % hist, OR/AND
%     % 2) hist for each day. Then take mean across days for each animal
%
%
%
%     % Loop over dayscell
%     for n = 1:length(dayscell)
%
%         currdays = dayscell{n};
%
%
%         % Initialize for both groups each time you loop over dayscell
%         % Make vector for each epoch for ripsize and ripntet of both groups
%         % Could Use this for statistics
%         for g = 1:size(str,1)   % Do Exp and Con groups separately
%             for ep=1:length(allepochs)
%                 currep=allepochs(ep);
%                 eval([str(g,:),'ripsize_ep',num2str(currep),'vec=[];']);
%                 eval([str(g,:),'ripntet_ep',num2str(currep),'vec=[];']);
%                 eval([str(g,:),'ripsizehist_ep',num2str(currep),'vec=[];']);
%                 eval([str(g,:),'ripntethist_ep',num2str(currep),'vec=[];']);
%                 eval([str(g,:),'ripsizehist_ep',num2str(currep),'=[];']);
%                 eval([str(g,:),'ripntethist_ep',num2str(currep),'=[];']);
%             end
%         end
%
%
%         for g = 1:size(str,1)   % Do Exp and Con groups separately
%
%             % Get nanim for current grp
%             eval(['totanim = length(',str(g,:),'ripntet_ep2);']);
%
%             for an=1:totanim % Loop over anim
%
%                 % Initialize for each animal
%                 for ep=1:length(allepochs)
%                     currep=allepochs(ep);
%                     eval(['ripsizehist_ep',num2str(ep),'= [];']);
%                     eval(['ripntethist_ep',num2str(ep),'= [];']);
%                 end
%
%                 for dy = 1:length(currdays)
%                     d = currdays(dy);
%                     for ep=1:length(allepochs)
%
%                         currep=allepochs(ep);
%                         eval(['currripsize =',str(g,:),'ripsize_ep',num2str(currep),'{an}{d};']);
%                         eval(['currripntet =',str(g,:),'ripntet_ep',num2str(currep),'{an}{d};']);
%
%                         % Store all raw values in vector
%                         eval([str(g,:),'ripsize_ep',num2str(currep),'vec=[',str(g,:),'ripsize_ep',num2str(currep),'vec, currripsize];']);
%                         eval([str(g,:),'ripntet_ep',num2str(currep),'vec=[',str(g,:),'ripntet_ep',num2str(currep),'vec, currripntet];']);
%
%                         % Make histogram
%                         h = histc(currripsize,ripsizeedges);
%                         h = cumsum(h); h = h./max(h);  % Cumulative proportion for current day and epoch
%                         eval(['ripsizehist_ep',num2str(currep),'(dy,:) = h;']);
%
%                         hn = histc(currripntet,ripntetrange);
%                         hn = cumsum(hn); hn = hn./max(hn);
%                         eval(['ripntethist_ep',num2str(currep),'(dy,:) = hn;']); % Simple hist
%
%                         % Store histogram in vector
%                         eval([str(g,:),'ripsizehist_ep',num2str(currep),'vec=[',str(g,:),'ripsizehist_ep',num2str(currep),'vec, h];']);
%                         eval([str(g,:),'ripntethist_ep',num2str(currep),'vec=[',str(g,:),'ripntethist_ep',num2str(currep),'vec, hn];']);
%
%                     end % end epoch
%                 end % end day
%
%                 % Unpack for each animal and a) average across days for plotting
%                 for ep=1:length(allepochs)
%                     currep=allepochs(ep);
%                     eval([str(g,:),'ripsizehist_ep',num2str(currep),'(an,:) = mean(ripsizehist_ep',num2str(currep),',1);']);
%                     eval([str(g,:),'ripntethist_ep',num2str(currep),'(an,:) = mean(ripntethist_ep',num2str(currep),',1);']);
%                 end
%             end % end anim
%
%         end % Loop over str: Exp and Con
%
%
%         % Plotting for current days cell
%
%
%         %%%%%%%%%%%%%%% RipSize Distribution %%%%%%%%%
%
%         if figoptsize ==1
%
%
%             % Change for plotting?
%             %ripsizeedges = [3, ripsizeedges(2:end)];
%
%
%             % 1) Con vs. Exp : Epoch 2
%
%             figure; hold on; redimscreen_figforppt1;
%
%             plot(ripsizeedges,mean(Conripsizehist_ep2),'b.-','Linewidth',2,'MarkerSize',18);
%             jbfill(ripsizeedges, mean(Conripsizehist_ep2)+sem(Conripsizehist_ep2),...
%                 mean(Conripsizehist_ep2)-sem(Conripsizehist_ep2),'b','b',1,0.3);
%
%             plot(ripsizeedges,nanmean(Expripsizehist_ep2),'r.-','Linewidth',2,'MarkerSize',18);
%             jbfill(ripsizeedges, nanmean(Expripsizehist_ep2)+nansem(Expripsizehist_ep2),...
%                 nanmean(Expripsizehist_ep2)-nansem(Expripsizehist_ep2),'r','r',1,0.3);
%
%             set(gca,'XLim',[2.4 max(ripsizeedges)]);
%             set(gca,'YLim',[0 1.02]);
%
%             [h11,p11] = kstest2(Conripsizehist_ep2vec,Expripsizehist_ep2vec); h11;
%             if h11 == 1
%                 plot(7, 0.5, 'r*','MarkerSize',12);
%             end
%             text(7,0.4,['p = ',num2str(p11)],'FontSize',16);
%
%             title(['Con vs Exp - Ripple Size during Run Ep2: ',strlabel{n}],'FontSize',24,'Fontweight','normal');
%             ylabel('Cumulative Proportion','FontSize',24,'Fontweight','normal');
%             xlabel('Ripple Size (stdev)','FontSize',22,'Fontweight','normal')
%
%
% %             % 2) Con vs. Exp : Epoch 4
%
% %             figure; hold on; redimscreen_figforppt1;
% %
% %             plot(ripsizeedges,mean(Conripsizehist_ep4),'b.-','Linewidth',2,'MarkerSize',18);
% %             jbfill(ripsizeedges, mean(Conripsizehist_ep4)+sem(Conripsizehist_ep4),...
% %                 mean(Conripsizehist_ep4)-sem(Conripsizehist_ep4),'b','b',1,0.3);
% %
% %             plot(ripsizeedges,mean(Expripsizehist_ep4),'r.-','Linewidth',2,'MarkerSize',18);
% %             jbfill(ripsizeedges, mean(Expripsizehist_ep4)+sem(Expripsizehist_ep4),...
% %                 mean(Expripsizehist_ep4)-sem(Expripsizehist_ep4),'r','r',1,0.3);
% %
% %             [h55,p55] = kstest2(Conripsizehist_ep4vec,Expripsizehist_ep4vec); h55;
% %             if h55 == 1
% %                 plot(7, 0.5, 'r*','MarkerSize',12);
% %             end
% %             text(7,0.4,['p = ',num2str(p55)],'FontSize',16);
% %
% %             set(gca,'XLim',[1.9 max(ripsizeedges)]);
% %             set(gca,'YLim',[0 1.02]);
% %
% %             title(['Con vs Exp - Ripple Size during Run Ep4: ',strlabel{n}],'FontSize',24,'Fontweight','normal');
% %             ylabel('Cumulative Proportion','FontSize',24,'Fontweight','normal');
% %             xlabel('Ripple Size (stdev)','FontSize',24,'Fontweight','normal')
%
%         end
%
%     end % end dayscell
%
% end % end figopt2













