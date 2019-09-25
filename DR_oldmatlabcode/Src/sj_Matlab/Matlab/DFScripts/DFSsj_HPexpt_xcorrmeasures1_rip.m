
% Also see DFSsj_placefields1 and mkarlsso/xcorrscript

% First use to get and plot overlap, trajdata and mapdata for run epochs: 
% Make a DFA_calcoverlap for this purpose
% Overlap: Dont Plot all at once for all days - too many plots!

% Then add calls to DFA_calcxcorrmeasures_linfields. Get xcorr measures for
% run and sleep epochs. Sleep filter can be speed

% Later on, can incorporate filters to include-only or exclude ripple times 


clear; %close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 1; % Figure Options - 

savedir = '/data25/sjadhav/HPExpt/ProcessedData/';

%savefile = [savedir 'HPa_xcorrmeasures_day2']; % CA1 vs PFC
%savefile = [savedir 'HPa_xcorrmeasures_day2_CA1only']; % 

%Using smoothing of 10 instead of 5 ms. bin =0.002; Usually, sw1=2.5*bin. Here, sw1=5*bin
savefile = [savedir 'HPa_xcorrmeasures_day2_CA1PFC_rip']; % 
%savefile = [savedir 'HPa_xcorrmeasures_day2_CA1only_rip']; % 

%savefile = [savedir 'HPa_xcorrmeasures_day8_CA1PFC_rip']; % 
%savefile = [savedir 'HPa_xcorrmeasures_day8_CA1only_rip']; % 


% Parameters
norm_overlap = 0;  % For place-field overlap
thresh_peakrate = 3; % For place-field overlap 

% Plot options
plotanimidx =  []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days



% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
        animals = {'HPa'};
        
    %Filter creation
    %-----------------------------------------------------
    
    % Epoch filter
    % -------------
    dayfilter = '2'; % Shantanu - I am adding day filter to parse out epoch filter
    runepochfilter = 'isequal($type, ''run'')';
    sleepepochfilter = 'isequal($type, ''sleep'')';
    
    % Cell filter
    % -----------
    cellpairfilter = {'allcomb','strcmp($tag, ''CA1Pyr'')|strcmp($tag, ''iCA1Pyr'')','strcmp($tag, ''PFC'')'};
    
    % Time filter
    % -----------
    
    % Using linfields for place field overlap - timefilter not being used 
    % Might be used for distance between peaks
    riptetfilter = '(isequal($descrip, ''riptet''))';
%     timefilter_place = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6},...
%         {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',2} }; % Dont care about stim here. Get rid of all ripples incl artifacts    

    % Get ripple times during run - with and without speed criterion Stim artifact removal implemented. Win around stimtimes is set to 0 in riptimes
    % Do as you do in DFSsj_getriprate_run
    % Exp
%     timefilterrun_onlyrip_speed = {{'DFTFsj_getvelpos', '(($absvel <= 10))'},...
%         {'DFTFsj_getriptimes','($nripples > 0)','tetfilter',riptetfilter,'minthresh',3}};

    timefilterrun_onlyrip = {{'DFTFsj_getriptimes','($nripples > 0)','tetfilter',riptetfilter,'minthresh',3}};
    
    % Get ripple times in sleep. Ripple on any tet of size > minthres. Stim times removed. No need for speed filter.
    timefiltersleep_onlyrip = {{'DFTFsj_getriptimes','($nripples > 0)','tetfilter',riptetfilter,'minthresh',3}};
    
    
    % Iterator
    % --------
    iterator = 'singlecellanal';
    
    % Filter creation
    % ----------------

    % Place field overlap - No timefilter. Trajs already are filtered for
    % no ripples and linvel>thresh
    xrun = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter,'cellpairs',...
        cellpairfilter,'iterator', iterator);
    
    % Sleep Corr - Timefilter: only during ripples
    sleepcorr_rip = createfilter('animal',animals,'days',dayfilter,'epochs',sleepepochfilter,'cellpairs',...
        cellpairfilter,'excludetime', timefiltersleep_onlyrip,'iterator', iterator);
    
    sleepcorr_all = createfilter('animal',animals,'days',dayfilter,'epochs',sleepepochfilter,'cellpairs',...
        cellpairfilter,'iterator', iterator);
        
    
%     % Run Corr - No timefilter for this. All spikes during run
%     runcorr = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter,'cellpairs',...
%         cellpairfilter,'iterator', iterator);
%     
%     % Run Corr - Place fields Only during run. linvel and noripple criterion
%     runcorr_pl = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter,'cellpairs',...
%         cellpairfilter,'excludetime', timefilter_place,'iterator', iterator);
%     
    % Run Corr - Only ripples during run. Ripple criterion only.
    runcorr_rip = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter,'cellpairs',...
        cellpairfilter,'excludetime', timefilterrun_onlyrip,'iterator', iterator);
    
    runcorr_all = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter,'cellpairs',...
        cellpairfilter,'iterator', iterator);
    
%         
%     % Run Corr - Only ripples during run. Ripple criterion and speed criterion.
%     runcorr_sp = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter,'cellpairs',...
%         cellpairfilter,'excludetime', timefilterrun_onlyrip_speed,'iterator', iterator);
    
    disp('Done Filter Creation');
    
    % Set analysis function
    % ----------------------
    xrun = setfilterfunction(xrun, 'DFAsj_calcoverlap', {'linfields', 'mapfields'},...
        'normalize',norm_overlap,'thresh',thresh_peakrate,'minbins',0.5);
   
    sleepcorr_all = setfilterfunction(sleepcorr_all, 'DFAsj_HPexpt_calcxcorrmeasures', {'spikes'}); % No time filter. But Corrln like ripples for 0.4 sec
    sleepcorr_rip = setfilterfunction(sleepcorr_rip, 'DFAsj_HPexpt_calcxcorrmeasures', {'spikes'}); % During ripples. Corrln for 0.4sec
    %sleepcorr_rip = setfilterfunction(sleepcorr_rip, 'DFAsj_calcxcorrmeasures_HpPFC', {'spikes'}); % For CA1-PFC during ripples, do not exclude PFC spikes
     
    runcorr_all = setfilterfunction(runcorr_all, 'DFAsj_HPexpt_calcxcorrmeasures', {'spikes'},'forripples',0); % To get long duration cross-corrln (1 sec instead of 0.4 sec)
    runcorr_rip = setfilterfunction(runcorr_rip, 'DFAsj_HPexpt_calcxcorrmeasures', {'spikes'}); % During ripples. Corrln for 0.4sec
    %runcorr_rip = setfilterfunction(runcorr_rip, 'DFAsj_calcxcorrmeasures_HpPFC', {'spikes'}); % For CA1-PFC during ripples, do not exclude PFC spikes
    
%     runcorr = setfilterfunction(runcorr, 'DFAsj_calcxcorrmeasures', {'spikes'},'forripples',0);
%     
%     runcorr_pl = setfilterfunction(runcorr_pl, 'DFAsj_calcxcorrmeasures', {'spikes'},'forripples',0);
%     
%     runcorr_sl = setfilterfunction(runcorr_sl, 'DFAsj_calcxcorrmeasures', {'spikes'});
%     
%     runcorr_sp = setfilterfunction(runcorr_sp, 'DFAsj_calcxcorrmeasures', {'spikes'});
    
    % Run analysis
    % ------------
    xrun = runfilter(xrun);  % 
    sleepcorr_all = runfilter(sleepcorr_all);
    sleepcorr_rip = runfilter(sleepcorr_rip);
    runcorr_all = runfilter(runcorr_all);
    runcorr_rip = runfilter(runcorr_rip);

    
%     runcorr = runfilter(runcorr);
%     runcorr_pl = runfilter(runcorr_pl);
%     runcorr_sl = runfilter(runcorr_sl);
%     runcorr_sp = runfilter(runcorr_sp);
    
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear figopt1 runscript plotdays plotanimidx savedata
        save(savefile);
    end
    
else
    
    load(savefile);
    
end  % end runscript

if ~exist('savedata')
    return
end

% ----------------------------------


% To Control Plotting, enter parameters here

if ~isempty(plotanimidx)
    useanim = plotanimidx;
else
    useanim = 1:length(xrun); % Use all animals
end
if ~isempty(plotdays)
    usedays = plotdays;
else
    usedays = [];   % Get for each animal separately
end

% ---------------------------------
% Parameters
corrwin=0.05; % Quantify - 2*50ms=100ms for ripples (50 ms on each side). For run-theta, 100ms on each side 

% ------------------------------
% Figure and Font Sizes

forppr = 1;
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/RippleInterruption/Figures/01AugSep11_RippleDisFigs/Corrln/';
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
savefig1=0;
%-----------------------------------------

% % Get trajdata and days and epochs
alloverlap = []; trajdata = []; mapdata=[]; index=[];
for anidx = 1:length(useanim)
    an = useanim(anidx);
    for i = 1:length(xrun(an).output{1}), % This should be no. of values X 2 epochs
        index{an}(i,:) = xrun(an).output{1}(i).index;
        alltrajdata1{an}{i} = xrun(an).output{1}(i).trajdata1;
        alltrajdata2{an}{i} = xrun(an).output{1}(i).trajdata2;
        allmapdata1{an}{i} = xrun(an).output{1}(i).mapdata1;
        allmapdata2{an}{i} = xrun(an).output{1}(i).mapdata2;
        
        alloverlap{an}(i) = xrun(an).output{1}(i).overlap;
        allpeakcomb{an}(i) = xrun(an).output{1}(i).peakcomb;
        alltrajpeakcomb{an}(i) = xrun(an).output{1}(i).trajpeakcomb;
    end
end

ripcorrtime=[]; % Get timebase for sleepcorr only once
for anidx = 1:length(useanim)
    an = useanim(anidx);
    for i = 1:length(sleepcorr_rip(an).output{1}), % This should be no. of values X 3 epochs
        index_sl{an}(i,:) = sleepcorr_rip(an).output{1}(i).index;
        
        % Sleep Corrln Ripples
        normsmoothcorr_slrip{an}(i,:) = sleepcorr_rip(an).output{1}(i).normsmoothcorr; 
        xc_slrip{an}{i} = sleepcorr_rip(an).output{1}(i).corr; 
        Neventscorr_slrip{an}(i) = sleepcorr_rip(an).output{1}(i).Neventscorr;
        coactivez_slrip{an}(i) = sleepcorr_rip(an).output{1}(i).coactivez;
        ec_slrip{an}(i) = sleepcorr_rip(an).output{1}(i).ec;      
%        Q_slrip{an}(i,:) = sleepcorr_rip(an).output{1}(i).corr.Q;
%        if i==1,
%           QZ = sleepcorr_rip(an).output{1}(i).corr.Z; % Z value is same for all
%        end
   
        % Sleep Corrln All
        normsmoothcorr_slall{an}(i,:) = sleepcorr_all(an).output{1}(i).normsmoothcorr; 
        xc_slall{an}{i} = sleepcorr_all(an).output{1}(i).corr; 
        Neventscorr_slall{an}(i) = sleepcorr_all(an).output{1}(i).Neventscorr;
        coactivez_slall{an}(i) = sleepcorr_all(an).output{1}(i).coactivez;
        ec_slall{an}(i) = sleepcorr_all(an).output{1}(i).ec;
%        Q_slall{an}(i,:) = sleepcorr_all(an).output{1}(i).corr.Q;
        
        if isempty(ripcorrtime)
            if isfield(sleepcorr_rip(an).output{1}(i).corr,'time');
                ripcorrtime =  sleepcorr_rip(an).output{1}(i).corr.time;
            end
        end
    end
end

runcorrtime=[]; % Get timebase for runcorr only once
for anidx = 1:length(useanim)
    an = useanim(anidx);
    for i = 1:length(runcorr_rip(an).output{1}), % This should be no. of values X 2 epochs               
        index_run{an}(i,:) = runcorr_rip(an).output{1}(i).index;
        
        % Run Corrln All
        normsmoothcorr_runall{an}(i,:) = runcorr_all(an).output{1}(i).normsmoothcorr; 
        xc_runall{an}{i} = runcorr_all(an).output{1}(i).corr; 
        Neventscorr_runall{an}(i) = runcorr_all(an).output{1}(i).Neventscorr;
        coactivez_runall{an}(i) = runcorr_all(an).output{1}(i).coactivez;
        ec_runall{an}(i) = runcorr_all(an).output{1}(i).ec;
        Q_runall{an}(i,:) = runcorr_all(an).output{1}(i).corr.Q;
                
        % Run Corrln Ripples
        if ~isnan(runcorr_rip(an).output{1}(i).normsmoothcorr)
            normsmoothcorr_runrip{an}(i,:) = runcorr_rip(an).output{1}(i).normsmoothcorr; 
        else
            normsmoothcorr_runrip{an}(i,:) = zeros(1,80);
        end
        Neventscorr_runrip{an}(i) = runcorr_rip(an).output{1}(i).Neventscorr;
        xc_runrip{an}{i} = runcorr_rip(an).output{1}(i).corr; 
        coactivez_runrip{an}(i) = runcorr_rip(an).output{1}(i).coactivez;
        ec_runrip{an}(i) = runcorr_rip(an).output{1}(i).ec;
%        if ~isnan(runcorr_rip(an).output{1}(i).corr.Q)
%           Q_runrip{an}(i,:) = runcorr_rip(an).output{1}(i).corr.Q;
%        else
%           Q_runrip{an}(i,:) = zeros(1,80);
%        end

%         % Run Corrln Place
%         normsmoothcorr_runpl{an}(i,:) = runcorr_pl(an).output{1}(i).normsmoothcorr; 
%         Neventscorr_runpl{an}(i) = runcorr_pl(an).output{1}(i).Neventscorr;

        
%         % Run Corrln Ripples+Speed Condition
%         if ~isnan(runcorr_sp(an).output{1}(i).normsmoothcorr)
%             normsmoothcorr_runsp{an}(i,:) = runcorr_sp(an).output{1}(i).normsmoothcorr; 
%         else
%             normsmoothcorr_runsp{an}(i,:) = zeros(1,400);
%         end
%         Neventscorr_runsp{an}(i) = runcorr_sp(an).output{1}(i).Neventscorr;
        
        if isempty(runcorrtime)
            if isfield(runcorr_all(an).output{1}(i).corr,'time');
                runcorrtime =  runcorr_all(an).output{1}(i).corr.time;
            end
        end
    end
end




% Initialize
clr = {'b','r','g','m','c','y','k','r'};
novel=1:3;
fam=4:8;

overlap_deltacorr = []; nooverlap_deltacorr = []; % Change in mean corrln pre-to-post-sleep (+/- 100ms)
overlap_meancorr_pre = []; nooverlap_meancorr_pre = []; % Mean corrln in sleep sessions (+/- 100ms)
overlap_meancorr_post = []; nooverlap_meancorr_post = []; % Mean corrln in sleep sessions (+/- 100ms)
overlap_runallcorr = []; nooverlap_runallcorr = []; % Mean corrln in run sessions (+/- 250ms or 100ms)
overlap_runplcorr = []; nooverlap_runplcorr = []; % Mean corrln in run sessions in place fields only (+/- 250ms pr 100ms)
overlap_runslcorr = []; nooverlap_runslcorr = []; % Mean corrln in run sessions during ripples only (+/- 100ms)
overlap_runspcorr = []; nooverlap_runspcorr = []; % Mean corrln in run sessions during ripples only with speed crit (+/- 100ms)
overlap_ov=[]; nooverlap_ov=[]; % Overlap values

midoverlap_deltacorr = []; 
midoverlap_meancorr_pre = []; 
midoverlap_meancorr_post = []; 
midoverlap_runallcorr = []; 
midoverlap_runplcorr = []; 
midoverlap_runslcorr = []; 
midoverlap_runspcorr = []; 
midoverlap_ov=[]; nooverlap_ov=[]; 




for anidx = 1:length(useanim)
    an = useanim(anidx);
    prefix = xrun(an).animal{1};
    animdirect = xrun(an).animal{2};
    if (animdirect(end) == '/')
        animdirect = animdirect(1:end-1);
    end
    
    % ******************
    % RUN DATA
    
    % Combine data for the two run epochs for each animal. 
    curridxs_antot = index{an};
    ep2idxs = find(index{an}(:,2)==2);
    ep4idxs = find(index{an}(:,2)==4); % Ep2 and Ep4 MUST be same lengths
    curridxs_an = curridxs_antot(ep2idxs,:); % Going to collapse overlap across epochs   
    %ep2idxs = 1:size(curridxs_an,1)/2;
    %ep4idxs = (size(curridxs_an,1)/2)+1:size(curridxs_an,1);
    %curridxs_an = curridxs_an(1:size(curridxs_an,1)/2,:);
    
    %1) Overlap
    overlap1 = alloverlap{an}(ep2idxs);
    overlap2 = alloverlap{an}(ep4idxs);
    %excludeidxs = unique([find(isnan(overlap1)), find(isnan(overlap2))]);
    
    excludeidxs = [find(isnan(overlap1)), find(isnan(overlap2))];
    
    %overlap1(excludeidxs)=[];
    %overlap2(excludeidxs)=[];
    overlapm{an} = mean([overlap1; overlap2]);
    
    % Also update curridxs_an, ep2idxs and ep4idxs
    
    %curridxs_an(excludeidxs,:)=[];
    %ep2idxs(excludeidxs)=[];
    %ep4idxs(excludeidxs)=[];
    
    %2) PeakRate
    peakcomb1 = allpeakcomb{an}(ep2idxs);
    peakcomb2 = allpeakcomb{an}(ep4idxs);
    peakcombm{an} = mean([peakcomb1; peakcomb2]);
    
    %3) Mapdata
    for i=1:length(ep2idxs)
        mapdata11 = allmapdata1{an}{ep2idxs(i)}.smoothedspikerate;
        mapdata12 = allmapdata1{an}{ep4idxs(i)}.smoothedspikerate;
        [x1,y1]  = size(mapdata11); [x2,y2] = size(mapdata12);
        xu = min([x1 x2]); yu = min([y1 y2]);
        mapdata1m{an}{i} = (mapdata11(1:xu,1:yu)+mapdata12(1:xu,1:yu))./2;
        mapdata11m{an}{i}=mapdata11;
        mapdata11m{an}{i}=mapdata11;
        
        mapdata21 = allmapdata2{an}{ep2idxs(i)}.smoothedspikerate;
        mapdata22 = allmapdata2{an}{ep4idxs(i)}.smoothedspikerate;
        [x1,y1]  = size(mapdata21); [x2,y2] = size(mapdata22);
        xu = min([x1 x2]); yu = min([y1 y2]);
        mapdata2m{an}{i} = (mapdata21(1:xu,1:yu)+mapdata22(1:xu,1:yu))./2;
        mapdata21m{an}{i}=mapdata21;
        mapdata21m{an}{i}=mapdata22;
        
    end
    
    % trajdata
    for i=1:length(ep2idxs)
        trajdata11 = alltrajdata1{an}{ep2idxs(i)};
        trajdata12 = alltrajdata1{an}{ep4idxs(i)};
        for traj=1:length(trajdata11)
            if length(trajdata11)>=traj && length(trajdata12)>=traj
                if ~isempty(trajdata11{traj}) && ~isempty(trajdata12{traj})
                    trajlengths = [length(trajdata11{traj}(:,5)) length(trajdata12{traj}(:,5))];
                    trajlength = min(trajlengths);
                    trajdata1m{an}{i}{traj} = mean([trajdata11{traj}(1:trajlength,5), trajdata12{traj}(1:trajlength,5)],2);
                end
            end
        end
        
        trajdata21 = alltrajdata2{an}{ep2idxs(i)};
        trajdata22 = alltrajdata2{an}{ep4idxs(i)};
        for traj=1:length(trajdata21)
            if length(trajdata21)>=traj && length(trajdata22)>=traj
                if ~isempty(trajdata21{traj}) && ~isempty(trajdata22{traj})
                    trajlengths = [length(trajdata21{traj}(:,5)) length(trajdata22{traj}(:,5))];
                    trajlength = min(trajlengths);
                    trajdata2m{an}{i}{traj} = mean([trajdata21{traj}(1:trajlength,5), trajdata22{traj}(1:trajlength,5)],2);
                end
            end
        end
    end
   
    % ******************
    % SLEEP CORR DATA
    
    curridxs_antot_sl = index_sl{an};
    ep1idxs = find(index_sl{an}(:,2)==1);
    ep3idxs = find(index_sl{an}(:,2)==3);
    ep5idxs = find(index_sl{an}(:,2)==5); % Eps 1, 3 and 5 MUST be same lengths
    curridxs_an_sl = curridxs_antot_sl(ep1idxs,:); 
    
%     % Update them based on overlap
%     % update curridxs_an_sl, ep1idxs and ep3idxs and ep5idxs

%     curridxs_an_sl(excludeidxs,:)=[];
%     ep1idxs(excludeidxs)=[];
%     ep3idxs(excludeidxs)=[];
%     ep5idxs(excludeidxs)=[];
    
    
    % ******************
    % RUN CORR DATA
    
%     curridxs_antot_run = index_run{an};
%     ep2idxs_corr = find(index_run{an}(:,2)==2);
%     ep4idxs_corr = find(index_run{an}(:,2)==4);  % Ep2 and Ep4 MUST be same lengths
%     curridxs_an_run = curridxs_antot_run(ep2idxs_corr,:); 
%     
%     % Update them based on overlap
%      % Also update curridxs_an, ep2idxs and ep4idxs
%     curridxs_an_run(excludeidxs,:)=[];
%     ep2idxs_corr(excludeidxs)=[];
%     ep4idxs_corr(excludeidxs)=[]; % These should be the same as ep2idxs and ep4idxs
    
    % All indexes should be the same for all run correlations
    
    % ******************
    % CALCULATION
    
    % Calculate for current animal
    bins = find(abs(ripcorrtime)<=corrwin); % Corrln window is Between -50 and 50 ms
    bins_slrip = find(abs(ripcorrtime)<=corrwin); % Corrln window is Between -50 and 50 ms
    bins_runall = find(abs(runcorrtime)<=2*corrwin); % Corrln window is between -100 and 100 ms
    bins_runrip = find(abs(ripcorrtime)<=corrwin); % % Corrln window is Between -50 and 50 ms for corrln during run ripples
    
    thrsev=10; thrsev_noov=thrsev;
    thrsov_high=0.2;
    thrsov_low=0.05;
    
    flagns=0;
    % Overlapping
    overlapidx = find(overlapm{an}>=thrsov_high);
 
%     for oi = 1:length(overlapidx),
%         curroi = overlapidx(oi);       
%         if flagns==0
%             % Sleep Correlations
%             Nev1 = Neventscorr_sl{an}(ep1idxs(curroi));
%             Nev3 = Neventscorr_sl{an}(ep3idxs(curroi));
%             Nev5 = Neventscorr_sl{an}(ep5idxs(curroi));
%             %if (((Nev3>=thrsev) || (Nev5>=thrsev)) && Nev1>=thrsev)
%             if ((Nev3>=thrsev) || (Nev5>=thrsev))
%                 % Save overlap value
%                 overlap_ov = [overlap_ov; overlapm{an}(curroi)];
%                 
%                 mean_ep1 = sum(normsmoothcorr_sl{an}(ep1idxs(curroi),bins)); % Take total prob in 2*corrwin window around 0
%                 mean_ep3 = sum(normsmoothcorr_sl{an}(ep3idxs(curroi),bins));
%                 mean_ep5 = sum(normsmoothcorr_sl{an}(ep5idxs(curroi),bins));
%                 
%                 % If either posts have less than thresh spikes, use the other one, otherwise use mean
%                 if ((Nev3<thrsev) || (Nev5<thrsev))
%                     if Nev3<thrsev
%                         usemean = mean_ep5;
%                     else
%                         usemean = mean_ep3;
%                     end
%                 else
%                     usemean = mean([mean_ep3, mean_ep5]); % mean of the two summed-probs
%                 end
%                 
%                 if isnan(mean_ep1);
%                     mean_ep1=0;
%                 end
%                 
%                 overlap_deltacorr = [overlap_deltacorr; usemean-mean_ep1];
%                 overlap_meancorr_post = [overlap_meancorr_post; usemean];
%                 overlap_meancorr_pre = [overlap_meancorr_pre; mean_ep1];
%                 %overlap_ec = [overlap_ec; mean([ec_sl{an}(ep3idxs(curroi)), ec_sl{an}(ep5idxs(curroi))])-ec_sl{an}(ep1idxs(curroi))];
%                 
%                 % Run Correlations: Implement under same condition for Nevs in sleep-corr - Window is 100ms each side
%                 meancorr_run_ep2 = sum(normsmoothcorr_run{an}(ep2idxs(curroi),bins_run));
%                 meancorr_run_ep4 = sum(normsmoothcorr_run{an}(ep4idxs(curroi),bins_run));
%                 overlap_runallcorr = [overlap_runallcorr; mean([meancorr_run_ep2, meancorr_run_ep4])];
%                 meancorr_runpl_ep2 = sum(normsmoothcorr_runpl{an}(ep2idxs(curroi),bins_run));
%                 meancorr_runpl_ep4 = sum(normsmoothcorr_runpl{an}(ep4idxs(curroi),bins_run));
%                 overlap_runplcorr = [overlap_runplcorr; mean([meancorr_runpl_ep2, meancorr_runpl_ep4])];
%                 % ripples during run - use bins/binsrunrip again
%                 meancorr_runsl_ep2 = sum(normsmoothcorr_runsl{an}(ep2idxs(curroi),bins_runrip));
%                 meancorr_runsl_ep4 = sum(normsmoothcorr_runsl{an}(ep4idxs(curroi),bins_runrip));
%                 overlap_runslcorr = [overlap_runslcorr; mean([meancorr_runsl_ep2, meancorr_runsl_ep4])];
%                 meancorr_runsp_ep2 = sum(normsmoothcorr_runsp{an}(ep2idxs(curroi),bins_runrip));
%                 meancorr_runsp_ep4 = sum(normsmoothcorr_runsp{an}(ep4idxs(curroi),bins_runrip));
%                 overlap_runspcorr = [overlap_runspcorr; mean([meancorr_runsp_ep2, meancorr_runsp_ep4])];
%             end
%         end
%     end
%     
%     % Non-overlapping
%     nooverlapidx = find(overlapm{an}<thrsov_low);
%     for oi = 1:length(nooverlapidx),
%         curroi = nooverlapidx(oi);
%         ep1idxs(curroi);
%         Nev1 = Neventscorr_sl{an}(ep1idxs(curroi));
%         Nev3 = Neventscorr_sl{an}(ep3idxs(curroi));
%         Nev5 = Neventscorr_sl{an}(ep5idxs(curroi));
%         %if (((Nev3>=thrsev) || (Nev5>=thrsev)) && Nev1>=thrsev)
%         if ((Nev3>=thrsev_noov) || (Nev5>=thrsev_noov))
%             % Save overlap value
%             nooverlap_ov = [nooverlap_ov; overlapm{an}(curroi)];
%             
%             mean_ep1 = sum(normsmoothcorr_sl{an}(ep1idxs(curroi),bins));
%             mean_ep3 = sum(normsmoothcorr_sl{an}(ep3idxs(curroi),bins));
%             mean_ep5 = sum(normsmoothcorr_sl{an}(ep5idxs(curroi),bins));
%             
%             % If either posts have less than thresh spikes, use the other one, otherwise use mean
%             if ((Nev3<thrsev) || (Nev5<thrsev))
%                 if Nev3<thrsev
%                     usemean = mean_ep5;
%                 else
%                     usemean = mean_ep3;
%                 end
%             else
%                 usemean = mean([mean_ep3, mean_ep5]); % mean of the two summed-probs
%             end
%             
%             if isnan(mean_ep1)
%                 mean_ep1=0;
%             end
%             
%             nooverlap_deltacorr = [nooverlap_deltacorr; usemean-mean_ep1];
%             nooverlap_meancorr_post = [nooverlap_meancorr_post; usemean];
%             nooverlap_meancorr_pre = [nooverlap_meancorr_pre; mean_ep1];
%             %nooverlap_ec = [nooverlap_ec; mean([ec_sl{an}(ep3idxs(curroi)), ec_sl{an}(ep5idxs(curroi))])-ec_sl{an}(ep1idxs(curroi))];
%             
%             % Run Correlations: Implement under same condition for Nevs in sleep-corr
%             meancorr_run_ep2 = sum(normsmoothcorr_run{an}(ep2idxs(curroi),bins_run));
%             meancorr_run_ep4 = sum(normsmoothcorr_run{an}(ep4idxs(curroi),bins_run));
%             nooverlap_runallcorr = [nooverlap_runallcorr; mean([meancorr_run_ep2, meancorr_run_ep4])];
%             meancorr_runpl_ep2 = sum(normsmoothcorr_runpl{an}(ep2idxs(curroi),bins_run));
%             meancorr_runpl_ep4 = sum(normsmoothcorr_runpl{an}(ep4idxs(curroi),bins_run));
%             nooverlap_runplcorr = [nooverlap_runplcorr; mean([meancorr_runpl_ep2, meancorr_runpl_ep4])];
%             meancorr_runsl_ep2 = sum(normsmoothcorr_runsl{an}(ep2idxs(curroi),bins_runrip));
%             meancorr_runsl_ep4 = sum(normsmoothcorr_runsl{an}(ep4idxs(curroi),bins_runrip));
%             nooverlap_runslcorr = [nooverlap_runslcorr; mean([meancorr_runsl_ep2, meancorr_runsl_ep4])];
%             meancorr_runsp_ep2 = sum(normsmoothcorr_runsp{an}(ep2idxs(curroi),bins_runrip));
%             meancorr_runsp_ep4 = sum(normsmoothcorr_runsp{an}(ep4idxs(curroi),bins_runrip));
%             nooverlap_runspcorr = [nooverlap_runspcorr; mean([meancorr_runsp_ep2, meancorr_runsp_ep4])];
% 
%         end
%     end
%     
%     
%     % Also get intermediate - for doing correlation between overlap and
%     % run/sleep correlation
%     
%     % Non-overlapping
%     midoverlapidx = find((overlapm{an}>=thrsov_low) & (overlapm{an}<thrsov_high) );
%     for oi = 1:length(midoverlapidx),
%         curroi = midoverlapidx(oi);
%         ep1idxs(curroi);
%         Nev1 = Neventscorr_sl{an}(ep1idxs(curroi));
%         Nev3 = Neventscorr_sl{an}(ep3idxs(curroi));
%         Nev5 = Neventscorr_sl{an}(ep5idxs(curroi));
%         %if (((Nev3>=thrsev) || (Nev5>=thrsev)) && Nev1>=thrsev)
%         if ((Nev3>=thrsev) || (Nev5>=thrsev))
%             % Save overlap value
%             midoverlap_ov = [midoverlap_ov; overlapm{an}(curroi)];
%             
%             mean_ep1 = sum(normsmoothcorr_sl{an}(ep1idxs(curroi),bins));
%             mean_ep3 = sum(normsmoothcorr_sl{an}(ep3idxs(curroi),bins));
%             mean_ep5 = sum(normsmoothcorr_sl{an}(ep5idxs(curroi),bins));
%             
%             % If either posts have less than thresh spikes, use the other one, otherwise use mean
%             if ((Nev3<thrsev) || (Nev5<thrsev))
%                 if Nev3<thrsev
%                     usemean = mean_ep5;
%                 else
%                     usemean = mean_ep3;
%                 end
%             else
%                 usemean = mean([mean_ep3, mean_ep5]); % mean of the two summed-probs
%             end
%             
%             if isnan(mean_ep1);
%                 mean_ep1=0;
%             end
%             
%             
%             midoverlap_deltacorr = [midoverlap_deltacorr; usemean-mean_ep1];
%             midoverlap_meancorr_post = [midoverlap_meancorr_post; usemean];
%             midoverlap_meancorr_pre = [midoverlap_meancorr_pre; mean_ep1];
%             %midoverlap_ec = [midoverlap_ec; mean([ec_sl{an}(ep3idxs(curroi)), ec_sl{an}(ep5idxs(curroi))])-ec_sl{an}(ep1idxs(curroi))];
%             
%             % Run Correlations: Implement under same condition for Nevs in sleep-corr
%             meancorr_run_ep2 = sum(normsmoothcorr_run{an}(ep2idxs(curroi),bins_run));
%             meancorr_run_ep4 = sum(normsmoothcorr_run{an}(ep4idxs(curroi),bins_run));
%             midoverlap_runallcorr = [midoverlap_runallcorr; mean([meancorr_run_ep2, meancorr_run_ep4])];
%             meancorr_runpl_ep2 = sum(normsmoothcorr_runpl{an}(ep2idxs(curroi),bins_run));
%             meancorr_runpl_ep4 = sum(normsmoothcorr_runpl{an}(ep4idxs(curroi),bins_run));
%             midoverlap_runplcorr = [midoverlap_runplcorr; mean([meancorr_runpl_ep2, meancorr_runpl_ep4])];
%             meancorr_runsl_ep2 = sum(normsmoothcorr_runsl{an}(ep2idxs(curroi),bins_runrip));
%             meancorr_runsl_ep4 = sum(normsmoothcorr_runsl{an}(ep4idxs(curroi),bins_runrip));
%             midoverlap_runslcorr = [midoverlap_runslcorr; mean([meancorr_runsl_ep2, meancorr_runsl_ep4])];
%             meancorr_runsp_ep2 = sum(normsmoothcorr_runsp{an}(ep2idxs(curroi),bins_runrip));
%             meancorr_runsp_ep4 = sum(normsmoothcorr_runsp{an}(ep4idxs(curroi),bins_runrip));
%             midoverlap_runspcorr = [midoverlap_runspcorr; mean([meancorr_runsp_ep2, meancorr_runsp_ep4])];
% 
%         end
%     end
    
    
    
    % ******************
    % PLOTTING 
    
    % Get days for current animal
    if isempty(usedays)
        days = unique(xrun(an).epochs{1}(:,1));
    else
        days = usedays;
    end
    
    % Initialize
    cnt_allpairs=0;  % Count pairs across days for current animal
    
    for d = 1:length(days)
        
        day = days(d);
        cnt_daypairs = 0;    % This is reset for each day
        
        dayidxs = find(curridxs_an(:,1)==day);
        dayidxs_sl = find(curridxs_an_sl(:,1)==day); % Should be same length as dayidxs
        
        for cellpair = 1:size(dayidxs,1)
            curridx = dayidxs(cellpair);
            
            % Sleep corrln during ripples
            plotcorr_slrip_ep1 = normsmoothcorr_slrip{an}(ep1idxs(curridx),:);
            plotcorr_slrip_ep3 = normsmoothcorr_slrip{an}(ep3idxs(curridx),:);
            plotcorr_slrip_ep5 = normsmoothcorr_slrip{an}(ep5idxs(curridx),:);                     
            Neventscorr_slrip_ep1 = Neventscorr_slrip{an}(ep1idxs(curridx));
            Neventscorr_slrip_ep3 = Neventscorr_slrip{an}(ep3idxs(curridx));
            Neventscorr_slrip_ep5 = Neventscorr_slrip{an}(ep5idxs(curridx));
%            plotQ_slrip_ep1 = Q_slrip{an}(ep1idxs(curridx),:);
%            plotQ_slrip_ep3 = Q_slrip{an}(ep3idxs(curridx),:);
%            plotQ_slrip_ep5 = Q_slrip{an}(ep5idxs(curridx),:);        
            
            % Sleep corrln during ripples
            plotcorr_slall_ep1 = normsmoothcorr_slall{an}(ep1idxs(curridx),:);
            plotcorr_slall_ep3 = normsmoothcorr_slall{an}(ep3idxs(curridx),:);
            plotcorr_slall_ep5 = normsmoothcorr_slall{an}(ep5idxs(curridx),:);                     
            Neventscorr_slall_ep1 = Neventscorr_slall{an}(ep1idxs(curridx));
            Neventscorr_slall_ep3 = Neventscorr_slall{an}(ep3idxs(curridx));
            Neventscorr_slall_ep5 = Neventscorr_slall{an}(ep5idxs(curridx));
%            plotQ_slall_ep1 = Q_slall{an}(ep1idxs(curridx),:);
%            plotQ_slall_ep3 = Q_slall{an}(ep3idxs(curridx),:);
%            plotQ_slall_ep5 = Q_slall{an}(ep5idxs(curridx),:);      
            
            
            % Run Corrln Ripples
            plotcorr_runrip_ep2 = normsmoothcorr_runrip{an}(ep2idxs(curridx),:); 
            plotcorr_runrip_ep4 = normsmoothcorr_runrip{an}(ep4idxs(curridx),:);
            plotcorr_runrip = mean([plotcorr_runrip_ep2;plotcorr_runrip_ep4]);
            Neventscorr_runrip_ep2 = Neventscorr_runrip{an}(ep2idxs(curridx));
            Neventscorr_runrip_ep4 = Neventscorr_runrip{an}(ep4idxs(curridx));
            Neventscorr_runripeps = Neventscorr_runrip_ep2 + Neventscorr_runrip_ep4;  
 %           plotQ_runrip_ep2 = Q_runrip{an}(ep2idxs(curridx),:);
 %           plotQ_runrip_ep4 = Q_runrip{an}(ep4idxs(curridx),:);
 %           plotQ_runrip = mean([plotQ_runrip_ep2;plotQ_runrip_ep4]);
          
            
            % Run Corrln All
            plotcorr_runall_ep2 = normsmoothcorr_runall{an}(ep2idxs(curridx),:); % same as ep2idxs_corr
            plotcorr_runall_ep4 = normsmoothcorr_runall{an}(ep4idxs(curridx),:);
            % Take mean
            plotcorr_runall = mean([plotcorr_runall_ep2;plotcorr_runall_ep4]);
            Neventscorr_runall_ep2 = Neventscorr_runall{an}(ep2idxs(curridx));
            Neventscorr_runall_ep4 = Neventscorr_runall{an}(ep4idxs(curridx));
            Neventscorr_runalleps = Neventscorr_runall_ep2 + Neventscorr_runall_ep4;
 %           plotQ_runall_ep2 = Q_runall{an}(ep2idxs(curridx),:); % same as ep2idxs_corr
 %           plotQ_runall_ep4 = Q_runall{an}(ep4idxs(curridx),:);
 %           plotQ_runall = mean([plotQ_runall_ep2;plotQ_runall_ep4]);
            
            
            
            
%             
%             % Run Corrln Place
%             plotcorr_runpl_ep2 = normsmoothcorr_runpl{an}(ep2idxs(curridx),:); 
%             plotcorr_runpl_ep4 = normsmoothcorr_runpl{an}(ep4idxs(curridx),:);
%             plotcorr_runpl = mean([plotcorr_runpl_ep2;plotcorr_runpl_ep4]);
%             Neventscorr_runpl_ep2 = Neventscorr_runpl{an}(ep2idxs(curridx));
%             Neventscorr_runpl_ep4 = Neventscorr_runpl{an}(ep4idxs(curridx));
%             Neventscorr_runpleps = Neventscorr_runpl_ep2 + Neventscorr_runpl_ep4;
%             
%             
%             % Run Corrln Ripples+Speed
%             plotcorr_runsp_ep2 = normsmoothcorr_runsp{an}(ep2idxs(curridx),:); 
%             plotcorr_runsp_ep4 = normsmoothcorr_runsp{an}(ep4idxs(curridx),:);
%             plotcorr_runsp = mean([plotcorr_runsp_ep2;plotcorr_runsp_ep4]);
%             Neventscorr_runsp_ep2 = Neventscorr_runsp{an}(ep2idxs(curridx));
%             Neventscorr_runsp_ep4 = Neventscorr_runsp{an}(ep4idxs(curridx));
%             Neventscorr_runspeps = Neventscorr_runsp_ep2 + Neventscorr_runsp_ep4;
            
            
            %if overlapm{an}(curridx) >= 0.2
                overlapm{an}(curridx)
                tet1 = curridxs_an(curridx,3); cell1 = curridxs_an(curridx,4);
                tet2 = curridxs_an(curridx,5); cell2 = curridxs_an(curridx,6);
                
                cnt_daypairs = cnt_daypairs + 1;
                cnt_allpairs = cnt_allpairs + 1;  % Pair Count across days - Not reset.
                pairidx{an}{day}(cnt_daypairs) = cnt_allpairs; % Very Imp - To index into processed and saved data later
                
                trajdata1 = trajdata1m{an}{curridx};
                trajdata2 = trajdata2m{an}{curridx};
                mapdata1 = mapdata11m{an}{curridx}; % 1st epoch only/ or 2nd epoch
                mapdata2 = mapdata21m{an}{curridx};
                %mapdata1 = mapdata1m{an}{curridx}; % mean of both epochs
                %mapdata2 = mapdata2m{an}{curridx};
                
                
                if figopt1==1,
                    figure(str2num([num2str(day) num2str(cnt_daypairs)]));
                    redimscreen;
                    hold on;
                    
                    % Plot rate map
                    subplot(2,3,1); hold on;
                    imagesc(flipud(mapdata1)); colorbar
                    set(gca,'YLim',[0 100]);
                    set(gca,'XLim',[0 100]);
%                     title(['Overlap ',num2str(roundn(overlapm{an}(curridx),-2)),'; Peakrate ',num2str(roundn(peakcombm{an}(curridx),-1)), ' Hz']...
%                         ,'FontSize',20,'Fontweight','normal');
                    title([prefix,' Day ',num2str(day),'  Tet ' num2str(tet1),'  Cell ',num2str(cell1)],'FontSize',20,'Fontweight','normal');
                    
                    subplot(2,3,4); hold on;
                    imagesc(flipud(mapdata2)); colorbar
                    set(gca,'YLim',[0 100]);
                    set(gca,'XLim',[0 100]);
                    xlabel ('X-position (cm)','FontSize',24,'Fontweight','bold');
                    ylabel ('Y-position (cm)','FontSize',24,'Fontweight','bold');
                    title(['Overlap ',num2str(roundn(overlapm{an}(curridx),-2)),' ; Tet ' num2str(tet2),'  Cell ',num2str(cell2)],'FontSize',20,'Fontweight','normal');
                    
                    
                    % Plot traj data / Or Corrln
                    
                    % Plot Ripple Corrln
                    subplot(2,3,2); hold on;
                    plot(ripcorrtime, plotcorr_slrip_ep1,'k','LineWidth',3);
                    %plot(ripcorrtime, plotcorr_sl_ep3,'b','LineWidth',3);
                    plot(ripcorrtime, plotcorr_slrip_ep5,'r','LineWidth',3);
                    plot(ripcorrtime, plotcorr_runrip,'g','LineWidth',3);
                   % line([0 0], [0 max([plotcorr_sl_ep1, plotcorr_runsl, plotcorr_sl_ep5])],'Color',[0.5 0.5 0.5],'LineWidth',2);
                    line([0 0], [0 max([plotcorr_slrip_ep1, plotcorr_slrip_ep5, plotcorr_runrip])],'Color',[0.5 0.5 0.5],'LineWidth',2);
                    xlabel('Time (sec)');
                    ylabel('Corrln in ripples','FontSize',20,'Fontweight','bold')
                    set(gca,'XLim',[-0.4 0.4]);
                    set(gca,'YLim',[0 max([plotcorr_slrip_ep1, plotcorr_slrip_ep5, plotcorr_runrip])+0.001]);
                    title(['Nev1(bla): ',num2str(Neventscorr_slrip_ep1),'; Nev5(red): ',num2str(Neventscorr_slrip_ep5),'; Nevrun(gre): ',num2str(Neventscorr_runripeps)]...
                        ,'FontSize',20,'Fontweight','normal');
                    %legend('Pre','Post1','Post2');
                    
                    % Plot run corrln
                    subplot(2,3,5); hold on;
                    
                    plot(runcorrtime, plotcorr_runall,'b','LineWidth',3);
                    plot(ripcorrtime, plotcorr_runrip,'g','LineWidth',3);
                    %plot(runcorrtime(bins_runrip), plotcorr_runsl(bins_runrip),'r','LineWidth',3);
                    %plot(ripcorrtime, plotcorr_runsl,'r','LineWidth',3);
                    %plot(runcorrtime(bins_runrip), plotcorr_runsp(bins_runrip),'g','LineWidth',3);
                    line([0 0], [0 max([plotcorr_runall,plotcorr_runrip])],'Color',[0.5 0.5 0.5],'LineWidth',2);
                    xlabel('Time (sec)');
                    ylabel('Run Corrln','FontSize',20,'Fontweight','bold')
                    set(gca,'XLim',[-0.4 0.4]);
                    set(gca,'YLim',[0 max([plotcorr_runall,plotcorr_runrip])+0.001]);
                    title(['Nevsrunrip(gre): ',num2str(Neventscorr_runripeps),'; Nevsrunall(blu): ',num2str(Neventscorr_runalleps)]...
                        ,'FontSize',20,'Fontweight','normal');
                    %legend('All Times','Place Only','Rip Only');
                    
                    
                    % Plot Sleep-all and Sleep Ripple Corrln on on eplot for Pre-Sleep and Post-Sleep
                    subplot(2,3,3); hold on;
                    plot(ripcorrtime, plotcorr_slrip_ep1,'k','LineWidth',3);
                    %plot(ripcorrtime, plotcorr_sl_ep3,'b','LineWidth',3);
                    plot(ripcorrtime, plotcorr_slall_ep1,'b','LineWidth',3);
                   
                   % line([0 0], [0 max([plotcorr_sl_ep1, plotcorr_runsl, plotcorr_sl_ep5])],'Color',[0.5 0.5 0.5],'LineWidth',2);
                    line([0 0], [0 max([plotcorr_slrip_ep1, plotcorr_slall_ep1])],'Color',[0.5 0.5 0.5],'LineWidth',2);
                    xlabel('Time (sec)');
                    ylabel('Pre-Sleep Corrln ','FontSize',20,'Fontweight','bold')
                    set(gca,'XLim',[-0.4 0.4]);
                    set(gca,'YLim',[0 max([plotcorr_slrip_ep1, plotcorr_slall_ep1])+0.001]);
                    title(['Nev1rip(bla): ',num2str(Neventscorr_slrip_ep1),'; Nev1all(blu): ',num2str(Neventscorr_slall_ep1)]...
                        ,'FontSize',20,'Fontweight','normal');
                    %legend('Pre','Post1','Post2');
                 
                    subplot(2,3,6); hold on;
                    plot(ripcorrtime, plotcorr_slrip_ep5,'r','LineWidth',3);
                    %plot(ripcorrtime, plotcorr_sl_ep3,'b','LineWidth',3);
                    plot(ripcorrtime, plotcorr_slall_ep5,'b','LineWidth',3);
                   
                   % line([0 0], [0 max([plotcorr_sl_ep1, plotcorr_runsl, plotcorr_sl_ep5])],'Color',[0.5 0.5 0.5],'LineWidth',2);
                    line([0 0], [0 max([plotcorr_slrip_ep5, plotcorr_slall_ep5])],'Color',[0.5 0.5 0.5],'LineWidth',2);
                    xlabel('Time (sec)');
                    ylabel('Post-Sleep Corrln ','FontSize',20,'Fontweight','bold')
                    set(gca,'XLim',[-0.4 0.4]);
                    set(gca,'YLim',[0 max([plotcorr_slrip_ep5, plotcorr_slall_ep5])+0.001]);
                    title(['Nev5rip(red): ',num2str(Neventscorr_slrip_ep5),'; Nev5all(blu): ',num2str(Neventscorr_slall_ep5)]...
                        ,'FontSize',20,'Fontweight','normal');
                    
                    % ------------------ Cross-Covariance------------------------------------------------
                    % -----------------------------------------------------------------------------------
                    
%                     figure(100*str2num([num2str(day) num2str(cnt_daypairs)]));
%                     redimscreen;
%                     hold on;
%                     
%                     % Plot All Corrln for Pre-sleep, Post-sleep and Run
%                     subplot(2,3,1); hold on;
%                     plot(ripcorrtime, plotQ_slall_ep1,'k','LineWidth',2);
%                     %plot(ripcorrtime, plotQ_sl_ep3,'b','LineWidth',2);
%                     plot(ripcorrtime, plotQ_slall_ep5,'r','LineWidth',2);
%                     plot(runcorrtime, plotQ_runall,'g','LineWidth',2);
%                    % line([0 0], [0 max([plotQ_sl_ep1, plotQ_runsl, plotQ_sl_ep5])],'Color',[0.5 0.5 0.5],'LineWidth',2);
%                     line([0 0], [0 max([plotQ_slall_ep1, plotQ_slall_ep5, plotQ_runall])],'Color',[0.5 0.5 0.5],'LineWidth',2);
%                     xlabel('Time (sec)');
%                     ylabel('COV in All','FontSize',20,'Fontweight','bold')
%                     set(gca,'XLim',[-0.4 0.4]);
%                     set(gca,'YLim',[0 max([plotQ_slall_ep1, plotQ_slall_ep5, plotQ_runall])+0.2]);
%                     title(['Nev1(bla): ',num2str(Neventscorr_slall_ep1),'; Nev5(red): ',num2str(Neventscorr_slall_ep5),'; Nevrun(gre): ',num2str(Neventscorr_runripeps)]...
%                         ,'FontSize',20,'Fontweight','normal');
%                     %legend('Pre','Post1','Post2');
%                     
%                     
%                     % Plot All Corrln for Pre-sleep, Post-sleep and Run
%                     subplot(2,3,4); hold on;
%                     plot(ripcorrtime, plotcorr_slall_ep1,'k','LineWidth',2);
%                     %plot(ripcorrtime, plotcorr_sl_ep3,'b','LineWidth',2);
%                     plot(ripcorrtime, plotcorr_slall_ep5,'r','LineWidth',2);
%                     plot(runcorrtime, plotcorr_runall,'g','LineWidth',2);
%                    % line([0 0], [0 max([plotcorr_sl_ep1, plotcorr_runsl, plotcorr_sl_ep5])],'Color',[0.5 0.5 0.5],'LineWidth',2);
%                     line([0 0], [0 max([plotcorr_slall_ep1, plotcorr_slall_ep5, plotcorr_runall])],'Color',[0.5 0.5 0.5],'LineWidth',2);
%                     line([-0.4 0.4], [QZ QZ],'Color',[0.5 0.5 0.5],'LineWidth',2);
%                     xlabel('Time (sec)');
%                     ylabel('Corrln in All','FontSize',20,'Fontweight','bold')
%                     set(gca,'XLim',[-0.4 0.4]);
%                     set(gca,'YLim',[0 max([plotcorr_slall_ep1, plotcorr_slall_ep5, plotcorr_runall])+0.001]);
%                     title(['Nev1(bla): ',num2str(Neventscorr_slall_ep1),'; Nev5(red): ',num2str(Neventscorr_slall_ep5),'; Nevrun(gre): ',num2str(Neventscorr_runripeps)]...
%                         ,'FontSize',20,'Fontweight','normal');
%                     %legend('Pre','Post1','Post2');
%          
%                     % Plot Ripple Corrln
%                     subplot(2,3,2); hold on;
%                     plot(ripcorrtime, plotQ_slrip_ep1,'k','LineWidth',2);
%                     %plot(ripcorrtime, plotQ_sl_ep3,'b','LineWidth',2);
%                     plot(ripcorrtime, plotQ_slrip_ep5,'r','LineWidth',2);
%                     plot(ripcorrtime, plotQ_runrip,'g','LineWidth',2);
%                    % line([0 0], [0 max([plotQ_sl_ep1, plotQ_runsl, plotQ_sl_ep5])],'Color',[0.5 0.5 0.5],'LineWidth',2);
%                     line([0 0], [0 max([plotQ_slrip_ep1, plotQ_slrip_ep5, plotQ_runrip])],'Color',[0.5 0.5 0.5],'LineWidth',2);
%                      line([-0.4 0.4], [QZ QZ],'Color',[0.5 0.5 0.5],'LineWidth',2);
%                     xlabel('Time (sec)');
%                     ylabel('COV in ripples','FontSize',20,'Fontweight','bold')
%                     set(gca,'XLim',[-0.4 0.4]);
%                     set(gca,'YLim',[0 max([plotQ_slrip_ep1, plotQ_slrip_ep5, plotQ_runrip])+0.2]);
%                     title(['Nev1(bla): ',num2str(Neventscorr_slrip_ep1),'; Nev5(red): ',num2str(Neventscorr_slrip_ep5),'; Nevrun(gre): ',num2str(Neventscorr_runripeps)]...
%                         ,'FontSize',20,'Fontweight','normal');
%                     %legend('Pre','Post1','Post2');
%                     
%                     % Plot run corrln
%                     subplot(2,3,5); hold on;
%                     
%                     plot(runcorrtime, plotQ_runall,'b','LineWidth',2);
%                     plot(ripcorrtime, plotQ_runrip,'g','LineWidth',2);
%                     %plot(runcorrtime(bins_runrip), plotQ_runsl(bins_runrip),'r','LineWidth',2);
%                     %plot(ripcorrtime, plotQ_runsl,'r','LineWidth',2);
%                     %plot(runcorrtime(bins_runrip), plotQ_runsp(bins_runrip),'g','LineWidth',2);
%                     line([0 0], [0 max([plotQ_runall,plotQ_runrip])],'Color',[0.5 0.5 0.5],'LineWidth',2);
%                      line([-0.4 0.4], [QZ QZ],'Color',[0.5 0.5 0.5],'LineWidth',2);
%                     xlabel('Time (sec)');
%                     ylabel('Run COV','FontSize',20,'Fontweight','bold')
%                     set(gca,'XLim',[-0.4 0.4]);
%                     set(gca,'YLim',[0 max([plotQ_runall,plotQ_runrip])+0.2]);
%                     title(['Nevsrunrip(gre): ',num2str(Neventscorr_runripeps),'; Nevsrunall(blu): ',num2str(Neventscorr_runalleps)]...
%                         ,'FontSize',20,'Fontweight','normal');
%                     %legend('All Times','Place Only','Rip Only');
%                     
%                     
%                     % Plot Sleep-all and Sleep Ripple Corrln on on eplot for Pre-Sleep and Post-Sleep
%                     subplot(2,3,3); hold on;
%                     plot(ripcorrtime, plotQ_slrip_ep1,'k','LineWidth',2);
%                     %plot(ripcorrtime, plotQ_sl_ep3,'b','LineWidth',2);
%                     plot(ripcorrtime, plotQ_slall_ep1,'b','LineWidth',2);
%                    
%                    % line([0 0], [0 max([plotQ_sl_ep1, plotQ_runsl, plotQ_sl_ep5])],'Color',[0.5 0.5 0.5],'LineWidth',2);
%                     line([0 0], [0 max([plotQ_slrip_ep1, plotQ_slall_ep1])],'Color',[0.5 0.5 0.5],'LineWidth',2);
%                      line([-0.4 0.4], [QZ QZ],'Color',[0.5 0.5 0.5],'LineWidth',2);
%                     xlabel('Time (sec)');
%                     ylabel('Pre-Sleep COV ','FontSize',20,'Fontweight','bold')
%                     set(gca,'XLim',[-0.4 0.4]);
%                     set(gca,'YLim',[0 max([plotQ_slrip_ep1, plotQ_slall_ep1])+0.2]);
%                     title(['Nev1rip(bla): ',num2str(Neventscorr_slrip_ep1),'; Nev1all(blu): ',num2str(Neventscorr_slall_ep1)]...
%                         ,'FontSize',20,'Fontweight','normal');
%                     %legend('Pre','Post1','Post2');
%                  
%                     subplot(2,3,6); hold on;
%                     plot(ripcorrtime, plotQ_slrip_ep5,'r','LineWidth',2);
%                     %plot(ripcorrtime, plotQ_sl_ep3,'b','LineWidth',2);
%                     plot(ripcorrtime, plotQ_slall_ep5,'b','LineWidth',2);
%                    
%                    % line([0 0], [0 max([plotQ_sl_ep1, plotQ_runsl, plotQ_sl_ep5])],'Color',[0.5 0.5 0.5],'LineWidth',2);
%                     line([0 0], [0 max([plotQ_slrip_ep5, plotQ_slall_ep5])],'Color',[0.5 0.5 0.5],'LineWidth',2);
%                      line([-0.4 0.4], [QZ QZ],'Color',[0.5 0.5 0.5],'LineWidth',2);
%                     xlabel('Time (sec)');
%                     ylabel('Post-Sleep COV ','FontSize',20,'Fontweight','bold')
%                     set(gca,'XLim',[-0.4 0.4]);
%                     set(gca,'YLim',[0 max([plotQ_slrip_ep5, plotQ_slall_ep5])+0.2]);
%                     title(['Nev5rip(red): ',num2str(Neventscorr_slrip_ep5),'; Nev5all(blu): ',num2str(Neventscorr_slall_ep5)]...
%                         ,'FontSize',20,'Fontweight','normal');
                    
                    
                  
                    
                    
                    
                    % Plot Traj data
%                     subplot(2,2,2); hold on;
%                     for i=1:length(trajdata1),
%                         plot(trajdata1{i},[clr{i} '.-'],'Linewidth',2);
%                     end
%                     title([prefix,' Day ',num2str(day),'  Tet ' num2str(tet1),'  Cell ',num2str(cell1)],'FontSize',20,'Fontweight','normal');
%                     legend('OutRealLeft','InRealLeft','OutRealRight','InRealRight');
%                     
%                     subplot(2,2,4); hold on;
%                     for i=1:length(trajdata2),
%                         plot(trajdata2{i},[clr{i} '.-'],'Linewidth',2);
%                     end
%                     title([prefix,' Day ',num2str(day),'  Tet ' num2str(tet2),'  Cell ',num2str(cell2)],'FontSize',20,'Fontweight','normal');
%                     xlabel ('Position along linear trajectory (cm)','FontSize',20,'Fontweight','bold');
%                     ylabel ('Firing Rate (Hz)','FontSize',24,'Fontweight','bold');
                    
                     
                     keyboard; % pause;   
                    
                end % end figopt
                
            %end % end if overlapm>thresh
        
        end % end cell pair
        
    end % end day
    
end % end animal




%*********************************************
%   FIGURES
%********************************************

if iscon==0,
    grp='Exp';
else
    grp='Con';
end



%***************************
% Distribution of Overlaps
% -------------------------
overlapvec = [];
for an=1:length(overlapm)
    overlapvec = [overlapvec, overlapm{an}];
end

figure; hold on;
hist(overlapvec,[0:0.05:0.5]);


%***************************************
% Delta Corrln for Overlap vs No Overlap
% --------------------------------------
figure; hold on; 
if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
end
    
meanov = mean(overlap_deltacorr);
errov = sem(overlap_deltacorr);
meanno = nanmean(nooverlap_deltacorr);
errno = nansem(nooverlap_deltacorr);

plot(1.1,meanno,'ro','MarkerSize',12,'LineWidth',2);
plot(2.1,meanov,'ro','MarkerSize',12,'LineWidth',2);
errorbar(1.1,meanno,errno,'r','LineWidth',2);
errorbar(2.1,meanov,errov,'r','LineWidth',2);
line([1.1 2.1], [meanno meanov],'Color','r','LineWidth',2);

[hExp,pExp] = ttest2(overlap_deltacorr,nooverlap_deltacorr);
if hExp==1,
    mul = sign(mean(overlap_deltacorr));
    plot(2.1, mean(overlap_deltacorr)+sem(overlap_deltacorr)+0.02, 'r*','MarkerSize',12);
end
title([grp,': Reactivn during rest'],'FontSize',tfont,'Fontweight','normal');
ylabel(['Delta Correlation (Mean Prob in ' num2str(corrwin*1000) 'ms win)'],'FontSize',yfont,'Fontweight','normal');
set(gca,'XTick',[1:2],'XTickLabel',{'No Overlap';'Overlap'},'FontSize',xfont,'Fontweight','normal');
text(1.4,0.15,['p = ',num2str(roundn(pExp,-3))],'FontSize',tfont,'Color','r');


%***************************************
% Corrln Values for Overlap vs No Overlap
% --------------------------------------

% Sleep Corrln Pre Beh (100ms), Run All Corr (250ms), Run Pl Corr (250ms),
% Run Ripple Corr (100ms), Sleep COrrln Post Beh
figure; hold on; 
if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
end
bar([0.8,1.8,2.8,3.8],[mean(nooverlap_meancorr_pre),mean(nooverlap_runallcorr),...
    mean(nooverlap_runplcorr),mean(nooverlap_meancorr_post)],'BarWidth',0.25,'FaceColor','b');
bar([1.2,2.2,3.2,4.2],[mean(overlap_meancorr_pre),mean(overlap_runallcorr),...
    mean(overlap_runplcorr),mean(overlap_meancorr_post)],'BarWidth',0.25,'FaceColor','m');

% bar([0.8,1.8,2.8,3.8,4.8],[mean(nooverlap_meancorr_pre),mean(nooverlap_runallcorr),...
%     mean(nooverlap_runplcorr),nanmean(nooverlap_runspcorr),mean(nooverlap_meancorr_post)],'BarWidth',0.25,'FaceColor','b');
% bar([1.2,2.2,3.2,4.2,5.2],[mean(overlap_meancorr_pre),mean(overlap_runallcorr),...
%     mean(overlap_runplcorr),mean(overlap_runspcorr),mean(overlap_meancorr_post)],'BarWidth',0.25,'FaceColor','m');

errorbar(0.8,mean(nooverlap_meancorr_pre),sem(nooverlap_meancorr_pre),'k');
errorbar(1.2,mean(overlap_meancorr_pre),sem(overlap_meancorr_pre),'k');
errorbar(1.8,mean(nooverlap_runallcorr),sem(nooverlap_runallcorr),'k');
errorbar(2.2,mean(overlap_runallcorr),sem(overlap_runallcorr),'k');
errorbar(2.8,mean(nooverlap_runplcorr),sem(nooverlap_runplcorr),'k');
errorbar(3.2,mean(overlap_runplcorr),sem(overlap_runplcorr),'k');
% errorbar(3.8,nanmean(nooverlap_runspcorr),nansem(nooverlap_runspcorr),'k');
% errorbar(4.2,mean(overlap_runspcorr),sem(overlap_runspcorr),'k');
errorbar(3.8,mean(nooverlap_meancorr_post),sem(nooverlap_meancorr_post),'k');
errorbar(4.2,mean(overlap_meancorr_post),sem(overlap_meancorr_post),'k');


% No overlap vs Overlap in Post Sleep
[hpost,ppost] = ttest2(nooverlap_meancorr_post,overlap_meancorr_post);
if hpost==1,
    plot(3.8, mean(nooverlap_meancorr_post)+sem(nooverlap_meancorr_post)+0.04, 'k+','MarkerSize',12);
    plot(4.2, mean(overlap_meancorr_post)+sem(overlap_meancorr_post)+0.04, 'k+','MarkerSize',12);
end
text(3.8,mean(overlap_meancorr_post)+sem(overlap_meancorr_post)+0.1,...
    ['p = ',num2str(roundn(ppost,-3))],'Color','k','FontSize',14);


% Overlap Pre-Sleep vs Overlap Post-Sleep
[ho,po] = ttest2(overlap_meancorr_pre,overlap_meancorr_post);
if ho==1,
    plot(1.2, mean(overlap_meancorr_pre)+sem(overlap_meancorr_pre)+0.04, 'r*','MarkerSize',12);
    plot(4.4, mean(overlap_meancorr_post)+sem(overlap_meancorr_post)+0.04, 'r*','MarkerSize',12);
end
text(4.3,mean(overlap_meancorr_post)+sem(overlap_meancorr_post)+0.01,...
    ['p = ',num2str(roundn(po,-3))],'Color','r','FontSize',16);


title([grp,': Reactivn(1,4,5) ',num2str(corrwin*1000),'ms win;  RunCorrln(2,3) ',num2str(corrwin*1000),'ms win'],'FontSize',tfont,'Fontweight','normal');
ylabel('Mean Correlation','FontSize',yfont,'Fontweight','normal');
%set(gca,'XTick',[1:4],'XTickLabel',{'Pre-Sl';'Run All';'Run Place';'Post-Sl'},'FontSize',14,'Fontweight','normal');
set(gca,'XTick',[1:5],'XTickLabel',{'Pre-Sl';'Run All';'Run Place';'Run Rip';'Post-Sl'},'FontSize',xfont,'Fontweight','normal');
set(gca,'XLim', [0.5 4.9]);
legend(' No Overlap',' Overlap','Location','NorthWest');



%****************************************************
% Overlap during run vs Correlation during run (Obvious: Confirmn) and Corrln during post-sleep
% ---------------------------------------------------

all_ov = [overlap_ov; nooverlap_ov; midoverlap_ov];  % Overlap value
allruncorr = [overlap_runallcorr; nooverlap_runallcorr; midoverlap_runallcorr]; % Corrln during all run
allruncorrpl = [overlap_runplcorr; nooverlap_runplcorr; midoverlap_runplcorr]; % Corrln during place run
allpostslcorr = [overlap_meancorr_post; nooverlap_meancorr_post; midoverlap_meancorr_post]; % Corrln during post-sleep
allpreslcorr = [overlap_meancorr_pre; nooverlap_meancorr_pre; midoverlap_meancorr_pre]; % Corrln during pre-sleep

alldiffslcorr = [allpostslcorr-allpreslcorr];

% ----------------------------------------------------------------------------------
% Prob of Correlated Firing vs Overlap for Run-All&Place and SleepRipples-PostandPre
% ----------------------------------------------------------------------------------
% 1) Scatter Plot
% ----------------
figure; hold on; 
if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
end

plot(all_ov,allruncorr,'k.','MarkerSize',16);
plot(all_ov,allruncorrpl,'b.','MarkerSize',16);
plot(all_ov,allpostslcorr,'r.','MarkerSize',16);
plot(all_ov,allpreslcorr,'m.','MarkerSize',16);

[r_runall,p_runall] = corrcoef(all_ov,allruncorr);
[r_runpl,p_runpl] = corrcoef(all_ov,allruncorrpl);
[r_postsl,p_postsl] = corrcoef(all_ov,allpostslcorr);
[r_presl,p_presl] = corrcoef(all_ov,allpreslcorr);

% bint=95% conf intervals, r=residuals, 
% stats: R^2 statistic, F statistic, p-value of F-statistic, estimate of error variance
[b,bint,r,rint,stats] = regress(allruncorr,[ones(size(all_ov)) all_ov]);         % Run All 
[b1,bint1,r1,rint1,stats1] = regress(allruncorrpl,[ones(size(all_ov)) all_ov]);  % Run Place
[b2,bint2,r2,rint2,stats2] = regress(allpostslcorr,[ones(size(all_ov)) all_ov]); % Post Sleep
[b3,bint3,r3,rint3,stats3] = regress(allpreslcorr,[ones(size(all_ov)) all_ov]);  % Pre-Sleep

xpts = 0:0.01:max(all_ov);
bfit = b(1)+b(2)*xpts;
bfit1 = b1(1)+b1(2)*xpts;
bfit2 = b2(1)+b2(2)*xpts;
bfit3 = b3(1)+b3(2)*xpts;

plot(xpts,bfit,'k-','LineWidth',2);  % Run All
plot(xpts,bfit1,'b-','LineWidth',2); % Run Place
plot(xpts,bfit2,'r-','LineWidth',2); % Post Sleep
plot(xpts,bfit3,'m-','LineWidth',2); % Pre-Sleep

set(gca,'XLim',[0 0.5]);
set(gca,'YLim',[0 0.8]);

title([grp,'ScatterPlot: AllRunCorr, PlaceCorr,PostSlCorr,PreSlCorr vs PFOverlap'],'FontSize',12,'Fontweight','normal');
xlabel(['Place Field Overlap'],'FontSize',xfont,'Fontweight','normal');
ylabel(['Total Corr Firing in ' num2str(2*corrwin*1000) ' ms win'],'FontSize',yfont,'Fontweight','normal');
%set(gca,'XTick',[1:2],'XTickLabel',{'No Overlap';'Overlap'},'FontSize',xfont,'Fontweight','normal');
%text(1.5,0.0025,['p = ',num2str(roundn(pExp,-3))],'FontSize',tfont,'Color','r');

% 2) Binned Plot
% -------------
range = [0:0.1:0.5];
ovrangex = range(1:end-1); % x-axis for plot
for i=1:length(range)-1
    rangest=range(i); rangeend=range(i+1);
    rangeidxs = find((all_ov>=rangest) & (all_ov<rangeend));
    allruncorr_range(i) = mean(allruncorr(rangeidxs)); allruncorr_rangeerr(i) = sem(allruncorr(rangeidxs));
    allruncorrpl_range(i) = mean(allruncorrpl(rangeidxs)); allruncorrpl_rangeerr(i) = sem(allruncorrpl(rangeidxs));
    allpostslcorr_range(i) = mean(allpostslcorr(rangeidxs)); allpostslcorr_rangeerr(i) = sem(allpostslcorr(rangeidxs));
    allpreslcorr_range(i) = mean(allpreslcorr(rangeidxs)); allpreslcorr_rangeerr(i) = sem(allpreslcorr(rangeidxs));
end

figure; hold on;
if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
end
errorbar(ovrangex,allruncorr_range,allruncorr_rangeerr,'k','MarkerSize',16,'LineWidth',2);
errorbar(ovrangex,allruncorrpl_range,allruncorrpl_rangeerr,'b','MarkerSize',16,'LineWidth',2);
%errorbar(ovrangex,allpostslcorr_range,allpostslcorr_rangeerr,'r','MarkerSize',16,'LineWidth',2);
%errorbar(ovrangex,allpreslcorr_range,allpreslcorr_range,'m','MarkerSize',16,'LineWidth',2);
title([grp,': AllRunCorr, PlaceCorr,PostSlCorr,PreSlCorr vs PFOverlap'],'FontSize',12,'Fontweight','normal');
xlabel(['Place Field Overlap'],'FontSize',xfont,'Fontweight','normal');
ylabel(['Total Corr Firing in ' num2str(2*corrwin*1000) ' ms win'],'FontSize',yfont,'Fontweight','normal');
set(gca,'XLim',[-0.02 0.42]);
set(gca,'YLim',[-0.05 2.5]);

if savefig1==1,
    figfile = [figdir,grp,'_PlaceCorrvsOverlap_Binned'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end

% ------------------------------------------------------------------------------------------------
% Run-Place/Run-All Prob of Correlated Firing vs SleepRipples-PostandPre Prob of Correlated Firing
% ------------------------------------------------------------------------------------------------
% 1) Scatter Plot
% ----------------
figure; hold on; 
if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
end

plot(allruncorrpl,allpostslcorr,'r.','MarkerSize',16); % Run on x-axis
plot(allruncorrpl,allpreslcorr,'m.','MarkerSize',16);
[r_runplpost,p_runplpost] = corrcoef(allruncorrpl,allpostslcorr);
[r_runplpre, p_runplpre] = corrcoef(allruncorrpl,allpreslcorr);
[b4,bint4,r4,rint4,stats4] = regress(allpostslcorr,[ones(size(allruncorrpl)) allruncorrpl]);  % Post vs RunPl
[b5,bint5,r5,rint5,stats5] = regress(allpreslcorr,[ones(size(allruncorrpl)) allruncorrpl]);  % Pre vs RunPl
xpts = 0:0.001:max(allruncorrpl);
bfit4 = b4(1)+b4(2)*xpts;
bfit5 = b5(1)+b5(2)*xpts;
plot(xpts,bfit4,'r-','LineWidth',2); % Post vs RunPl
plot(xpts,bfit5,'m-','LineWidth',2); % Pre vs RunPl
set(gca,'XLim',[0 0.8]);
set(gca,'YLim',[0 0.8]);
title([grp,': Post vs RunPl; Pre vs RunPl'],'FontSize',14,'Fontweight','normal');
xlabel(['RunPl Total Corr Firing in',num2str(2*2*corrwin*1000),' ms win'],'FontSize',xfont,'Fontweight','normal');
ylabel(['Sleep Total Corr Firing in',num2str(2*corrwin*1000),' ms win'],'FontSize',yfont,'Fontweight','normal');


% ------------------------------------------------------------------------------------------------
% The Real Thing: [Post-PreSl] vs Overlap and [Post-PreSl] vs Run-Place. Scatter as well as binned 
% ------------------------------------------------------------------------------------------------

% 1) [Post-PreSl] vs Overlap
% --------------------------
figure; hold on; 
if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
end
[r_diffsl,p_diffsl] = corrcoef(all_ov,alldiffslcorr);
% OR, [r,p] = corr(all_ov,alldiffslcorr);
% Exp: r=0.35, p=1.5e-06
% Con: r=0.31, p=2.3e-04
[b0,bint0,r0,rint0,stats0] = regress(alldiffslcorr,[ones(size(all_ov)) all_ov]);
xpts = 0:0.01:max(all_ov);
bfit0 = b0(1)+b0(2)*xpts;
plot(all_ov,alldiffslcorr,'ro','Markersize',8,'LineWidth',2);
plot(xpts,bfit0,'r-','LineWidth',4);  % diffsl vs overlap
if strcmp(grp,'Exp')    
    set(gca,'XLim',[-0.02 0.455]);
    set(gca,'YLim',[-0.5 0.7]);
else
    set(gca,'XLim',[-0.02 0.37]);
    set(gca,'YLim',[-0.4 0.5]);
end
title([grp,'ScatterPlot: AllRunCorr, PlaceCorr,PostSlCorr,PreSlCorr vs PFOverlap'],'FontSize',12,'Fontweight','normal');
xlabel(['Place Field Overlap'],'FontSize',xfont,'Fontweight','normal');
ylabel(['Diff in Rest SWR Corr Firing in ' num2str(2*corrwin*1000) ' ms win'],'FontSize',yfont,'Fontweight','normal');

if savefig1==1,
    figfile = [figdir,grp,'_RestSWRCorrDiffvsOverlap_Scatter'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end

if strcmp(grp,'Exp')    
    %range = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.5]; % combine the last 2 ranges
    range = [0:0.05:0.35];
    ovrangex = range(1:end-1); % x-axis for plot
else
    range=[0,0.05,0.1,0.15,0.2,0.4]; % combine last 2 ranges
    ovrangex = range(1:end-1);
end
for i=1:length(range)-1
    rangest=range(i); rangeend=range(i+1);
    rangeidxs = find((all_ov>=rangest) & (all_ov<rangeend));
    alldiffslcorr_range(i) = mean(alldiffslcorr(rangeidxs)); alldiffslcorr_rangeerr(i) = sem(alldiffslcorr(rangeidxs));
end    

figure; hold on;
if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
end
errorbar(ovrangex,alldiffslcorr_range,alldiffslcorr_rangeerr,'r','MarkerSize',16,'LineWidth',3);
if strcmp(grp,'Exp')    
    set(gca,'XLim',[-0.03 0.32]);
    set(gca,'YLim',[-0.03 0.3]);
else
    set(gca,'XLim',[-0.03 0.23]);
    set(gca,'YLim',[-0.03 0.38]);
end
title([grp,': Change in Rest SWR Corr Firing  vs PFOverlap'],'FontSize',12,'Fontweight','normal');
xlabel(['Place Field Overlap'],'FontSize',xfont,'Fontweight','normal');
ylabel(['Diff in Rest SWR Corr Firing in ' num2str(2*corrwin*1000) ' ms win'],'FontSize',yfont,'Fontweight','normal');

if savefig1==1,
    figfile = [figdir,grp,'_RestSWRCorrDiffvsOverlap_Binned'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end

% Do regression after shifting data to make intercept 0
alldiffslcorr_0 = alldiffslcorr-mean(alldiffslcorr); 
all_ov_0 = all_ov-mean(all_ov); 
[b0,bint0,r0,rint0,stats0] = regress(alldiffslcorr_0,[ones(size(all_ov_0)) all_ov_0]);
% figure; hold on; 
% if forppr==1
%         redimscreen_figforppr1;
%     else
%         redimscreen_figforppt1;
% end
% xpts = min(all_ov_0):0.01:max(all_ov_0);
% bfit0 = b0(1)+b0(2)*xpts;
% plot(all_ov_0,alldiffslcorr_0,'ro','Markersize',8,'LineWidth',2);
% plot(xpts,bfit0,'r-','LineWidth',4);  % diffsl vs overlap

% Regression Bar Graph
figure; hold on;
if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
end
bar(1,b0(2),'r');
errorbar(1,b0(2),b0(2)-bint0(2,1),bint0(2,2)-b0(2),'k','LineWidth',2);
%bar(2,b0(2),'b');
%errorbar(2,b0(2),b0(2)-bint0(2,1),bint0(2,2)-b0(2),'k','LineWidth',2);
%set(gca,'XLim',[-0.03 0.32]);
set(gca,'YLim',[0 1.3]);
if savefig1==1,
    figfile = [figdir,'ExpvsCon_Regression_RestSWRCorrDiffvsOverlap'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end


% 2) [Post-PreSl] vs Run-Place
% -----------------------------
figure; hold on; 
if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
end
[r_runpldiff,p_runpldiff] = corrcoef(allruncorrpl,alldiffslcorr);
% OR, [r,p] = corr(allruncorrpl,alldiffslcorr);
% Exp: r=0.50, p=1.2e-12
% Con: r=0.31, p=2.9e-04
[b00,bint00,r00,rint00,stats00] = regress(alldiffslcorr,[ones(size(allruncorrpl)) allruncorrpl]);
xpts = 0:0.01:max(allruncorrpl);
bfit00 = b00(1)+b00(2)*xpts;
plot(allruncorrpl,alldiffslcorr,'ro','Markersize',8,'LineWidth',2); % Run on x-axis
plot(xpts,bfit00,'r-','LineWidth',4);  % diffsl vs place
if strcmp(grp,'Exp')  
    set(gca,'XLim',[-0.04 1.15]);
    set(gca,'YLim',[-0.4 0.7]);
else
    set(gca,'XLim',[-0.04 1.15]);
    set(gca,'YLim',[-0.37 0.45]);
end
title([grp,': Change in Rest SWR Corr Firing  vs Run-Place Corr Firing'],'FontSize',12,'Fontweight','normal');
xlabel(['Run-Place Corr Firing'],'FontSize',xfont,'Fontweight','normal');
ylabel(['Diff in Rest SWR Corr Firing in ' num2str(2*corrwin*1000) ' ms win'],'FontSize',yfont,'Fontweight','normal');
if savefig1==1,
    figfile = [figdir,grp,'_RestSWRCorrDiffvsPlaceCorr_Scatter'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end

% Do regression after shifting data to make intercept 0
alldiffslcorr_0 = alldiffslcorr-mean(alldiffslcorr); 
allruncorrpl_0 = allruncorrpl-mean(allruncorrpl); 
[b00,bint00,r00,rint00,stats00] = regress(alldiffslcorr_0,[ones(size(allruncorrpl_0)) allruncorrpl_0]);

range = [0:0.1:0.8];
plcorrrangex = range(1:end-1); % x-axis for plot
for i=1:length(range)-1
    rangest=range(i); rangeend=range(i+1);
    rangeidxs = find((allruncorrpl>=rangest) & (allruncorrpl<rangeend));
    alldiffslcorr_range2(i) = mean(alldiffslcorr(rangeidxs)); alldiffslcorr_range2err(i) = sem(alldiffslcorr(rangeidxs));
end
figure; hold on;
if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
end
errorbar(plcorrrangex,alldiffslcorr_range2,alldiffslcorr_range2err,'r','MarkerSize',16,'LineWidth',3);
set(gca,'XLim',[-0.04 0.75]);
set(gca,'YLim',[-0.04 0.6]);
title([grp,': Change in Rest SWR Corr Firing  vs Run-Place Corr Firing'],'FontSize',12,'Fontweight','normal');
xlabel(['Run-Place Corr Firing'],'FontSize',xfont,'Fontweight','normal');
ylabel(['Diff in Rest SWR Corr Firing in ' num2str(2*corrwin*1000) ' ms win'],'FontSize',yfont,'Fontweight','normal');
if savefig1==1,
    figfile = [figdir,grp,'_RestSWRCorrDiffvsPlaceCorr_Binned'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end

% Regression Bar Graph
figure; hold on;
if forppr==1
        redimscreen_figforppr1;
    else
        redimscreen_figforppt1;
end
bar(1,b00(2),'r');
errorbar(1,b00(2),b00(2)-bint00(2,1),bint00(2,2)-b00(2),'k','LineWidth',2);
%bar(2,b00(2),'b');
%errorbar(2,b00(2),b00(2)-bint00(2,1),bint00(2,2)-b00(2),'k','LineWidth',2);
%set(gca,'XLim',[-0.03 0.32]);
set(gca,'YLim',[0 0.35]);
if savefig1==1,
    figfile = [figdir,'ExpvsCon_Regression_RestSWRCorrDiffvsPlaceCorr'];
    print('-dpdf', figfile);
    print('-djpeg', figfile);
    saveas(gcf,figfile,'fig');
end

% ----------------------------------------------
i=1;
keyboard;

