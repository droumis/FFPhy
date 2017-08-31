

% rundecode_four (kk_fourdecode)  >> produces output plotted here

% Trial by trial plot of Outbound to CP

load('/opt/data13/kkay/Unitgroups/Unitgroup_Principal_16_animals.mat','Principal_16_animals');
load('/opt/data13/kkay/Unitgroups/Unitgroup_Principal_CA1_16_animals.mat','Principal_CA1_16_animals');
%load('/home/kkay/Src/Matlab/_KennyFunctions/BBR_3_cmap.mat','cmap')
    %colormap(cmap);

    
FLIP = 1;
OFFSET_ANG = {'Frank',[0 0]};    

decode_dir      = '/opt/data13/kkay/___Superfourcode_data/Decode';
decode_dir_2    =  '/opt/data13/kkay/___Superfourcode_data/Decode2';

Codes_toplot = [1 5];

BREAKPLOT = 0;

ade = [7 5 6];  % an day ep  %% fantastic
ade = [7 4 6];
%ade = [6 11 4 ];  % 18 R -- 3 representations  -- lots of triplets
%ade = [1 9 2];  % 30 L -- great example
%ade = [1 9 3];
%ade = [6 8 4];  % ** Premier 11 L, 17 

%ade = [6 10 2];  % 17 L -- check headdir
%ade = [6 10 4];  % 3 representations.. incredible
%ade = [6 11 2];  % 18 R -- 3 representations

%ade = [6 11 6];  % quadruplet ensemble

%ade = [2 8 2];  % 17 L

%ade = [7 10 4];
%ade = [8 4 2];

%ade = [6 11 4];  % an day ep
%ade = [6 12 2];  % an day ep  %% fantastic

%ade = [7 4 6];  % an day ep  %% fantastic

%ade = [14 2 4];

% Plot parameters
adtc_toplot = [];  % specify adtc or leave blank to plot all
units_tofilter = Principal_16_animals;
windowsize = 10;
principal_flag = 0;
lopass = [];
place_thresh_fr = 2;
place_thresh_size = 6;
CELLORDER = 1;  % sort by spatial index

plot_altspeed = 1;
plot_clusterless = 0;
plotDECODE = 1;
plotEEG = 1;

% Scales
HDSCALE = 850;
EEGSCALE = 10 ;  %-0.2; -.2; %-.4 %.5;
VELSCALE = 30; %35;
ANGVELSCALE = 300;
plot_sleepblock = 0;
ymax = 2143;
eeg_ylevel = 1700;
    yoffset_eeg = 100;  % shifts the EEG upwards (to avoid crowding)
rip_ylevel = 1300; %1025;
rasterwidth = 2; %.6;
eeglinewidth = .8;
riplinewidth = .8;
CPclr = [.87 .95 .87];

% Region settings
reginfo{1} = {'CA1','CA2','CA3','DG','unknown'};
reginfo{2} = {[0 0 0],[0 .7 0],[1 0 0],[1 0 1],[.7 .7 .7]};

%%%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander',...
                        'Iceland','Justice','Kapital','Laplace','Eight','Ten','Conley','Miles'};
an = ade(1);
    animalname = animals_order{an};
    animalinfo = animaldef(animalname);
        animaldir = animalinfo{2};
        animalprefix = animalinfo{3};
d = ade(2);
ep = ade(3);

% basic data %%
cellinfo = loaddatastruct(animaldir, animalprefix, 'cellinfo');
spikes = loaddatastruct(animaldir, animalprefix, 'spikes', d); 
linpos = loaddatastruct(animaldir, animalprefix, 'linpos', d);
pos = loaddatastruct(animaldir, animalprefix, 'pos', d);
    postimevec = pos{d}{ep}.data(:,1);
        pos_fs = 1/(postimevec(2)-postimevec(1));
        windowsize_possamp = round(windowsize * pos_fs);
    epstart = postimevec(1);    
    if size(pos{d}{ep}.data,2) == 7
        vel = pos{d}{ep}.data(:,6);
        hdir = pos{d}{ep}.data(:,4); % in radians
        xpos  = pos{d}{ep}.data(:,[2]); % in radians
        ypos  = pos{d}{ep}.data(:,[3]); % in radians
    elseif size(pos{d}{ep}.data,2) == 9
         vel = pos{d}{ep}.data(:,9);
         hdir = pos{d}{ep}.data(:,8); % in radians
         xpos = pos{d}{ep}.data(:,6); % in radians
         ypos = pos{d}{ep}.data(:,7); % in radians
    elseif size(pos{d}{ep}.data,2) == 5
         vel = pos{d}{ep}.data(:,5);
         hdir = pos{d}{ep}.data(:,4); % in radians
         xpos = pos{d}{ep}.data(:,2); % in radians
         ypos = pos{d}{ep}.data(:,3); % in radians        
    end
    if 0
        kernel = gaussian(15,15*8);
        vels = smoothvect(vels,kernel);
    end
    
if plotDECODE
    
    cd(decode_dir)
    filename = dir(sprintf('%s_%d_%d.mat',animalname(1:3),d,ep));
    load(filename.name,'P')    
    
    tbins = loaddatastruct(animaldir,animalprefix,'thetabins',d);
        thetatet = tbins{d}.thetatet;
        if plotEEG
            dwaves = loaddatastruct(animaldir, animalprefix, 'dwaves', d);
                deltatet = dwaves{d}{ep}.reftet;  % DG ref (or CA3)        
        % load eeg used to create decoding bins
        theta = loadeegstruct(animaldir,animalprefix,'theta',d,ep,thetatet);
            thetatrace = theta{d}{ep}{thetatet}.data(:,1);
            thetatimevec = geteegtimes(theta{d}{ep}{thetatet});
        eeg = loadeegstruct(animaldir,animalprefix,'eeg',d,ep,thetatet);
            eegtrace = eeg{d}{ep}{thetatet}.data;
            eegtimevec = geteegtimes(eeg{d}{ep}{thetatet});
        end
        
    binlist = P.binlist;
    binvec_c = P.binvec_c;
    numbins = P.numbins;
    bindurs = P.bindurs;
    
    twovec      = cell(1,5);
    sivec       = cell(1,5);
    si_clr      = cell(1,5);
    altspeed    = cell(1,5);
    
    % Sorted %%%%%%%%%%%%%%%%%%%%%%%%%% 
    for C = 1:5
        
        % "Place map" vector
        twovec{C} = P.posteriors{1}{C};
        
        % Place map Distribution (Bernoulli)
        sivec{C} = twovec{C}(:,1)./sum(twovec{C},2);
        
        % Alternation speed
        altspeed{C} = [nan ; diff(sivec{C})];
            % abs + normalize
            altspeed{C} = abs(altspeed{C});   % 1 is max change
        
        % Color vector
        if C == 1
            load('/home/kkay/Src/Matlab/_KennyFunctions/BYR_cmap.mat','cmap')
            cmap2 = flipud(cmap);  % for head dir
            colormap(cmap);
        else
            cmap = colormap('gray');
            cmap = flipud(cmap);
        end
           
        
        si_clr{C} = nan(size(sivec{C},1),3); % get color vector
        numclr = size(cmap,1);
        for zz = 1:size(sivec{C},1)
            if ~isnan(sivec{C}(zz))
                jj = round(numclr * sivec{C}(zz));
                if jj == 0
                    jj = 1;
                end
                si_clr{C}(zz,:) = cmap(jj,:);
            end
        end
        
    end
   
    % (optional) Clusterless %%%%%%%%%%%%%%%%%%%%%%%%%%
    if plot_clusterless
        cd(decode_dir_2)
        filename = dir(sprintf('%s__%d_%d_fullep.mat',animalname(1:3),d,ep));
        load(filename.name,'P')
        
        % time vector
        tvec2 = P.xvecms/1000;
        % get SI vector
        twovec = P.posteriors{1}';
        sivec2 = twovec(:,1)./sum(twovec,2);
        % get color vector
        si_clr2 = nan(size(sivec2,1),3);
        numclr = size(cmap,1);
        for zz = 1:size(sivec2,1)
            if ~isnan(sivec2(zz))
                jj = round(numclr * sivec2(zz));
                if jj == 0
                    jj = 1;
                end
                si_clr2(zz,:) = cmap(jj,:);
            end
        end
    end
    
end
    
    
 % angular velocity %
    % get offset

    %OFFSET_ANG{find(strcmp(animalname,OFFSET_ANG))+1}

    % interpolate hdir
 naninds = isnan(hdir);
 hdir = interp1(postimevec(~naninds),hdir(~naninds),postimevec,'linear');
 
 hdir2 = hdir;
 posinds = hdir2 > 0;
 neginds = hdir2 < 0;
        hdir2(posinds) = hdir2(posinds) - pi;
        hdir2(neginds) = hdir2(neginds) + pi;

 hdir_revo_1 = hdir  / (2*pi);
 hdir_revo_2 = hdir2 / (2*pi);   % plotted on top
        % get corresponding colormap points
        hdir_amp = hdir_revo_2 + 0.5;
        hdir_clr = nan(length(hdir_amp),3);
            mx = max(hdir_amp);  % just 1 right
            numclr = size(cmap,1);
        for pn = 1:size(hdir_clr,1)         
           jj = round(numclr * hdir_amp(pn) / mx);
           if jj == 0
               jj = 1;
           end
           hdir_clr(pn,:) = cmap(jj,:);
        end        
 angvel1 = abs(diff(hdir_revo_1(:))) * pos_fs;
 angvel2 = abs(diff(hdir_revo_2(:))) * pos_fs;
    angvel = min([angvel1(:) angvel2(:)],[],2); 
    if 0
       figure; hold on
       plot(angvel1,'r')
       plot(angvel2,'b')
       plot(angvel,'k');
    end
 
 % linear position data %   
CPbuff = choicepointer(linpos{d}{ep});
lindist = linpos{d}{ep}.statematrix.lindist;
CPvec = (lindist > CPbuff-10) & (lindist < CPbuff );   % 10 cm zone
CAvec = (lindist < CPbuff) ;        % Center Arm
CPperiods = vec2list(CPvec,postimevec);
CAperiods = vec2list(CAvec,postimevec);

% place field data %% 
if ~exist('superlin','var') || ~strcmp(superlin.animalname,animalname)
    superlinfields_dir = '/opt/data13/kkay/Superlin_data/';
    filename = dir(sprintf('%sSuperlin_%s.*',superlinfields_dir,animalname(1:3)));
    load([superlinfields_dir filename(end).name],'superlin')
end

   
% Obtain unit data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 0. Identify units
if isempty(adtc_toplot)
    [adtc_list] = clusteredunits(spikes,cellinfo,an,d,ep,[],principal_flag);
else
    adtc_list = adtc_toplot;
end

if ~isempty(units_tofilter)
    inds_filt = ismember(adtc_list,units_tofilter,'rows');
    adtc_list = adtc_list(inds_filt,:);
    numfilteredout = sum(~inds_filt);
    disp(sprintf('not plotting %d units as specified by filter',numfilteredout))
end

    % initialize outputs
numcells = size(adtc_list,1);
    spiketimes      = cell(numcells,1);
    regs            = nan(numcells,1);
    spatialindex    = nan(numcells,1);   % -1 is LEFT, 0 is BIDIR, +1 is RIGHT
        spatialclr  = nan(numcells,3);
    
for uu = 1:numcells
    
    tet = adtc_list(uu,3);
    cellnum = adtc_list(uu,4);
    
    % 1. Spike times %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    spkdata = spikes{d}{ep}{tet}{cellnum}.data;
    if ~isempty(spkdata)
        spiketimes{uu} = spikes{d}{ep}{tet}{cellnum}.data(:,1);
    else
        spiketimes{uu} = []; 
    end
    
    % 2. Area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cellentry = cellinfo{d}{ep}{tet}{cellnum};
    if ~isfield(cellentry,'area')
        regs(uu) = 5;
    else
        regstr = cellentry.area;
        regnum = find(strcmp(regstr,reginfo{1}));
        if isempty(regnum)
            regnum = 5;
        end
        regs(uu) = regnum;
    end
    
    % 3. Spatial index %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SLSTATE = 2;
    [ind_sl,regnum] = superlinfind([d ep tet cellnum],superlin,SLSTATE);
    if ind_sl == 0
         spatialclr(uu,:) = [1 1 1];  % White means can't find in superlin
        continue
    end
    
    dataentry = superlin.data{regnum}.detc_output{SLSTATE}(ind_sl);

    occnorm = [nan nan];
    placefield_flag = [0 0];
    for RL = 1:2   % 1: R, 2: L
        if RL == 1
            tr = 1;
        elseif RL == 2
            tr = 3;
        end
        if tr <= length(dataentry.trajdata) && ...
                   ~isempty(dataentry.trajdata{tr})
            
            datamat = dataentry.trajdata{tr};   % 1st col: location, 2nd col: raw occ, 3rd col: raw spk count, 5: smoothed place map
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
                    placefield_flag(RL) = 1;  % flag that p.f. detected
                end
            end
            
            % Tabulation of spikes and time occupied
            totalspk = sum(spkcnt);
            totaltime = sum(occ);
                occnorm(RL) = totalspk / totaltime;
        end
    end
    
    % Spatial index is R / (R + L)  and scaled to -1 (L) to 0 (bidir) to 1 (R)
    if any(placefield_flag)
        spatialratio = occnorm(1) / (occnorm(1) + occnorm(2));
        spatialindex(uu) = 2 * (spatialratio  - 0.5 ) ;
            % assign an RGB value from a 3-color scheme 
            ind_cmap = round(spatialratio * size(cmap,1));
            if ind_cmap == 0
                ind_cmap = 1;
            end
            spatialclr(uu,:) = cmap2(ind_cmap,:);
    else
        spatialindex(uu) = -99999;
        spatialclr(uu,:) = [0 .7 0];  % Green means doesn't have place field on Trajs 1, 3
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

% if sort specified
if CELLORDER == 1  % sort by spatialindex
   [spatialindex,ord] = sort(spatialindex,'descend');
   spatialclr = spatialclr(ord,:);
   spiketimes = spiketimes(ord);
   regs = regs(ord);
   adtc_list(ord,:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine Prospective Trajectory periods  (Right and Left) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wtrackallsegs = loaddatastruct(animaldir, animalprefix, 'wtrackallsegs', d);
inds_all = ismember(wtrackallsegs{d}{ep}.trajseg(:,3),[1 3]); 
inds{1} = ismember(wtrackallsegs{d}{ep}.trajseg(:,3),[1]); 
inds{2} = ismember(wtrackallsegs{d}{ep}.trajseg(:,3),[3]);

travs = [];                 
travs2 = cell(1,2);

for T = 1:2             % 1: Right, 2: Left
    
    segblock = wtrackallsegs{d}{ep}.trajseg(inds{T},[1 2]) ;
        numtrv = size(segblock,1);
    
    % Determine time when animal crosses choice point (CP) for the first
    % time this traversal
    crosstimes = [];
    for s = 1:numtrv
        
        time_a = segblock(s,1);  % beginning of traversal
        ind_b = segblock(s,2);  % end of traversal
            a = lookup(time_a,postimevec);
            ind_b = lookup(ind_b,postimevec);
        lindist_trav = lindist(a:ind_b);
        ptimevec_trav = postimevec(a:ind_b);
        
        pastCPinds = find(lindist_trav > CPbuff);
            % c: ind of crossing
            c = min(pastCPinds);            
            % time elapsed since beginning of traversal to time of crossing
            time_elapsed = ptimevec_trav(c) - ptimevec_trav(1);
            
        crosstime = time_a + time_elapsed;
        crosstimes = [crosstimes ; crosstime];
    end
    
    travs =     [travs     ;   segblock  crosstimes  T*ones(numtrv,1) ];
    travs2{T} = [travs2{T} ;   segblock  crosstimes                    ];
    
end

[~,ord] = sort(travs(:,1),'ascend');
travs = travs(ord,:);
numtravs = size(travs,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





  

%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotcount = 1;
trv = 1;

K = figure('units','normalized','outerposition',[.1 .05 .65 .9]);

maxnumplot = 3;

YMAX = 10000;

while trv < numtravs
    
    tstart = travs(trv,1);
    tend = travs(trv,2);
        travdur = tend-tstart;
    tcross = travs(trv,3);
    RL = travs(trv,4);
        % R vs. L colors + background
        if RL == 1       % Right is blue-based
            RLstr = 'R';
            trajclr1 = [1 1 1];
            trajclr2 = [.2 .5 1];
            trajclr3 = 'b';
        elseif RL == 2   % Left is red-based
            RLstr = 'L';
            trajclr1 = [1 1 1];
            trajclr2 = [1 .3 .3];
            trajclr3 = 'r';
        end

    if plotcount > maxnumplot
        plotcount = 1;
        if BREAKPLOT
            break
        end
        %keyboard
        %break
        %pause
        %close all
        K = figure('units','normalized','outerposition',[.1 .05 .65 .9]);
    end
    
    pos_a = lookup(tstart,postimevec);
    pos_b = lookup(tend,postimevec);
    xvec_pos = postimevec(pos_a:pos_b) - postimevec(pos_a);

    % Main subplot
    mainplots{1} = [1 2 3 4];
    mainplots{2} = [7 8 9 10];
    mainplots{3} = [13 14 15 16];
    subplot(maxnumplot,6,mainplots{plotcount})
        set(gca,'color',trajclr1); hold on
    
    % CP line
    cpl = plot([tcross tcross]-tstart,[0 YMAX],'-','color',[0 .7 0],'linewidth',2);
    uistack(cpl,'bottom')
    
    % Center junction ("crossing") periods 
    inds_cross = (CPperiods(:,1) > tstart) & (CPperiods(:,1) < tend);
    crosspers = CPperiods(inds_cross,:);
    for cps = 1:size(crosspers,1)
        cpstart = crosspers(cps,1)-tstart;
        cpend = crosspers(cps,2)-tstart;
        patch([cpstart cpend cpend cpstart],[0 0 YMAX YMAX],CPclr,'edgecolor','none'); hold on
    end
    
    % Speed 
    area(xvec_pos-xvec_pos(1),...
        VELSCALE*vel(pos_a:pos_b),'linewidth',1,'facecolor',[.9 .9 .9]); hold on
    plot(xvec_pos-xvec_pos(1),...
        VELSCALE*vel(pos_a:pos_b),'-','color',[.6 .6 .6],'linewidth',2); hold on
    
    % Angular speed 
    if 0
    plot(xvec_pos-xvec_pos(1),...
        ANGVELSCALE*angvel(pos_a:pos_b),'-','color',[.4 .25 .4],'linewidth',3); hold on
    end
    
    % Spikes 
    spk_y = 100;  % dimension
    YOFFSET = 400;
    for ccc = 1:numcells

        spikeinds = logical(isExcluded(spiketimes{ccc}, [tstart tend]));
        
        ycoord = YOFFSET + (ccc)*spk_y;
        
        if any(spikeinds)
            for spt = spiketimes{ccc}(spikeinds)'
                st = spt-tstart;  
                h1 = plot([st st], [ycoord  ycoord+spk_y],'linewidth',rasterwidth,'Color',spatialclr(ccc,:));
            end
        end
        
        regclr = reginfo{2}{regs(ccc)};
        
        xl = xlim;
        
        % print unit # on plot
        xpltcoord = .01*(xl(2)-xl(1))-xl(1);
        text(xpltcoord,ycoord + spk_y/3,sprintf('%d %d',adtc_list(ccc,[3 4])),'fontsize',7,'fontweight','bold','color',regclr)        
    end
  
    ycoord = ycoord + 800;
    
    % Decodes %%%
    if plotDECODE
        
        % Sorted %%%%%%%%%%%%%%%%%%%%
        numCTP = length(Codes_toplot);
        for CNUM = 1:numCTP
            
            Code = Codes_toplot(CNUM);
            
            % guide lines
            plot([tstart tend]-tstart,ycoord + [HDSCALE HDSCALE],'-','Color',[.8 .8 .8],'linewidth',1)
            plot([tstart tend]-tstart,ycoord - [0 0],'-','Color',[.8 .8 .8],'linewidth',1)
            %plot([tstart tend]-tstart,ycoord + HDSCALE*[0.5 0.5],'--','color',[.8 .8 .8],'linewidth',1)
            
            % plot
            aaa = lookup(tstart,binvec_c);
            bbb = lookup(tend,binvec_c);
            
            % Plot decodes as patches
            for bn = aaa:bbb
                ind_a = binlist(bn,1) - tstart;
                ind_b = binlist(bn,2) - tstart;
                
                % if Cen dir decode, then omit plotting decodes
                if Code == 5
                    if ~any( logical ( isExcluded( binlist(bn,:),CAperiods)))
                        continue
                    end
                end
                
                if any(isnan(si_clr{Code}(bn,:)))
                    continue
                else
                    ylo = ycoord+0*HDSCALE;
                    yhi = ycoord+1*HDSCALE;
                    patch([ind_a ind_b ind_b ind_a],[ylo ylo yhi yhi],si_clr{Code}(bn,:),'edgecolor','k'); hold on
                end
            end
            
            % Plot decodes as dots
            % scatter(binlist(aaa:bbb) - binlist(aaa),...
            %        ycoord + HDSCALE*sivec{Code}(aaa:bbb),300,si_clr{Code}(aaa:bbb,:),'.'); hold on
            
            ycoord = ycoord + 1100;            
            
        end
            
%         % Clusterless %%%%%%%%%%%%%%%%%%%%
%         if plot_clusterless
%             ycoord = ycoord + 1100;
%             % guide lines
%             plot([tstart tend]-tstart,ycoord + [HDSCALE HDSCALE],'-r','linewidth',1)
%             plot([tstart tend]-tstart,ycoord - [0 0],'-b','linewidth',1)
%             plot([tstart tend]-tstart,ycoord + HDSCALE*[0.5 0.5],'--','color',[.8 .8 .8],'linewidth',1)
%             % plot
%             aaa = lookup(tstart,tvec2);
%             bbb = lookup(tend,tvec2);
%             scatter(tvec2(aaa:bbb)-tvec2(aaa),...
%                 ycoord + HDSCALE*sivec2(aaa:bbb),300,si_clr2(aaa:bbb,:),'.'); hold on
%         end

    end
    
   
    % Plot Alternation "Speed" %%%%%%%%%%%%%%%%%
    
    if plot_altspeed
        
        altspeed;
        binvec_c;
        
        ALTSCALE = 1000;
        
        ACODE = 1;
        
        % guide lines
        plot([tstart tend]-tstart,ycoord + [ALTSCALE ALTSCALE],'-','Color',[.8 .8 .8],'linewidth',1)
        plot([tstart tend]-tstart,ycoord - [0 0],'-','Color',[.8 .8 .8],'linewidth',1)
        
        % plot
        aaa = lookup(tstart,binvec_c);
        bbb = lookup(tend,binvec_c);
        
        % Plot speed as black dots
        for bn = aaa:bbb
            
            if ~isnan(altspeed{ACODE}(bn))
                %ylo = ycoord+0*HDSCALE;
                %yhi = ycoord+1*HDSCALE;
                %patch([ind_a ind_b ind_b ind_a],[ylo ylo yhi yhi],altspeed{Code}(bn,:),'edgecolor','none'); hold on
                x_coord = binvec_c(bn) - tstart;
                y_coord = ycoord + 1000*altspeed{ACODE}(bn);
                scatter(x_coord,y_coord,350,'k','.');
               
            end
            
        end
        
        ycoord = ycoord + 800;
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Head direction trace %%%
    TL = 0.05;
    ycoord = ycoord + 800;
        % headdir bars
        headmax = HDSCALE*max(hdir_revo_2);
        headmin = HDSCALE*min(hdir_revo_2);
            headmid = (headmax + headmin) / 2 ;
            % side bars
    plot([0 TL],ycoord + [headmax headmax],'linewidth',1,'Color','k'); hold on;
    plot([0 TL],ycoord + [headmin headmin],'linewidth',1,'Color','k'); hold on;
            % long bars
    plot([0 travdur],ycoord + [headmax headmax],'-','linewidth',1,'Color',[.8 .8 .8]); hold on;
    plot([0 travdur],ycoord + [headmin headmin],'-','linewidth',1,'Color',[.8 .8 .8]); hold on;
    plot([0 travdur],ycoord + [headmid headmid],'--','linewidth',1,'Color',[.8 .8 .8]); hold on;     
        % headdir
    if 0
        plot(xvec_pos-xvec_pos(1),...
            ycoord + HDSCALE*hdir_revo_2(pos_a:pos_b),'-','color',[.2 .2 .2],'linewidth',2); hold on
    elseif 0
        plot(xvec_pos-xvec_pos(1),...
            ycoord + HDSCALE*hdir_revo_2(pos_a:pos_b),'.','color',[.2 .2 .2],'markersize',15); hold on
    elseif 1  %hdir_clr
        scatter(xvec_pos-xvec_pos(1),...
                ycoord + HDSCALE*hdir_revo_2(pos_a:pos_b),300,hdir_clr(pos_a:pos_b,:),'.'); hold on
    end
    axis tight
    
    ycoord = ycoord + 550;
    
    % plot EEG %%%%%%%%%%%%%%%%%%%%%%%%
    if plotEEG
        
        ind_a = lookup(tstart,eegtimevec)      ;
        ind_b = lookup(tend,eegtimevec)      ;
        a = lookup(tstart,thetatimevec)      ;
        b = lookup(tend,thetatimevec)      ;
        
        xvec_eeg = eegtimevec(ind_a:ind_b)-eegtimevec(ind_a)  ;
        etrace = eegtrace(ind_a:ind_b);
        
        xvec = thetatimevec(a:b)-thetatimevec(a)  ;
        ttrace = thetatrace(a:b)   ;
        
%         if ~isempty(lopass)
%             eegtrace=eegfilt(eegdata(ind_a:ind_b)',1500,lopass,0); % this is just filtering  
%             %efilt = designeegfilt(1500,10,400);
%             %etrace = filtfilt(efilt, 1 , eegdata(a:b));
%         end

        yoffset_eeg = -ycoord / 1.4;

        % plot eeg
        if 0
            h2 = plot(xvec_eeg, ycoord + yoffset_eeg +  EEGSCALE * double(etrace),...
                      'Color',[.8 .8 .8],'linewidth',eeglinewidth); hold on
        end
        % plot theta
        h1 = plot(xvec, ycoord + yoffset_eeg + EEGSCALE*ttrace,'Color',[.8 .8 .8],'linewidth',eeglinewidth); hold on
 
        uistack(h1,'bottom');
 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    % buffer on Y
    ylim([0 ycoord])
    
    
    % Scale bars
        % velocity scale bar
    TL1 = 0.05;
    plot([0 travdur],[20*VELSCALE 20*VELSCALE],'--','linewidth',1,'Color',[.4 .4 .4]); hold on;
    plot([0 0],[0 20*VELSCALE],'linewidth',1,'Color','k'); hold on;
    plot([0 TL1],[20*VELSCALE 20*VELSCALE],'linewidth',1,'Color','k'); hold on;
    plot([0 0],[0 40*VELSCALE],'linewidth',1,'Color','k'); hold on;
    plot([0 TL1],[40*VELSCALE 40*VELSCALE],'linewidth',1,'Color','k'); hold on;
    plot([travdur travdur],[0 20*VELSCALE],'linewidth',1,'Color','k'); hold on;
    plot([travdur-TL1 travdur],[20*VELSCALE 20*VELSCALE],'linewidth',1,'Color','k'); hold on;
    plot([travdur travdur],[0 40*VELSCALE],'linewidth',1,'Color','k'); hold on;
    plot([travdur-TL1 travdur],[40*VELSCALE 40*VELSCALE],'linewidth',1,'Color','k'); hold on;
   
    
    % main plot formatting
    set(gca,'fontsize',14)
    set(gca,'xtick',0:1:travdur)
    set(gca,'ticklength',get(gca,'ticklength'))
    set(gca,'Layer','top')
    set(gca,'ytick',[0 20*VELSCALE 40*VELSCALE])
    set(gca,'yticklabel',[0 20 40])

    
    % main plot title
    title(sprintf('%s %d %d (Trav: %d %s)  (epstart: %d)',...
            animalname(1:3),d,ep,trv,RLstr,round(tstart)),'fontsize',16,'fontweight','bold')
    colormap(cmap)
    % % colorbar
    %cb = colorbar;
    %    set(cb,'YTick',[])
    %    set(cb,'position',[.62 .765 .005 .09])
            
    if trv == numtravs || plotcount > maxnumplot
        xlabel('Time (s)','fontweight','bold','fontsize',15)
    end    

    
    % Side plot: 2D map 
    sideplots{1} = [5 6];
    sideplots{2} = [11 12];
    sideplots{3} = [17 18];
    subplot(maxnumplot,6,sideplots{plotcount})
        % background color
        
        % plot all positions
        plot(xpos,ypos,'.','markersize',5,'color',[.8 .8 .8]); hold on
        % plot CP positions 
        plot(xpos(CPvec),ypos(CPvec),'.','markersize',5,'color',CPclr); hold on
        % plot traversal positions
        plot(xpos(pos_a:pos_b),...
              ypos(pos_a:pos_b),'.','markersize',15,'color',trajclr2);
        % plot head direction arrows  
        for ii = pos_a:10:pos_b  
            xdir = 5*cos(hdir(ii));
            ydir = 5*sin(hdir(ii));
            % plot "arrow line"
            plot([xpos(ii)  xpos(ii)+xdir],...
                 [ypos(ii)  ypos(ii)+ydir],'linewidth',2,'color','k');
             % plot point at base of arrow
             scatter(xpos(ii),ypos(ii),20,'k','.');
             % plot point at tip of arrow
             plot(xpos(ii)+xdir,ypos(ii)+ydir,'.','markersize',20,'color','k');
        end
         % formatting
         set(gca,'xticklabel',[],'yticklabel',[],'xtick',[],'ytick',[])
         set(gca,'color',trajclr1)
         axis tight
         axis square
         
     plotcount = plotcount + 1;
    trv = trv + 1;       

end





























