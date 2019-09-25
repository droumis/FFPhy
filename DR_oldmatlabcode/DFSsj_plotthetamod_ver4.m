% Ver4 : Starting 10Feb2014 - Sync codes with everyone

% Theta modulation of cells. 
% Also see sj_plotthetamod
% Gather data like DFSsj_getcellinfo and DFSsj_HPexpt_xcorrmeasures2

clear; %close all;
runscript = 0;
savedata = 0; % save data option - only works if runscript is also on
figopt1 = 0; % Figure Options - Indicidual cells

%savedir = '/data15/gideon/ProcessedData/';
savedir = '/data25/sjadhav/HPExpt/ProcessedData/';
[y, m, d] = datevec(date);

if runscript==1
    %val=1; savefile = [savedir 'HP_thetamod_PFC_alldata_',num2str(m),'-',num2str(d),'-',num2str(y)]; area = 'PFC'; clr = 'b'; % PFC
    %val=2; savefile = [savedir 'HP_thetamod_CA1_alldata_',num2str(m),'-',num2str(d),'-',num2str(y)]; area = 'CA1';  clr = 'r';% CA1
else
    val=1; savefile = [savedir 'HP_thetamod_PFC_alldata_2-12-2014']; area = 'PFC'; clr = 'b'; % PFC - 153 Theta-Phase Locked, 253 Total  
    %val=2; savefile = [savedir 'HP_thetamod_CA1_alldata_2-12-2014']; area = 'CA1';  clr = 'r';% CA1 -  391 Theta-Phase Locked, 493 Total
end

savefig1=0;


% Plot options
plotanimidx =  []; % To pick animals for plotting
plotdays = []; % If you only load data when runscript=0 and savedata=0, then this field will supplant days


% If runscript, run Datafilter and save data
if runscript == 1
    
    %Animal selection
    %-----------------------------------------------------
    %animals = {'HPa','HPb','HPc','Nadal'};
    animals = {'HPa','HPb', 'HPc','Ndl'};
    
    %Filter creation
    %-----------------------------------------------------
    
    % Epoch filter
    % -------------
    %dayfilter = '1:8'; % Shantanu - I am adding day filter to parse out epoch filter
    % Either Only do 1st w-track. 2 or 1 epochs per day
    % Or do Wtr1 and Wtr2, 2 epochs per day
    runepochfilter = 'isequal($environment, ''wtr1'') || isequal($environment, ''wtr2'') || isequal($environment, ''ytr'')'; 
    
    % Cell filter
    % -----------
    switch val
        case 1
            cellfilter = 'strcmp($area, ''PFC'')'; % This includes all, including silent cells. Spike criterion later
        case 2    
            cellfilter = '(strcmp($area, ''CA1'') || strcmp($area, ''iCA1'')) '; % This includes all, including silent cells. Spike criterion later
    end
    % For more exclusive choosing, use $tag
    % Take care of number of spikes while gathering data
    
    % Time filter
    % -----------
    riptetfilter = '(isequal($descrip, ''riptet''))';
    
    % Change on 02-12-2014, for ver4. Change to absvel criterion
    % ------------------------------------------------------------
    timefilter_theta = {{'DFTFsj_getvelpos', '(($absvel >= 5))'},...
        {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',2} };     

    %     timefilter_theta = { {'DFTFsj_getlinstate', '(($state ~= -1) & (abs($linearvel) >= 5))', 6},...
%         {'DFTFsj_getriptimes','($nripples == 0)','tetfilter',riptetfilter,'minthresh',2} };

    % Can also use gethighthetatimes2 -
    
    % EEG Filter
    % ----------
    % Options are sametet, tetfilter, maxvar, file. kk also added "static tet"
    % I am going to use CA1Ref tetrode. Can pass along the file, or use $descrip
    ca1reftetfilter = '( isequal($descrip, ''CA1Ref'') || isequal($descrip, ''reftet'') )';
    % Both the normal and the Gnd files are the same for Ref electrode
    % theta or thetagnd. Make sure to add thetagnd to "iseeg" field
    eegfilter = {'sj_geteegtet', 'theta', 'tetfilter',ca1reftetfilter};
    
    % Iterator
    % --------
    iterator = 'singlecelleeganal';
    
    % Filter creation
    % ----------------
    %modf = createfilter('animal',animals,'days',dayfilter,'epochs',runepochfilter, 'excludetimefilter', timefilter_theta, 'cells',...
    %    cellfilter,'eegtetrodes', eegfilter, 'iterator', iterator);
    modf = createfilter('animal',animals,'epochs',runepochfilter, 'excludetimefilter', timefilter_theta, 'cells',...
        cellfilter,'eegtetrodes', eegfilter, 'iterator', iterator);
    
    disp('Done Filter Creation');
    
    % Set analysis function
    % ----------------------
    
    modf = setfilterfunction(modf,'DFAsj_plotthetamod4',{'spikes','theta'}); % Corresponding function is sj_plotthetamod
    %modf = setfilterfunction(modf,'DFAsj_plotthetamod4',{'spikes','theta'},'nbins',36); % Options can be "threshspks" or "nbins"
    
    % Run analysis
    % ------------
    modf = runfilter(modf);
    disp('Finished running filter script');
    %--------------------- Finished Filter Function Run -------------------
    
    if savedata == 1
        clear figopt1 runscript plotdays plotanimidx savedata
        save(savefile);
    end
    
else
    
    load(savefile);
    
end % end runscript

if ~exist('savedata')
    return
end


% -------------------------  Filter Format Done -------------------------



% ----------------------------------
% Whether to gather data or to load previously gathered data
% --------------------------------------------------------------------
gatherdata = 0; savegatherdata = 0;
[y, m, d] = datevec(date);

switch val
    case 1
        if gatherdata
            gatherdatafile = [savedir 'HP_thetamod_PFC_alldata_Nspk50_gather_',num2str(m),'-',num2str(d),'-',num2str(y)]; % PFC cells to Hipp theta
        else
            gatherdatafile = [savedir 'HP_thetamod_PFC_alldata_Nspk50_gather_2-12-2014']; % PFC cells to Hipp theta
        end
    case 2
        if gatherdata
            gatherdatafile = [savedir 'HP_thetamod_CA1_alldata_Nspk50_gather_',num2str(m),'-',num2str(d),'-',num2str(y)]; % CA1 cells to Hipp theta
        else
            gatherdatafile = [savedir 'HP_thetamod_CA1_alldata_Nspk50_gather_2-12-2014']; % CA1 cells to Hipp theta
            
        end
end


if gatherdata
    
    % Parameters if any
    % -----------------
    nbins = 50;
    % Circ Stats Box Von Mises pdf uses a default of 100 angles/nbin
    % -------------------------------------------------------------
    
    cnt=0; % Count how man cells will be kept based on nspikes in output
    allanimindex=[]; alldata=[]; all_Nspk=[];
    
    
    for an = 1:length(modf)
        for i=1:length(modf(an).output{1})
            % Check for empty output
            if (modf(an).output{1}(i).Nspikes > 0)
                cnt=cnt+1;
                anim_index{an}(cnt,:) = modf(an).output{1}(i).index;
                % Only indexes
                animindex=[an modf(an).output{1}(i).index]; % Put animal index in front
                allanimindex = [allanimindex; animindex]; % Collect all Anim Day Epoch Tet Cell Index
                % Data
                alldata{cnt} = modf(an).output{1}(i).sph; % Only get spike phases. Compute evrything else after combining epochs
                all_Nspk(cnt) = modf(an).output{1}(i).Nspikes;
            end
        end
        
    end
    
    % Consolidate single cells across epochs. 2 methods: see DFSsj_getcellinfo
    % and DFSsj_xcorrmeasures2
    % Here, I will use combined animal-index. Prevents having to loop over animals
    % ----------------------------------------------------------------------------
    allthetamod = struct;
    dummyindex=allanimindex;  % all anim-day-epoch-tet-cell indices
    cntcells=0;
    
    for i=1:size(alldata,2)
        animdaytetcell=allanimindex(i,[1 2 4 5]);
        ind=[];
        while rowfind(animdaytetcell,dummyindex(:,[1 2 4 5]))~=0          % collect all rows (epochs)
            ind = [ind rowfind(animdaytetcell,dummyindex(:,[1 2 4 5]))];        % finds the first matching row
            dummyindex(rowfind(animdaytetcell,dummyindex(:,[1 2 4 5])),:)=[0 0 0 0 0]; % after adding index, remove the corresponding row
            % so you could find the next one
        end
        
        % Gather everything for the current cell across epochs
        currsph=[]; currNspk=0;
        for r=ind
            currNspk = currNspk + all_Nspk(r);
            currsph = [currsph; alldata{r}];
        end
        
        if currNspk >= 50 % based on Siapas 2005; maybe it should be 100
            cntcells = cntcells + 1;
            % DONT SHIFT ANY DAYS 
            % For Ndl-GIdeon;s animal, shift days
%             if animdaytetcell(1)==4
%                 animdaytetcell(2)=animdaytetcell(2)-7; % Day starts from no. 8
%             end
            allthetamod_idx(cntcells,:)=animdaytetcell;
            allthetamod(cntcells).index=animdaytetcell;
            allthetamod(cntcells).sph=currsph; % Remember, this is seperated across epochs. Can take mean
            allthetamod(cntcells).Nspk=currNspk;
        end
    end
    
    
    % Calculations - can incorporate in loop above as well
    % ------------
    for i=1:cntcells
        sph = allthetamod(i).sph;
        % Rayleigh test and Modulation:
        stats = rayleigh_test(sph); % stats.p and stats.Z, and stats.n
        [m, ph] = modulation(sph);
        phdeg = ph*(180/pi);
        % Von Mises Distribution - From Circular Stats toolbox
        [thetahat, kappa] = circ_vmpar(sph); % Better to give raw data. Can also give binned data.
        thetahat_deg = thetahat*(180/pi);
        [prayl, zrayl] = circ_rtest(sph); % Rayleigh test for non-uniformity of circular data
        % Von Mises Fit - Use nbins defined above
        bins = -pi:(2*pi/nbins):pi;
        count = histc(sph, bins);
        % Make Von Mises Fit
        alpha = linspace(-pi, pi, 50)';
        [pdf] = circ_vmpdf(alpha,thetahat,kappa);
        % Another way of doing mean phase
        A=mean(exp(j*sph)); % problem with i - clashes with variable. so use j.
        meancos=real(A);
        meansin=imag(A);
        meanphase=atan2(meansin,meancos); % Exactly the same as von mises fit
        
        % Save
        allthetamod(i).thetahist = count; % Theta histogram plot
        allthetamod(i).thetahistnorm = count./max(count); % Normalized - Theta histogram plot
        allthetamod(i).thetaper = (count./sum(count))*100; % Histogram in units of percentage of spikes
        allthetamod(i).stats = stats;
        allthetamod(i).modln = m;
        allthetamod(i).phdeg = phdeg;
        allthetamod(i).kappa = kappa;
        allthetamod(i).thetahat = thetahat; % From von Mises fit - use this
        allthetamod(i).thetahat_deg = thetahat_deg; % From von Mises fit - use this
        allthetamod(i).prayl = prayl;
        allthetamod(i).zrayl = zrayl;
        allthetamod(i).alpha = alpha;
        allthetamod(i).vmpdf = pdf;
        allthetamod(i).meanphase = meanphase;
        allthetamod(i).anim = allthetamod(i).index(1); allanim(i) = allthetamod(i).index(1);
        allthetamod(i).days = allthetamod(i).index(2); alldays(i) = allthetamod(i).index(2);
    end
    
    % Save
    % -----
    if savegatherdata == 1
        save(gatherdatafile);
    end
    
else % gatherdata=0
    
    load(gatherdatafile);
    
end % end gather data

%add to "super" variable
% super{reg}=modf;


% ------------------------------
% Plotting for individual cells
% ------------------------------


figdir = '/data25/sjadhav/HPExpt/Figures/ThetaMod/Egs/';
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
forppr=0;
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

figopt1=0;
if (figopt1)
    for i=1:cntcells
        
        curridx = allthetamod(i).index;
        switch curridx(1)
            case 1
                prefix = 'HPa';
            case 2
                prefix = 'HPb';
            case 3
                prefix = 'HPc';
        end
        day = curridx(2); tet = curridx(3); cell = curridx(4);
        % To control plotting
        %if curridx(1)==3 && sig_shuf==1  
        if curridx(1)==2 && curridx(2)==2 && curridx(3)==12 && curridx(4)==1
            
            
            sph = allthetamod(i).sph;
            Nspk = allthetamod(i).Nspk;
            thetahist = allthetamod(i).thetahist; % Theta histogram plot
            thetahistnorm = allthetamod(i).thetahistnorm; % Normalized - Theta histogram plot
            thetaper = allthetamod(i).thetaper; % Histogram in units of percentage of spikes
            stats = allthetamod(i).stats;
            m = allthetamod(i).modln;
            phdeg = allthetamod(i).phdeg;
            kappa = allthetamod(i).kappa;
            thetahat_deg = allthetamod(i).thetahat_deg; % From von Mises fit - use this
            prayl = allthetamod(i).prayl;
            zrayl = allthetamod(i).zrayl;
            alpha = allthetamod(i).alpha;
            pdf = allthetamod(i).vmpdf;
            meanphase = allthetamod(i).meanphase;
            
            countper = thetaper;
            bins = -pi:(2*pi/nbins):pi;
            ph = phdeg*(pi/180)
            
            figure; hold on;
            set(gcf,'Position',[970 86 925 1000]); % Vertical figure at right edge of screen
            %set(gcf,'Position',[750 90 1145 1000]); % For 2 X 2 plot
            %redimscreen_2versubplots;
            % subplot(2,2,1); hold on; % Raster
            % Cant do raster without getting data for each individual theta cycle
            
            %         subplot(3,1,1); hold on; % Hist with Nspikes
            %         out = bar(bins, thetahist, 'hist'); set(gca, 'XTick', [-pi, -pi/2, 0, pi/2, pi]);
            %         set(gca,'XLim',[-pi pi]); ylabel('NSpikes');
            %         set(out,'FaceColor','r'); set(out,'EdgeColor','r');
            %         %pdf = pdf.*(max(count)/max(pdf));
            %         % Instead of maximum - match values at a bin, maybe close to peak
            %         binnum = lookup(thetahat,alpha);
            %         pdf = pdf.*(count(binnum)/pdf(binnum));
            %         plot(alpha,pdf,'k','LineWidth',3,'Color','k');
            
            
            %subplot(2,1,1); hold on; % Hist with percentage of spikes
            out = bar(bins, thetaper, 'hist'); set(gca, 'XTick', [-pi, -pi/2, 0, pi/2, pi]);
            set(gca,'XLim',[-pi pi]); ylabel('NSpikes');
            set(gca,'XLim',[-pi pi]); ylabel('% of Spikes');
            set(out,'FaceColor','r'); set(out,'EdgeColor','r');
            binnum = lookup(thetahat,alpha);
            pdf = pdf.*(countper(binnum)/pdf(binnum));
            plot(alpha,pdf,'k','LineWidth',3,'Color','k');
            title(sprintf('%d %d %d %d, Kappa%f, prefang%g, pval%f', curridx, kappa, thetahat_deg, prayl));
            %title(sprintf('Nspikes %d', totalspks));
            
%             subplot(2,1,2); hold on; % Polar plot
%             [t,r] = rose(sph);
%             polar(t,r,'r'); hold on;
%             % The peak phase angle
%             lims = get(gca,'XLim');
%             radius = lims(2);
%             xx = radius .* cos(ph); yy = radius .* sin(ph);
%             line([0 xx], [0 yy],'LineWidth',4,'Color','k');
%             title(sprintf('%s Day %d Tet %d Cell %d', prefix, day, tet, cell),'FontSize',tfont,'Fontweight','normal');

            
            figfile = [figdir,area,'EgThetamod_',num2str(i)];
            keyboard;
            
        end
        
        %print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
        
        
        
    end % end cntcells
    
end % end if figopt


% ------------------
% Population Figures
% ------------------

% Population Data
% ----------------------------

allsigphases = []; cntsig = 0;
allkappas = []; allsigkappas = [];
allZ = []; allsigZ = []; allm=[];
allsph =[]; allthetahist=[]; allthetahistnorm=[];

days = unique(alldays);
anim = unique(allanim);
% Force days to got from 1-10. For Ndl, you will push days by 7
days = 1:10;
ncells_days = zeros(length(days),1);
ncells_days_sig = zeros(length(days),1);

for i = 1:length(allthetamod)
    allkappas(i) = allthetamod(i).kappa;
    allZ(i) = allthetamod(i).zrayl;
    
    curranim = allthetamod(i).anim;
    currday = allthetamod(i).days;
    if curranim==4
        currday = currday-7; % For Ndl, days start from 8
    end
    ncells_days(currday) = ncells_days(currday)+1;
    
    if (allthetamod(i).prayl < 0.05)
        cntsig = cntsig+1;
        ncells_days_sig(currday) = ncells_days_sig(currday)+1;
        allsigphases(cntsig) = allthetamod(i).meanphase; % Can use mean phase (atan) or thetahat_deg (from von Mises distribution)
        allsigkappas(cntsig) = allthetamod(i).kappa;
        allsigZ(cntsig) = allthetamod(i).zrayl;
        allm(cntsig) = allthetamod(i).modln;
        allsph = [allsph; allthetamod(i).sph]; % All spikes pooled
        
        allthetahist(cntsig,:) = allthetamod(i).thetahist;
        allthetaper(cntsig,:) = allthetamod(i).thetaper;
        allthetahistnorm(cntsig,:) = allthetamod(i).thetahistnorm;
               
    end
end


forppr = 0; 
% If yes, everything set to redimscreen_figforppr1
% If not, everything set to redimscreen_figforppt1

figdir = '/data25/sjadhav/HPExpt/Figures/31Oct/';
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


if 1
    % 1) Histogram of mean phases
    % ----------------------------
    
    figure; hold on;
    if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end
    
    nbins=24; % For population plots, reduce nbins
    binsp = -pi:(2*pi/nbins):pi;
    
    N = histc(allsigphases,binsp);
    N=N(1:(end-1));
    bins_plot = binsp(1:(end-1));
    bins_plot = bins_plot + (binsp(2)-binsp(1))/2;
    h = bar([bins_plot bins_plot+2*pi],[N , N],'histc');
    set(h(1),'facecolor',clr)
    set(h(1),'edgecolor',clr)
    xlabel(['Phase'],'FontSize',xfont,'Fontweight','normal');
    ylabel(['No. of cells'],'FontSize',yfont,'Fontweight','normal');

    title(sprintf('Mean phases of sig. locked units: %d',cntsig),'FontSize',tfont,'FontWeight','normal')
    axis tight
    hold on
    plot([pi,pi],[0 999],'k--','LineWidth',1.5)
    plot([-pi,-pi],[0 999],'k--','LineWidth',1.5)
    plot([3*pi,3*pi],[0 999],'k--','LineWidth',1.5)
    plot([0,0],[0 999],'k:','LineWidth',1.5)
    plot([2*pi,2*pi],[0 999],'k:','LineWidth',1.5)
    
    %set(gca, 'XTick', [-pi:pi:3*pi], 'XTickLabel',num2str([-180,0,180,0,-180]'));
    
    a = num2str([-180,0,180,0,180]'); a(3,:) = '+-pi';
    set(gca, 'XTick', [-pi:pi:3*pi], 'XTickLabel',a);
    set(gca,'XLim',[-pi 3*pi]);
    %set(gca, 'XTick', [-pi], 'XTickLabel',sprintf('%s','-pi'));
    %set(gca, 'XTick', [0], 'XTickLabel',sprintf('%s','0'));

    if savefig1==1,
        figfile = [figdir,area,'_Thetamod_MeanPhases'];
        print('-depsc2', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
end


if 1
    % 2)plot phase histogram of aggregate spikes, sig units
     % ----------------------------------------------------
    
    norm = 1;
    
    figure; hold on;
    if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end
    phasehist=histc(allsph,bins);
    phasehist = phasehist(1:(end-1));
    phasehist_norm = phasehist/(sum(phasehist))*100; % percentage of spikes
    bins_plot = bins(1:(end-1));
    bins_plot = bins_plot + (bins(2)-bins(1))/2;
    if norm == 1
        h = bar([bins_plot bins_plot+2*pi],[phasehist_norm; phasehist_norm],'histc');
        set(h(1),'facecolor',clr)
        set(h(1),'edgecolor',clr)
        axis tight
        %ylim([0 .1])
        set(gca,'YLim',[0 max(phasehist_norm)+0.5])
        ylabel(['%tage of spikes'],'FontSize',yfont,'Fontweight','normal');
    else
        h = bar([bins_plot bins_plot+2*pi],[phasehist; phasehist],'histc');
        set(h(1),'facecolor',clr)
        set(h(1),'edgecolor',clr)
        axis tight
        ylim([0 max(phasehist)+1500])
        ylabel(['No. of cells'],'FontSize',yfont,'Fontweight','normal');

    end
    xlabel(['Phase'],'FontSize',xfont,'Fontweight','normal');
    title(sprintf('Phase hist of all spikes, sig units: %d',cntsig),'FontSize',tfont,'FontWeight','normal');
    hold on
    plot([pi,pi],[0 99999],'k--','LineWidth',1.5)
    plot([-pi,-pi],[0 99999],'k--','LineWidth',1.5)
    plot([3*pi,3*pi],[0 99999],'k--','LineWidth',1.5)
    plot([0,0],[0 99999],'k:','LineWidth',1.5)
    plot([2*pi,2*pi],[0 99999],'k:','LineWidth',1.5)
    
    %set(gca, 'XTick', [-pi:pi:3*pi], 'XTickLabel',num2str([-180,0,180,0,-180]'));
    a = num2str([-180,0,180,0,180]'); a(3,:) = '+-pi';
    set(gca, 'XTick', [-pi:pi:3*pi], 'XTickLabel',a);
    set(gca,'XLim',[-pi 3*pi]);
    
    figfile = [figdir,area,'_Thetamod_PhaseHistAllSpks'];
    if savefig1==1,
        
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end

    
end



if 1
   % 3) Matrix of normalized histograms aligned by peak phase 
    figure; hold on;
    if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end
   
   [sortedph, order] = sort(allsigphases,2,'descend');
   sort_histnorm = allthetahistnorm(order,:);
   sort_histper = allthetaper(order,:);
   
   smsort=[];    bins_plot = bins(1:(end-1));
   % nstd=1: gaussian of length 4. nstd = 2: gaussian of length 7, nstd=3: gaussian of length 10.
    nstd = 3; g1 = gaussian(nstd, 3*nstd+1);

   for n =1:length(order),
       curr = sort_histnorm(n,1:50); % Last bin should be skipped
       curr = smoothvect(curr,g1);
       smsort(n,:) = curr(2:end-1);
       smsort(n,:) = smsort(n,:)./max(smsort(n,:)); % Renormalize
   end
   bins_plot = bins_plot(2:end-1);
   
   imagesc(bins_plot,1:n,smsort); colorbar;
   set(gca,'XLim',[-pi pi]); set(gca,'YLim',[0 n]);
   a = num2str([-180,0,180]');
   set(gca, 'XTick', [-pi:pi:pi], 'XTickLabel',a);
   xlabel(['Phase'],'FontSize',xfont,'Fontweight','normal'); 
   ylabel(['Cell no'],'FontSize',yfont,'Fontweight','normal'); 
   title(sprintf('Phase-locked units aligned by Pref Phase: %d',cntsig),'FontSize',tfont,'FontWeight','normal');
   
   figfile = [figdir,area,'_Thetamod_MatrixAlignPrefPhase'];
   if savefig1==1,       
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
   
end


if 1
   % 4) Matrix of normalized histograms aligned by concentration parameter
    figure; hold on;
    if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end
   
   [sortedk, order] = sort(allsigkappas);
   sort_histnorm = allthetahistnorm(order,:);
   sort_histper = allthetaper(order,:);
   
   smsort=[];    bins_plot = bins(1:(end-1));
   % nstd=1: gaussian of length 4. nstd = 2: gaussian of length 7, nstd=3: gaussian of length 10.
    nstd = 3; g1 = gaussian(nstd, 3*nstd+1);

   for n =1:length(order),
       curr = sort_histnorm(n,1:50); % Last bin should be skipped
       curr = smoothvect(curr,g1);
       smsort(n,:) = curr(2:end-1);
       smsort(n,:) = smsort(n,:)./max(smsort(n,:)); % Renormalize
   end
   bins_plot = bins_plot(2:end-1);
   
   imagesc(bins_plot,1:n,smsort); colorbar;
   set(gca,'XLim',[-pi pi]); set(gca,'YLim',[0 n]);
   a = num2str([-180,0,180]');
   set(gca, 'XTick', [-pi:pi:pi], 'XTickLabel',a);
   xlabel(['Phase'],'FontSize',xfont,'Fontweight','normal'); 
   ylabel(['Cell no'],'FontSize',yfont,'Fontweight','normal'); 
   title(sprintf('Phase-locked units aligned by Conc Parm : %d',cntsig),'FontSize',tfont,'FontWeight','normal');
   
   figfile = [figdir,area,'_Thetamod_MatrixAlignModlnStrength'];
   if savefig1==1,
        
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
   end
    
end


if 1
    % 6) Plot distribution of Kappas
    % ------------------------------------
    figure; hold on;
    if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end
    currbins = [0:0.2:max(allsigkappas)]
    N = histc(allsigkappas,currbins);
    h = bar(currbins,N,'histc');
    set(h(1),'facecolor',clr);
    set(h(1),'edgecolor',clr);
    xlabel(['Conc parameter (Kappa)'],'FontSize',xfont,'Fontweight','normal');
    ylabel(['No. of cells'],'FontSize',yfont,'Fontweight','normal');
    title(sprintf('Conc par (kappa) of sig. locked units: %d',cntsig),'FontSize',tfont,'FontWeight','normal')
    set(gca,'XLim',[-0.15 max(allsigkappas)+0.1])
    if savefig1==1,
        figfile = [figdir,area,'_Thetamod_KappaDistr'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
   end
end


if 1
    % 7) Plot distribution of Rayleigh Z
    % ------------------------------------
    figure; hold on;
    if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end
    
    logZ = log(allZ);
    logsigZ = log(allsigZ);
    currbins = min(logZ):0.8:max(logZ);
    N = histc(logZ,currbins); Nsig = histc(logsigZ,currbins);
    h = bar(currbins,N,'histc');
    set(h(1),'facecolor',clr);
    set(h(1),'edgecolor',clr);
    h2 = bar(currbins,Nsig,'histc');
    set(h2(1),'facecolor','k');
    set(h2(1),'edgecolor','k');
    xlabel(['log(Rayleigh Z)'],'FontSize',xfont,'Fontweight','normal');
    ylabel(['No. of cells'],'FontSize',yfont,'Fontweight','normal');
    title(sprintf('Rayleigh Z distribution'),'FontSize',tfont,'FontWeight','normal')
    plot([1 1],[0 max(N)],'k--','LineWidth',1.5);
    %set(gca,'XLim',[-0.15 max(allsigkappas)+0.1])
    if savefig1==1,
        figfile = [figdir,area,'_Thetamod_RayleighZDistr'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
   end
    
end




if 1
    % 7) plot normalized histogram sig units, SEM
     % ------------------------------------------
     
    figure; hold on;
    if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end
  
    phasehist_mean=mean(allthetahist,1);
    phasehist_sem=std(allthetahist,1)/sqrt(size(allthetahist,1));
    phasehist_mean = phasehist_mean(1:(end-1));
    bins_plot = bins(1:(end-1));
    bins_plot = bins_plot + (bins(2)-bins(1))/2;
    
    h = bar([bins_plot bins_plot+2*pi],[phasehist_mean phasehist_mean],'histc');
    set(h(1),'facecolor',clr)
    set(h(1),'edgecolor',clr)
    axis tight
    if strcmp(area,'PFC')
        ylim([0 90])
    else
        ylim([0 200]);
    end
    hold on
    
    % plot sem bars
    for jj=1:length(bins_plot)
        plot([bins_plot(jj),bins_plot(jj)],[phasehist_mean(jj)-phasehist_sem(jj) phasehist_mean(jj)+phasehist_sem(jj)],'k','LineWidth',1.5)
        plot([bins_plot(jj)+2*pi,bins_plot(jj)+2*pi],[phasehist_mean(jj)-phasehist_sem(jj) phasehist_mean(jj)+phasehist_sem(jj)],'k','LineWidth',1.5)
    end
    
    titlestring=sprintf('Phase hist of all spikes, sig units: %d',cntsig);
    title(titlestring,'FontSize',14,'FontWeight','bold')
    hold on
    plot([pi,pi],[0 99999],'k--','LineWidth',1.5)
    plot([-pi,-pi],[0 99999],'k--','LineWidth',1.5)
    plot([3*pi,3*pi],[0 99999],'k--','LineWidth',1.5)
    plot([0,0],[0 99999],'k:','LineWidth',1.5)
    plot([2*pi,2*pi],[0 99999],'k:','LineWidth',1.5)
    %set(gca, 'XTick', [-pi:pi:3*pi], 'XTickLabel',num2str([-180,0,180,0,-180]'));
    a = num2str([-180,0,180,0,180]'); a(3,:) = '+-pi';
    set(gca, 'XTick', [-pi:pi:3*pi], 'XTickLabel',a);
    set(gca,'XLim',[-pi 3*pi]);
    
     if savefig1==1,
        figfile = [figdir,area,'_Thetamod_PhaseHistAllSpksWithErr'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
end



if 1
    % 8) plot distribution of modulation depths
     % ------------------------------------
    figure; hold on;
    if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end

    hist(allm,10,clr);
    title(sprintf('Distribution of modulation depths, nunits: %d',cntsig))
    xlabel(['Modulation Depth'],'FontSize',xfont,'Fontweight','normal');
    ylabel(['No. of cells'],'FontSize',yfont,'Fontweight','normal');
    
    if savefig1==1,
        figfile = [figdir,area,'_Thetamod_ModlnDepthDistr'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
   
end



if 1
    % 9) No of sig phase locked cells over days: %tage and number
     % ------------------------------------
    figure; hold on;
    if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end

    persig_days = 100*ncells_days_sig./ncells_days;
    plot(persig_days,[clr 'o'],'MarkerSize',18,'LineWidth',2);
    title(sprintf('%tage of sig phase locked units'));
    xlabel(['Day'],'FontSize',xfont,'Fontweight','normal');
    ylabel(['Percentage of cells'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'YLim',[0 max(persig_days)+5]);
    
    if savefig1==1,
        figfile = [figdir,area,'_Thetamod_PerSigDays'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
   
    
    figure; hold on;
    if forppr==1, redimscreen_figforppr1; else redimscreen_figforppt1; end

    Nsig_days = ncells_days_sig;
    plot(Nsig_days,[clr 'o'],'MarkerSize',18,'LineWidth',2);
    title(sprintf('No. of sig phase locked units'));
    xlabel(['Day'],'FontSize',xfont,'Fontweight','normal');
    ylabel(['Number of cells'],'FontSize',yfont,'Fontweight','normal');
    set(gca,'YLim',[0 max(Nsig_days)+2]);

    if savefig1==1,
        figfile = [figdir,area,'_Thetamod_NSigDays'];
        print('-dpdf', figfile); print('-djpeg', figfile); saveas(gcf,figfile,'fig');
    end
    
    
end























% 
% if 0
%     % plot individual phase histogram of all units
%     
%     norm = 1;
%     
%     figure
%     titlestring=sprintf('%s %s phase hist of individual units // %s',animals{1},region,referencestring);
%     title(titlestring,'FontSize',14,'FontWeight','bold')
%     counter=1;
%     for k=1:length(caf.celloutput)
%         if counter==81
%             counter=1;
%             figure
%             titlestring=sprintf('%s %s phase hist of individual units // %s',animals{1},region,referencestring);
%             title(titlestring,'FontSize',14,'FontWeight','bold')
%         end
%         subplot(8,10,counter)
%         bins_plot = caf.celloutput(k).bins(1:(end-1));
%         bins_plot = bins_plot + (bins(2)-bins(1))/2;
%         phasehist = caf.celloutput(k).phasehist(1:(end-1));
%         phasehist_norm = phasehist/sum(phasehist);
%         if norm == 1
%             if size(phasehist_norm,1) < size(phasehist_norm,2)
%                 phasehist_norm = phasehist_norm';
%             end
%             %plot
%             h = bar([bins_plot bins_plot+2*pi],[phasehist_norm ; phasehist_norm],'histc');
%             title(num2str(caf.celloutput(k).index))
%             axis tight
%             ylim([0 .2])
%         else
%             if size(phasehist,1) < size(phasehist,2)
%                 phasehist = phasehist';
%             end
%             %plot
%             h = bar([bins_plot bins_plot+2*pi],[phasehist ; phasehist],'histc');
%             title(num2str(caf.celloutput(k).index),'FontSize',12,'FontWeight','bold')
%             axis tight
%             ylim([0 250])
%         end
%         
%         set(h(1),'facecolor',clr)
%         set(h(1),'edgecolor',clr)
%         
%         % plot guide lines
%         hold on
%         plot([pi,pi],[0 9999],'k--','LineWidth',1.5)
%         plot([-pi,-pi],[0 9999],'k--','LineWidth',1.5)
%         plot([3*pi,3*pi],[0 9999],'k--','LineWidth',1.5)
%         plot([0,0],[0 9999],'k:','LineWidth',1.5)
%         plot([2*pi,2*pi],[0 9999],'k:','LineWidth',1.5)
%         
%         counter=counter+1;
%     end
% end
% 



% % bar
% count = histc(allspikephases, bins);
% out = bar(bins, count, 'hist');
% set(out,'facecolor','k')
% title('aggregate theta modulation');
%
% % lineplot
% dischargeprob=count./sum(count);
% plot(bins(1:(end-1)),dischargeprob(1:(end-1)),'k','LineWidth',2);
%
% [m ph] = modulation(allspikephases);






