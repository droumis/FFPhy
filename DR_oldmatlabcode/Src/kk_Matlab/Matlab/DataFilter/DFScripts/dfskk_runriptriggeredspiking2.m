
% Unlike dfskk_runriptriggeredspiking, this code preserves ripple-by-ripple
% firing (abandons use of xcorr).

runscript = 0;

% pre parameters
animals = {'Egypt'};
window = [.5 .5]          % size of psth window (in sec)
binsize = .005         % size of bins (in sec)
minthresh_rip = 3 ;     % how big your ripples are
exclusion_dur = 0.5;   % seconds within which consecutive ripples are eliminated
exclusion_nrip = 2;    % this should be the same as your nripples you specify when you call kk_getriptimes 

% analyses
postprocessing = 0;
unnorm_unit_psth = 0;
    smooth_flag = 1;
    baseline = [-0.5 -.4];
    ripple_time = [0 .1];
    numboots = 3000;
        
regional_aggregate = 1;
    smooth_flag = 1;
    negative_mod = 0;  % also plot units with negative modulation

    
smoothed_aggregate_latency = 0;     % plots smoothed aggregates + calculates latencies
    calculate_p_flag = 0;           % if calculate, set binsize very small!
        numperms = 1000; 

% post parameters
sig_level = 0.05;
smoothing_length = 10;   % std of smoothing gaussian (in ms)
regions = [1 3 2];
pairs = [1 3 ; 1 2 ; 3 2];
    latencies = nan(size(pairs));
min_numspikes = 5;    % don't process unit if too few total spikes in window around ripples


    

if runscript
    

% Animal Selection


% Epoch Filter
%epochfilter = '(isequal($type, ''run''))';
epochfilter = '(isequal($type, ''run'') || isequal($type, ''sleep''))';


% Time Filter
timefilter = { {'kk_get2dstate', '$velocity < 4'}  ;  ...
               {'kk_getriptimes', '($nripples >= 2)', [], 'cellfilter', '(isequal($area, ''CA1''))','minthresh',minthresh_rip,'exclusion',exclusion_dur,exclusion_nrip}};   % IMPORTANT HERE (ripples chosen here)
               
% 
ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && isequal($type, ''principal''))';    %% ($meanrate < 3))';
ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100) && isequal($type, ''principal''))';   
ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && isequal($type, ''principal''))';

%ca1cellfilter = '(isequal($area, ''CA1'') && ($numspikes > 100) && ($meanrate < 3))';
%ca2cellfilter = '(isequal($area, ''CA2'') && ($numspikes > 100) && ($meanrate < 3))';
%ca3cellfilter = '(isequal($area, ''CA3'') && ($numspikes > 100) && ($meanrate < 3))';

% Iterator
iterator = 'singlecellanal';

% Filter Creation
ca1f = createfilter('animal', animals, 'epochs', epochfilter, 'cells', ca1cellfilter, 'excludetime', timefilter, 'iterator', iterator);
ca2f = createfilter('animal', animals, 'epochs', epochfilter, 'cells', ca2cellfilter, 'excludetime', timefilter, 'iterator', iterator);
ca3f = createfilter('animal', animals, 'epochs', epochfilter, 'cells', ca3cellfilter, 'excludetime', timefilter, 'iterator', iterator);

% Set Analysis Function
ca1f = setfilterfunction(ca1f, 'dfakk_getriptriggeredspiking2', {'spikes','ripples','tetinfo'},'tetfilter','(isequal($area, ''CA1''))','window',window,'binsize',binsize);
ca2f = setfilterfunction(ca2f, 'dfakk_getriptriggeredspiking2', {'spikes','ripples','tetinfo'},'tetfilter','(isequal($area, ''CA1''))','window',window,'binsize',binsize);
ca3f = setfilterfunction(ca3f, 'dfakk_getriptriggeredspiking2', {'spikes','ripples','tetinfo'},'tetfilter','(isequal($area, ''CA1''))','window',window,'binsize',binsize);


% Run Analysis
ca1f = runfilter(ca1f);
ca2f = runfilter(ca2f);
ca3f = runfilter(ca3f);

end





% Post processing.
% Consolidate cells across epochs, within a day
    % collect all cell indices, ignoring specific epochs, into daytetcell

if postprocessing

for reg=regions       % iterate over regions
    if reg==1
        region='CA1';
        clr=[0 0 0];
        f=ca1f;
    elseif reg==2
        region='CA2';
        clr=[0 1 0];
        f=ca2f;
    else
        region='CA3';
        clr=[1 0 0];
        f=ca3f;
    end
    
indices = [];
for i=1:length(f.output{1})  % iterate over epochs
    indices = [indices ; f.output{1}(i).index];
end
daytetcell = unique(indices(:,[1 3 4]),'rows');

f.celloutput = struct;

for ind=1:size(daytetcell,1)
    f.celloutput(ind).daytetcell = daytetcell(ind,:);
    f.celloutput(ind).psth = [];
    f.celloutput(ind).nospikes = 0;
    f.celloutput(ind).noripples = 0;
    f.celloutput(ind).noepochs = 0;
    % consolidate over all epochs
    for c=1:length(f.output{1})
        if ~isempty(f.output{1}(c).time)   
            % if day tet cell matches
            if rowfind(daytetcell(ind,:),f.output{1}(c).index([1 3 4]))
                f.celloutput(ind).psth = [f.celloutput(ind).psth ; f.output{1}(c).psth];
                f.celloutput(ind).nospikes = f.celloutput(ind).nospikes + f.output{1}(c).nospikes;
                f.celloutput(ind).noripples = f.celloutput(ind).noripples + f.output{1}(c).noripples;
                f.celloutput(ind).noepochs = f.celloutput(ind).noepochs + 1;
                f.celloutput(ind).time = f.output{1}(c).time;   % (all time vectors are the same..)
            end
        end
    end
    
    f.celloutput(ind).psthsum = sum(f.celloutput(ind).psth,1);
    
    % test direction and significance of modulation using bootstraps
      % get baseline data
      if isempty(f.celloutput(ind).psthsum)
        f.celloutput(ind).pvalue = [];
        f.celloutput(ind).modulation = [];          
      else
      
      bl_startind = lookup(baseline(1),f.celloutput(ind).time);
      bl_endind = lookup(baseline(2),f.celloutput(ind).time);
      rip_startind = lookup(ripple_time(1),f.celloutput(ind).time);
      rip_endind = lookup(ripple_time(2),f.celloutput(ind).time);
      ripple_data = f.celloutput(ind).psth(:,rip_startind:rip_endind);
      baseline_data = f.celloutput(ind).psth(:,bl_startind:bl_endind);  % rows: trials (ripples), columns: bins
      % bootstrap to get dist'n of differences in spikerate
      baseline_spikesperbin = sum(baseline_data,2)/size(baseline_data,2);
      ripple_spikesperbin = sum(ripple_data,2)/size(ripple_data,2);
      N = size(baseline_spikesperbin,1);
      M = size(ripple_spikesperbin,1);
      for n = 1:numboots
          boot_diffspikesperbin(n) = mean(ripple_spikesperbin(ceil(M*rand(M,1)))) - ...
                                     mean(baseline_spikesperbin(ceil(N*rand(N,1))));
      end
      upperbound = prctile(boot_diffspikesperbin,95);
      lowerbound = prctile(boot_diffspikesperbin,5);
      
      pvalue = 0;
      modulation = 0;
      
      % plot test graph
      if 0
          figure
          hold on
          H = hist(boot_diffspikesperbin,100)
          hist(boot_diffspikesperbin,100)
          plot([upperbound upperbound],[0 max(H)+20],'k','LineWidth',3)
          plot([lowerbound lowerbound],[0 max(H)+20],'k','LineWidth',3)
          axis tight
      end
      
      % positive modulation
      if lowerbound > 0
          modulation = 1;
          textcolor = [1 .1 0];
          bgcolor = [250 210 210] / 255;
          pvalue = lookup(0,sort(boot_diffspikesperbin))/numboots;
      end
      % negative modulation
      if upperbound < 0
          modulation = -1;
          textcolor = [.05 .5 1];
          bgcolor = [176 224 230]/255;
          pvalue = 1-lookup(0,sort(boot_diffspikesperbin))/numboots;
      end

      f.celloutput(ind).pvalue = pvalue;
      f.celloutput(ind).modulation = modulation;
    
      end
end
    
    super{reg} = f;
  
end
end


% Plot unnormalized raw unit PSTH.
    
    if unnorm_unit_psth
        
        for reg=regions
            if reg==1
                region='CA1';
                clr=[0 0 0];
            elseif reg==2
                region='CA2';
                clr=[0 .8 0];
            else
                region='CA3';
                clr=[1 0 0];
            end
            figure
            counter = 1;

            for c=1:length(super{reg}.celloutput)
                
                if counter > 60
                    
                    %title of previous figure
                    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
                    titlestring=sprintf('%s %s ripple triggered spiking',animals{1},region);
                    title(titlestring)
                    text(0.5, .99,titlestring,'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',18,'FontWeight','bold')
                    
                    figure
                    counter = 1;
                    titlestring=sprintf('%s %s ripple triggered spiking',animals{1},region);
                    title(titlestring)
                end
                subplot(6,10,counter)
       
                % skip plotting if not enough spikes
                numpsthspikes = sum(super{reg}.celloutput(c).psthsum);
                if numpsthspikes < min_numspikes
                    disp(sprintf('skipped %d %d %d, only %d spikes in psth',super{reg}.celloutput(c).daytetcell,numpsthspikes))
                    continue
                end
                
                % no modulation plot settings
                textcolor = [0 0 0];
                bgcolor = [1 1 1];
                % positive modulation plot settings
                if (super{reg}.celloutput(c).modulation == 1) && (super{reg}.celloutput(c).pvalue <= sig_level)
                    textcolor = [1 .1 0];
                    bgcolor = [250 210 210] / 255;
                end
                % negative modulation plot settings
                if (super{reg}.celloutput(c).modulation == -1) && (super{reg}.celloutput(c).pvalue <= sig_level)
                    textcolor = [.05 .5 1];
                    bgcolor = [176 224 230]/255;
                end

                % plot psth (smoothed or not)
                if ~smooth_flag
                    h = bar(super{reg}.celloutput(c).time,super{reg}.celloutput(c).psthsum);
                else
                    smoothing_width = round(smoothing_length*.001/binsize);   % smoothing width in number of bins
                        kernel = gaussian(smoothing_width,smoothing_width*8);
                    h = bar(super{reg}.celloutput(c).time,smoothvect(super{reg}.celloutput(c).psthsum,kernel));
                end
                axis tight
                set(h(1),'facecolor',clr)
                set(h(1),'edgecolor',clr)
                % print title
                    titlestring=sprintf('%d %d %d',super{reg}.celloutput(c).daytetcell);
                    title(titlestring,'FontSize',12,'FontWeight','bold','Color',textcolor)
                    set(gca,'Color',bgcolor)

                counter = counter+1;
                
            end
            
            %title of previous figure
                    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
                    titlestring=sprintf('%s %s ripple triggered spiking',animals{1},region);
                    title(titlestring)
                    text(0.5, .99,titlestring,'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',18,'FontWeight','bold')
            
        end
    end
    


% Plot regional spike aggregates.

       
if regional_aggregate


hold on

for reg=regions
    
    %grab a time vector
    c = 0;
    while isempty(time)
        c=c+1;
        time = super{reg}.celloutput(c).time;
    end
    
    % collect all spikes
    psthsum{reg} = zeros(size(time));
    for c=1:length(super{reg}.celloutput)
        if ~isempty(super{reg}.celloutput(c).psthsum)
            if (negative_mod == 1 && super{reg}.celloutput(c).modulation >= 0)
                continue
            end
            psthsum{reg} = psthsum{reg} + super{reg}.celloutput(c).psthsum;    
        end
    end
    
    %plot (smoothed or not)
    if smooth_flag
        smoothing_width = round(smoothing_length*.001/binsize);   % smoothing width in number of bins
            kernel = gaussian(smoothing_width,smoothing_width*8);
        h = bar(time,smoothvect(psthsum{reg},kernel));
        bar(time,smoothvect(psthsum{reg},kernel))
        title(sprintf('%s regional aggregate PSTH, smoothed in %d ms',animals{1},smoothing_length),'fontweight','bold','fontsize',16)
    else
        h = bar(time,psthsum{reg});
        bar(time,psthsum{reg})
    end
    
    % color
    if reg==1
        set(h(1),'facecolor',[.5 .5 1],'edgecolor',[.5 .5 1])
    elseif reg==3
        set(h(1),'facecolor',[1 0 0],'edgecolor',[1 0 0])
    else
        set(h(1),'facecolor',[0 1 0],'edgecolor',[0 1 0])
    end
    axis tight
    patch=findobj(h,'Type','patch');
    set(patch,'facealpha',.4)
    set(patch,'edgealpha',.4)

    
end

    
end


% Plot SMOOTHED regional aggregates (check if normalized or not)
    % Measure latency difference.
    % Then check for significance test using Permutation test. 

if smoothed_aggregate_latency
    
figure
hold on

    peakheight = 0;

    % plot PSTH
    for reg=regions
        
        time = super{reg}.celloutput(end).time;
        c1vsc2_reg{reg} = [];
                
        % consolidate spikes over all cells
        for c=1:length(super{reg}.celloutput)
            numpsthspikes = sum(super{reg}.celloutput(c).c1vsc2total);
            % check if at least minimum before including
            if numpsthspikes < min_numspikes
                disp(sprintf('skipped %d %d %d, only %d spikes in psth',super{reg}.celloutput(1).daytetcell,numpsthspikes))
                continue
            else
                c1vsc2_reg{reg} = [c1vsc2_reg{reg} ; super{reg}.celloutput(c).c1vsc2total];
            end
        end
        
        % smooth
        kernel = gaussian(smoothing_width,smoothing_width*8);
        c1vsc2_reg_smooth{reg} = smoothvect(sum(c1vsc2_reg{reg},1),kernel);
        
        % plot PSTH
        h = bar(time,c1vsc2_reg_smooth{reg},1);
        title({sprintf('%s smoothed aggregate PSTHs',animals{1}), ...
               sprintf('%d SD ripples',minthresh_rip)},'fontweight','bold','fontsize',16)
        
        xlabel('time (s)','fontsize',14,'fontweight','bold')
        ylabel('# spikes','fontsize',14,'fontweight','bold')
        
        if reg==1
            set(h(1),'facecolor',[0 0 1],'edgecolor',[0 0 1])
        elseif reg==3
            set(h(1),'facecolor',[1 0 0],'edgecolor',[1 0 0])
        else
            set(h(1),'facecolor',[0 1 0],'edgecolor',[0 1 0])
        end
        
        if reg == 1 | reg == 3
            patch=findobj(h,'Type','patch');
            set(patch,'facealpha',.2)
            set(patch,'edgealpha',.2)
        end
        
        peaktime(reg) = time(find(max(c1vsc2_reg_smooth{reg})==c1vsc2_reg_smooth{reg}));
        peakheight = max([peakheight max(c1vsc2_reg_smooth{reg})]);
        
    end

        
    % draw line at peak time and label with time   
    axis(axis)
    for reg=regions

          line([peaktime(reg) peaktime(reg)],[0 9999],'linewidth',1,'color',[.8 .8 .8]);
            label(1) = {sprintf('%d',ceil(1000*peaktime(reg)))};    % ms
            label(2) = {sprintf('%d',size(c1vsc2_reg{reg},1))}; % # cells
               % use peakheight of region 1 to stick the text
        text(peaktime(reg),-5,label,'fontsize',8,'fontweight','bold')
    end
    
    
    
    % plot amplitude-scaled & smoothed PSTH (use -250 to -200 baseline like
    % Csicsvari--Buzsaki-2000)

    figure
    hold on
    
    for reg=regions

        % scale
        timevector = super{reg}.celloutput(end).time;
        startbaseline = lookup(-.3,timevector);
        endbaseline = lookup(-.25,timevector);
        baseline_amplitude = mean(c1vsc2_reg_smooth{reg}(startbaseline:endbaseline));
        peak_amplitude = max(c1vsc2_reg_smooth{reg});

        c1vsc2_reg_smooth_scaled{reg} = c1vsc2_reg_smooth{reg} / peak_amplitude ;
        
        % plot
        h = bar(time,c1vsc2_reg_smooth_scaled{reg},1);
        title({sprintf('%s smoothed and scaled aggregate PSTHs',animals{1}), ...
               sprintf('%d SD ripples',minthresh_rip)},'fontweight','bold','fontsize',16)
        
        xlabel('time (s)','fontsize',14,'fontweight','bold')
        ylabel('peak normalized # spikes','fontsize',14,'fontweight','bold')
        
        if reg==1
            set(h(1),'facecolor',[0 0 1],'edgecolor',[0 0 1])
            
        elseif reg==3
            set(h(1),'facecolor',[1 0 0],'edgecolor',[1 0 0])
            
        else
            set(h(1),'facecolor',[0 1 0],'edgecolor',[0 1 0])
            
        end
        
            patch=findobj(h,'Type','patch');
            set(patch,'facealpha',.1)
            set(patch,'edgealpha',.1)
            
       

    end
    

    
    
 
    % calculate p-values of latencies
    
    if calculate_p_flag
    
    for p=1:size(pairs,1)
        
        % measure latency difference in ms
        latencies(p,1) = binsize * 1000 * (find(max(c1vsc2_reg_smooth{pairs(p,1)})==c1vsc2_reg_smooth{pairs(p,1)}) - ...
                                           find(max(c1vsc2_reg_smooth{pairs(p,2)})==c1vsc2_reg_smooth{pairs(p,2)}) );
        
        % run permutations between regions -- then calculate latency differences
        % for smoothed histograms to get distribution given null hypothesis
        
        % version 1: permute units' psth
        
        pool = [ c1vsc2_reg{pairs(p,1)} ; c1vsc2_reg{pairs(p,2)} ];
        numregion1 = size(c1vsc2_reg{pairs(p,1)},1);                    % number of xcorr that belong to region 1
        permlatencies{p} = nan(1,numperms);
        
        % permute
        for ii = 1:numperms
            permutation = pool(randperm(size(pool,1)),:);
        
            region1perm = permutation(1:numregion1,:);
            region2perm = permutation((numregion1+1):end,:);
        
            region1smooth = smoothvect(sum(region1perm,1),kernel);
            region2smooth = smoothvect(sum(region2perm,1),kernel);
        
            permlatencies{p}(ii) = find(max(region1smooth)==region1smooth) - ...
                                   find(max(region2smooth)==region2smooth);
                               
            
            
        end
        
        figure
        hold on
        hist(permlatencies{p},50); axis(axis)
            xlabel('# bins','fontsize',14,'fontweight','bold')
            ylabel('# permutations','fontsize',14,'fontweight','bold')

        title(sprintf('perm distribution of latencies, %d %d',pairs(p,:)),'fontweight','bold','fontsize',16)
        % draw line at measured latency
        line([latencies(p,1) latencies(p,1)],[0 9999],'linewidth',3,'color','r');
        
        pvalue = 1 - lookup(latencies(p,1),sort(permlatencies{p}))/length(permlatencies{p});
        latencies(p,2) = pvalue;
        
        % stick pvalue and latency on the graph
            label(1) = {sprintf('latency: %d ms',latencies(p,1))};
            label(2) = {sprintf('p = %0.5g',pvalue)};
            label(3) = {sprintf('(numperms %d)',numperms)};
            N = hist(permlatencies{p},50); axis(axis);
        text(latencies(p,1)-2,0.75*max(N),label,'fontsize',14,'fontweight','bold')
        
    end
    
    end
    
latencies
    
    
end















