
% this function can be used for consensus ripples



runscript = 0;


% runscript parameters        
window = [.5 .5]          % size of psth window (in sec)
binsize = .005         % size of bins (in sec)
minthresh_rip = 3 ;     % how big your ripples are
exclusion_dur = 0.5;   % seconds within which consecutive ripples are eliminated
exclusion_nrip = 2;    % this should be the same as your nripples you specify when you call kk_getriptimes 

% analyses
unnorm_unit_psth = 1;
    significance_test = 1;
smoothed_unit_psth = 0;
    normalized_xcorr = 0;       % toggle if you want
regional_aggregate = 0;
    normalized_xcorr = 0;       % toggle if you want
smoothed_aggregate_latency = 0;     % plots smoothed aggregates + calculates latencies
    calculate_p_flag = 0;           % if calculate, set binsize very small!
        numperms = 1000; 

% analysis parameters
smoothing_width = 0;   % std of smoothing gaussian (in # bins)
regions = [1 3 2];
pairs = [1 3 ; 1 2 ; 3 2];
    latencies = nan(size(pairs));
min_numspikes = 25;    % don't process unit if too few total spikes in window around ripples


    

if runscript
    

% Animal Selection
animals = {'Egypt'};

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
iterator = 'kk_multicellanal';

% Filter Creation
ca1f = createfilter('animal', animals, 'epochs', epochfilter, 'cells', ca1cellfilter, 'excludetime', timefilter, 'iterator', iterator);
ca2f = createfilter('animal', animals, 'epochs', epochfilter, 'cells', ca2cellfilter, 'excludetime', timefilter, 'iterator', iterator);
ca3f = createfilter('animal', animals, 'epochs', epochfilter, 'cells', ca3cellfilter, 'excludetime', timefilter, 'iterator', iterator);

% Set Analysis Function
ca1f = setfilterfunction(ca1f, 'dfakk_getriptriggeredspiking', {'ripple','ripples','spikes','tetinfo'},'tetfilter','(isequal($area, ''CA1''))','window',window,'binsize',binsize);
ca2f = setfilterfunction(ca2f, 'dfakk_getriptriggeredspiking', {'ripple','ripples','spikes','tetinfo'},'tetfilter','(isequal($area, ''CA1''))','window',window,'binsize',binsize);
ca3f = setfilterfunction(ca3f, 'dfakk_getriptriggeredspiking', {'ripple','ripples','spikes','tetinfo'},'tetfilter','(isequal($area, ''CA1''))','window',window,'binsize',binsize);


% Run Analysis
ca1f = runfilter(ca1f);
ca2f = runfilter(ca2f);
ca3f = runfilter(ca3f);

end





% Plot.  

% CA1-CA1

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

% Consolidate cells across epochs, within a day

%collect all cell indices, ignoring specific epochs, into daytetcell
indices = [];
for i=1:length(f.output{1})  % iterate over epochs
    indices = [indices ; f.output{1}(i).cellindices];
end
daytetcell = unique(indices(:,[1 3 4]),'rows');

f.celloutput = struct;

for ind=1:size(daytetcell,1)
    f.celloutput(ind).daytetcell = daytetcell(ind,:);
    f.celloutput(ind).c1vsc2 = [];
    f.celloutput(ind).nospikes = 0;
    f.celloutput(ind).noripples = 0;
    f.celloutput(ind).noepochs = 0;
    %consolidate data across epochs
    for ep=1:length(f.output{1})
        if ~isempty(f.output{1}(ep).cellindices)   % epoch may not have any data (fully excluded)
            jj=rowfind(daytetcell(ind,:),f.output{1}(ep).cellindices(:,[1 3 4]));
            if jj
                f.celloutput(ind).c1vsc2 = [f.celloutput(ind).c1vsc2 ; f.output{1}(ep).c1vsc2(jj,:)];
                f.celloutput(ind).nospikes = f.celloutput(ind).nospikes + f.output{1}(ep).nospikes(jj);
                f.celloutput(ind).noripples = f.celloutput(ind).noripples + f.output{1}(ep).noripples(jj);
                f.celloutput(ind).noepochs = f.celloutput(ind).noepochs + 1;
                f.celloutput(ind).time = f.output{1}(ep).time;   % (all time vectors are the same..)
            end
        end
    end
    f.celloutput(ind).c1vsc2total = sum(f.celloutput(ind).c1vsc2,1);
    f.celloutput(ind).c1vsc2totalnorm = f.celloutput(ind).c1vsc2total/sqrt(f.celloutput(ind).nospikes*f.celloutput(ind).noripples);
end
    
    
    super{reg} = f;
   
    
end


% Plot unnormalized raw unit PSTH.
    
    if unnorm_unit_psth
        
        for reg=regions
            if reg==1
                region='CA1';
                clr=[0 0 0];
            elseif reg==2
                region='CA2';
                clr=[0 1 0];
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
                numpsthspikes = sum(super{reg}.celloutput(c).c1vsc2total);
                if numpsthspikes < min_numspikes
                    disp(sprintf('skipped %d %d %d, only %d spikes in psth',super{reg}.celloutput(1).daytetcell,numpsthspikes))
                    continue
                end
                
                % determine whether significantly modulated
                if significance_test
                    keyboard
                        % first obtain distribution of firing from -450
                       
                end
                
                
                h = bar(super{reg}.celloutput(c).time,sum(super{reg}.celloutput(c).c1vsc2total(1,:),1));
                set(h(1),'facecolor',clr)
                set(h(1),'edgecolor',clr)
                titlestring=sprintf('%d %d %d',super{reg}.celloutput(c).daytetcell);
                title(titlestring)
                
                counter = counter+1;
            end
            
            %title of previous figure
                    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
                    titlestring=sprintf('%s %s ripple triggered spiking',animals{1},region);
                    title(titlestring)
                    text(0.5, .99,titlestring,'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',18,'FontWeight','bold')
            
        end
    end
    
% Plot smoothed individual unit PSTH.
    % Then measure latencies.
    % First smooth individual hist, then plot.
    % Then then find latencies near 0, then collect latencies and plot.
        % this last part is UNWRITTEN

    if smoothed_unit_psth
        
        for reg=regions
            if reg==1
                region='CA1';
                clr=[0 0 0];
            elseif reg==2
                region='CA2';
                clr=[0 1 0];
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
                    titlestring=sprintf('%s %s %d ripple PSTH // binsize %d, smooth %d ms',animals{1},region,minthresh_rip, 1000*binsize,1000*smoothing_width*binsize);

                    title(titlestring)
                    text(0.5, .99,titlestring,'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',18,'FontWeight','bold')
                    
                    figure
                    counter = 1;
                    titlestring=sprintf('%s %s %d ripple PSTH // binsize %d, smooth %d ms',animals{1},region,minthresh_rip, 1000*binsize,1000*smoothing_width*binsize);

                    title(titlestring)
                end
                subplot(6,10,counter)
       
                % SKIP plotting if not >= 50 spikes
                numpsthspikes = sum(super{reg}.celloutput(c).c1vsc2total);
                if numpsthspikes < 0
                    disp(sprintf('skipped %d %d %d, only %d spikes in psth',super{reg}.celloutput(1).daytetcell,numpsthspikes))
                    continue
                end
                
                %smooth and plot
                kernel = gaussian(smoothing_width,smoothing_width*8);
                if normalized_xcorr == 1
                    h = bar(super{reg}.celloutput(c).time,smoothvect(sum(super{reg}.celloutput(c).c1vsc2totalnorm(1,:),1),kernel));
                    ylim([0 0.001])
                else
                    h = bar(super{reg}.celloutput(c).time,smoothvect(sum(super{reg}.celloutput(c).c1vsc2total(1,:),1),kernel));
                end

                
                set(h(1),'facecolor',clr)
                set(h(1),'edgecolor',clr)
                titlestring=sprintf('%d %d %d',super{reg}.celloutput(c).daytetcell);
                title(titlestring)
                
                counter = counter+1;
            end
            
                    %title of previous figure
                    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
                    titlestring=sprintf('%s %s %d ripple PSTH // binsize %d, smooth %d ms',animals{1},region,minthresh_rip, 1000*binsize,1000*smoothing_width*binsize);
                    title(titlestring)
                    text(0.5, .99,titlestring,'HorizontalAlignment','center','VerticalAlignment', 'top','FontSize',18,'FontWeight','bold')
            
            
        end
    end


% Plot regional spike aggregates.

       
if regional_aggregate

figure
hold on

    
for reg=regions
    
    time = super{reg}.celloutput(end).time;
    c1vsc2_reg{reg} = [];
    
    for c=1:length(super{reg}.celloutput)
        if normalized_xcorr
            c1vsc2_reg{reg} = [c1vsc2_reg{reg} ; super{reg}.celloutput(c).c1vsc2totalnorm];
        else
            c1vsc2_reg{reg} = [c1vsc2_reg{reg} ; super{reg}.celloutput(c).c1vsc2total];    
        end
    end
    
    h = bar(time,sum(c1vsc2_reg{reg},1));
    
    title(sprintf('%s smoothed regional aggregate PSTHs',animals{1},'fontweight','bold','fontsize',16))
    
    if reg==1
        set(h(1),'facecolor',[0 0 1],'edgecolor',[0 0 1])
    elseif reg==3
        set(h(1),'facecolor',[1 0 0],'edgecolor',[1 0 0])
    else
        set(h(1),'facecolor',[0 1 0],'edgecolor',[0 1 0])   
    end
    
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















