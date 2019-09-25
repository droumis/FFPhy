global subplot_count;
subplot_count = 1;

Veqn = '>=0';
minV =  str2num(Veqn(end));
maxstage = 3; % [1 2 3]
minVPF = 3; %cm/sec
minPeakPF = 3;
lessthan=0;
includestates = 6;

%Animal selection
%-----------------------------------------------------
animals = {'I1'};
%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filterionno

days='[6:15]';%,'1:10';
%days = '[1:1]';
%days = '[9:9]';

%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

epochfilter{1} = ['isequal($epochtype, ''Run'')'];




cellfilter = '(isequal($area, ''CA1'') && ($meanrate >0))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%timefilter = { {'JY_getlinvelocity', '($velocity <0.07) & ($velocity >-0.07)'} };

timefilter = { {'JY_getriptimes','($nripples > 1)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
    {'JY_getlinvelocity', strcat('$velocity < ',num2str(minVPF))}};
%timefilter = { {'JY_getriptimes','($nripples > 1)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'}};

%timefilter = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}, {'JY_getlinvelocity', '$velocity <0.01'}};
%timefilter2 = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter2 = { {'getriptimes','($nripples > 0)',[],'minthresh',2,'cellfilter',cellfilter,'tetfilter',[1 2 3 4]}};

%timefilter = { {'JY_getlinvelocity', '(($velocity) >= 0))', 6} };
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);

%only include cells with placefields
%if minPeakPF>0
%    includecells = calcincludecells(minVPF, minPeakPF, animals, epochfPF, lessthan);
%5    f = excludecellsf(f, includecells);
%end
%-----------------------------------------------------------


%run function- single cells
%--------------------------------------------
iterator = 'JY_singleepochanal';

f = setfilteriterator(f,iterator);

f = setfilterfunction(f, 'JY_getposition', {'data','linpos'});
% out = plottrajdata(index, excludetimes, spikes, linpos, includestates, minV, varargin)

f = runfilter(f);

daylist=unique(f.epochs{1,1}(:,1));





for i=1:size(daylist,1)
    
    % number of epochs for i
    epochlist=f.epochs{1,1}(f.epochs{1,1}(:,1)==daylist(i),2);
    
    % set dimension of subplots
    figure;
    
    set(gcf,'position',[0 0 800 1200]);
set(gcf,'PaperPositionMode','auto');
    
    
    
    
    nrows=size(epochlist,1);
    columns=3;
    %gap_h=0.0000000001;
    %gap_w=0.00000001;
    gap_h=0.05;
    gap_w=0.05;
    
    marg_h=[0.05 0.10];
    marg_w=[0.01 0.01];
    ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
    ax=1;
    
    
    for j=1:size(epochlist,1);
    
        dataindex=find(f.epochs{1,1}(:,1)==daylist(i) & f.epochs{1,1}(:,2)==epochlist(j));
        
    % plot where ripples are       
    axes(ha(ax));
        
    plot(f.output{1,1}{1,dataindex}.allpositions(:,2),f.output{1,1}{1,dataindex}.allpositions(:,3),'.','Color',[0.8 0.8 0.8]);
    hold on;
    plot(f.output{1,1}{1,dataindex}.rippleposition(:,1),f.output{1,1}{1,dataindex}.rippleposition(:,2),'.r');
    set(gca,'xtick',[],'ytick',[]);
    set(gca,'xlim',[10 190],'ylim',[-40 140]);
    xlimval=get(gca,'xlim');
    ylimval=get(gca,'ylim');
    
    axis square;
    seglabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.10*(ylimval(2)-ylimval(1))+ylimval(1)];
                    % counts number of segments using size of existref, or
                    % times rat leaves the intersection zones
                    text(seglabel(1),seglabel(2),sprintf('epoch %s \n n=%s',num2str(epochlist(j)),num2str(size(f.output{1,1}{1,dataindex}.rippleposition(:,1),1))),'FontSize',10);
    
    ax=ax+1;
    
    % plot as cumulative histogram
    axes(ha(ax));
    edges=floor(min(f.output{1,1}{1,dataindex}.ripplevelocity)):1:ceil(max(f.output{1,1}{1,dataindex}.ripplevelocity));
    h=histc(f.output{1,1}{1,dataindex}.ripplevelocity,edges);
    velhist=cumsum(h);
    bar(edges,velhist./max(velhist));
    xlabel('Linear Velocity (cm/s)');
    ylabel('Proportion');
    set(gca,'FontSize',8);
    axis square;
    
    
    % calculate how many ripples occur in one place
    ax=ax+1;
    axes(ha(ax));
    
     minx = -5;
    maxx = 380*0.5;
    binx = (minx:1:maxx);
    miny = -5;
    maxy = 280*0.5;
    biny = (miny:1:maxy);
    [output.occupancy output.xticks output.yticks] = hist2(f.output{1,1}{1,dataindex}.rippleposition(:,1),f.output{1,1}{1,dataindex}.rippleposition(:,2) , binx, biny);
    std=1;
    g = gaussian2(std,(6*std));
        output.occupancy = filter2(g,(output.occupancy));
    nc = 32;
            cmap = jet(nc);
            cmap(1,:) = 1;
            colormap(cmap);
    imagesc(flipud(output.occupancy));
    colorbar;
    axis off;
    axis image;
    set(gca,'xtick',[],'ytick',[]);
    
    ax=ax+1;
    end
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
        1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 0.95,sprintf('Ripples of %s for day %s ',...
        animals{1,1}, num2str(daylist(i))),'HorizontalAlignment','center','VerticalAlignment', 'top');
     set(gcf,'PaperPositionMode','auto');
    
    
    
    set(gcf,'PaperSize',[10 15]);
    cd(strcat(f.animal{1,2},'Plot/'));
    figurename = strcat(animals{1,1},'_ripples',num2str(daylist(i)));
    saveas(gcf, figurename, 'pdf');
     close;
end
