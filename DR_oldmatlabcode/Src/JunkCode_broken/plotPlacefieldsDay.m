% plots all cells for all epochs one day at a time

Veqn = '>=2'
minV =  str2num(Veqn(end))
maxstage = 3% [1 2 3]
minVPF = 2 %cm/sec
minPeakPF = 3
lessthan=0
includestates = 6

%Animal selection
%-----------------------------------------------------
animals = {'I1'};
%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filter

days='[11:15]';%,'1:10';



%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

epochfilter{1} = ['isequal($epochtype, ''Run'')'];

velcutoff=2;


cellfilter = '(isequal($area, ''CA1'') && ($meanrate >0))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%timefilter = { {'JY_getlinvelocity', strcat('$velocity <',num2str(velcutoff))} };
timefilter = { {'JY_getriptimes','($nripples == 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},{'JY_getlinvelocity', strcat('$velocity < ',num2str(minVPF))}};
%timefilter = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}, {'JY_getlinvelocity', '$velocity <0.01'}};
%timefilter2 = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter2 = { {'getriptimes','($nripples > 0)',[],'minthresh',2,'cellfilter',cellfilter,'tetfilter',[1 2 3 4]}};

%timefilter = { {'JY_getlinvelocity', '(($velocity) >= 0))', 6} };
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);




animaldir='I1_';
datadir = '/data14/jai/';

daylistmain=unique(f.epochs{1,1}(:,1));
for d=1:length(daylistmain);
    % find of epochs per day
    [drow, dcol]=find(f.epochs{1,1}(:,1)==daylistmain(d)); % get all epochs for the day
    celltettable=[];
    for dataind=1:length(drow);
        dataindsize(dataind)=size(f.data{1,1}{1,drow(dataind)},1); % for each day, concatenate the no. of cells/tetrodes
        celltettable=[celltettable;f.data{1,1}{1,drow(dataind)}];
    end
    
    % get unique celltet pairs
    uniquecelltet=unique(celltettable,'rows');
    uniquetet=unique(uniquecelltet(:,1));
    
    maxcellind=size(uniquecelltet(:,1),1); % maximum no. of cells per day used as max. column number
    for uniqueind=1:length(uniquetet);
        [urow ucol]=find(uniquecelltet(:,1)==uniquetet(uniqueind));
        % first column indicates tetrode, second column specifies which
        % column on the subplot the data should be plotted
        uniqueref(uniqueind,1)=uniquetet(uniqueind);
        uniqueref(uniqueind,2)=length(urow);
    end
    
    figure;
    
    set(gcf,'position',[0 0 2000 500]);
    set(gcf,'PaperPositionMode','auto');
    
    dist=0;
    nrows=length(dcol);
    columns=maxcellind;
    %gap_h=0.0000000001;
    %gap_w=0.00000001;
    gap_h=-0.3;
    gap_w=0.000005;
    
    marg_h=[0.01 0.1];
    marg_w=[0.01 0.1];
    ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
    ax=1;
    
    
    for epochind=1:length(drow);
        % plot circle round rewarded wells
        
        
        % get all data for each day
        % find unique days
        
        i=drow(epochind);
        
        % epoch defined in f.epoch
        day=f.epochs{1,1}(i,1);
        epoch=f.epochs{1,1}(i,2);
        
        
        % get epoch index for the current day
        daylist=f.epochs{1,1}(f.epochs{1,1}(:,1)==day,:);
        [row col]=find(daylist(:,2)==epoch);
        
        tetlast=f.data{1,1}{1,i}(1,1);
        
        
        % cells defined in f.data
        for j=1:size(f.data{1,1}{1,i},1);
            tet=f.data{1,1}{1,i}(j,1);
            cell=f.data{1,1}{1,i}(j,2);
            
            % excluded time defined in f.excluded times
            excludetime=f.excludetime{1,1}{1,i};
            
            [imagedatas,output,pospossize,Wellpos]=PlotopenfieldrateRerouteNoSave(animaldir,tet,cell,day,epoch, 2,1,0,0,1,0,excludetime);
            
            % get current axis depending on tet
            % find which subplot the current data should be on according to uniqueref tet
            plotcol=uniqueref(uniqueref(:,1)==tet,2);
            if tet~=tetlast;
                [uniquerowref dump]=find(uniqueref(:,1)==tetlast); % get the row reference of previous tetrode
                % get actual column ref
                currcol=ax-(row-1)*columns;
                % test if current column is less than what it should be
                % done by comparing the column the current figure (currcol) should be
                % on with the uniqueref table
                shouldbecol=cumsum(uniqueref(:,2))-uniqueref(:,2)+1; % list of columns each tetrode should be plotted on
                if currcol<  shouldbecol(uniquerowref+1);
                    ax=ax+1;
                end
                %                 if j<cumsum(uniqueref(1:uniquerowref,2));
                %                     cumsum(uniqueref(1:uniquerowref+1,2));
                %                 end
            end
            
            axes(ha(ax));
            
            currax=ax;
            
            bmax = max(output.smoothedspikerate(:));
            bmin = min(output.smoothedspikerate(:));
            tett = num2str(tet);
            
            % Firing rate
            image(imagedatas,'CDataMapping', 'scaled');
            hold on;
            for wind=1:size(Wellpos,1);
                plot(Wellpos(wind,2),Wellpos(wind,1),'+w','MarkerSize',3);
            end
            
            
            nc = 1024;
            cmap = jet(nc);
            cmap(1,:) = 1;
            colormap(cmap);
            
            % Generates a title from information provided
            %title(sprintf('Smoothed occupancy normalised spike rate for %s: Day %s Epoch %s Tetrode %s cell %s \n %s spikes',...
            %    animals{1,1}, day, epoch, tet,cell,pospossize),'FontSize',5);
            
            %caxis([-1 bmax]);
            
            % Keeps the proportions of the data
            axis image;
            % no ticks
            set(gca,'xtick',[],'ytick',[]);
            %set(gca,'PlotBoxAspectRatio',[1 1 1]);
            % Labels the sides as 'cm'
            %xlabel('cm');
            %ylabel('cm');
            
            xlimval=[0 size(imagedatas,2)];
            ylimval=[0 size(imagedatas,1)];
            
            distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
            %titlelabel=[0.8*(xlimval(2)-xlimval(1))+xlimval(1) 0.3*(ylimval(2)-ylimval(1))+ylimval(1)];
            
            maxfiringrate=max(imagedatas(:));
            
            text(distlabel(1),distlabel(2),sprintf('ep %s t %s c %s, %s spks %s Hz',num2str(epoch),num2str(tet),num2str(cell),num2str(pospossize),num2str(maxfiringrate,2)),'FontSize',3,'Color','k');
            %text(titlelabel(1),titlelabel(2),sprintf('day %s',num2str(daylist(1,d))),'FontSize',8,'Color','w');
            
            %display(sprintf('tet %s cell %s ax= %s',num2str(tet),num2str(cell),num2str(currax)));
            
            % The tick marks are  apart
            %set(gca, 'XTick', [0:5:maxx-minx]);
            %set(gca, 'YTick', [0:5:maxy-miny]);
            
            % Adds color legend, tick names
            % Make the colorbar, specifying the location and the tick names
            
            %scaleint = [10];
            
            %datascale=[0:scaleint:bmax];
            
            %tickname = [0:scaleint:bmax];
            
            %c = colorbar('location','EastOutside','YTickLabel', tickname);
            
            %ylabel(c,strcat('Hz'));
            
            % Changes the tick intervals of color legend at intervals defined by scalteint
            
            %set(c, 'YTick',datascale);
            
            tetlast=tet;
            
            if j== size(f.data{1,1}{1,i},1) && ax==(row-1)*columns+columns;
                ax=ax+1;
                
            elseif j== size(f.data{1,1}{1,i},1) && ax<=(row-1)*columns+columns;
                ax=(row-1)*columns+columns+1;
            elseif j< size(f.data{1,1}{1,i},1);
                ax=ax+1;
                
            end
            
            
        end
    end
    % Print title for whole figure
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
        1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 0.9,sprintf('Occupancy normalised spike rate of %s for day %s \n excluding periods with < %s cm/s',...
        animals{1,1}, num2str(daylistmain(d)), num2str(velcutoff)),'HorizontalAlignment','center','VerticalAlignment', 'top');
    %                 step=scf/5;
    %                 colorbarlabel=[0:step:scf];
    %                 hcol=colorbar('YTickLabel', num2str(colorbarlabel(2:end)'),'Location','EastOutside');
    %                 %hcol=colorbar('Location','EastOutside');
    %                 cpos=get(hcol,'Position');
    %                 cpos(4)=cpos(4)/2;      % Halve the thickness
    %                 cpos(1)=cpos(1)-0.05;
    %                 cpos(2)=cpos(2)+0.3;   %  Move it down outside the plot
    %                 set(hcol,...
    %                     'Position',cpos)
    
    % Plot your figure
    
    
    
    set(gcf,'PaperPositionMode','auto');
    
    papersize = get(gcf, 'PaperSize');
    
    set(gcf,'PaperSize',[25 15]);
    
    
    % ----Saving----
    % Saves figure as pdf
    % First checks if a folder called Plot exists in the processed data folder,
    % if not, creates it.
    
    cd(f.animal{1,2});
    plotdir = dir('Plot');
    if (isempty(plotdir))
        %an a plot folder needs to be created
        !mkdir Plot
    end
    
    % change to that directory and saves the figure with file name
    % animal_day_epoch
    cd(strcat(f.animal{1,2},'Plot/'));
    figurename = strcat(animals{1,1},'_sspr_d',num2str(daylistmain(d)));
    
    saveas(gcf, figurename, 'pdf');
    
    % Closes the figure
    close;
    
end
clear all;