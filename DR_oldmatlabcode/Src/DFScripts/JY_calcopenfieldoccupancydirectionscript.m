global subplot_count;
subplot_count = 1;

Veqn = '>=3';
minV =  str2num(Veqn(end));
maxstage = 3; % [1 2 3]
minVPF = 0; %cm/sec
minPeakPF = 3;
lessthan=3;
includestates = 6;

%Animal selection
%-----------------------------------------------------
animals = {'I1'};
%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filterionno

days='[8]';%,'1:10';
%days = '[1:1]';
%days = '[9:9]';

%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

%epochfilter{1} = [''];

epochtype='Run';

%epochfilter{1} = ['isequal($epochtype, ''Run'')'];
epochfilter{1} = ['isequal($epoch,  4)'];

%cellfilter = '(isequal($area, ''CA1'') && ($meanrate >0 ) )'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
cellfilter = '(isequal($area, ''CA1'') )'  ;

%timefilter = { {'JY_getlinvelocity', '($velocity <0.07) & ($velocity >-0.07)'} };

%timefilter = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'}};
timefilter = { {'JY_getriptimes','($nripples == 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
    {'JY_getlinvelocity', strcat('$velocity > ',num2str(minVPF))}};
%timefilter = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}, {'JY_getlinvelocity', '$velocity <0.01'}};
%timefilter2 = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter2 = { {'getriptimes','($nripples > 0)',[],'minthresh',2,'cellfilter',cellfilter,'tetfilter',[1 2 3 4]}};

%timefilter = { 'JY_getlinvelocity', '$velocity > 4' };

f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);
%f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'excludetimefilter', timefilter);

%f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter);
%only include cells with placefields
%if minPeakPF>0
%    includecells = calcincludecells(minVPF, minPeakPF, animals, epochfPF, lessthan);
%5    f = excludecellsf(f, includecells);
%end
%-----------------------------------------------------------


%run function- single cells
%--------------------------------------------
iterator = 'singlecellanal';

f = setfilteriterator(f,iterator);

f=setfilterfunction(f, 'JY_calcopenfieldoccupancydirectional', {'spikes','data','linpos'});

f=runfilter(f);

%outdata=out.output{1,1};



cd(strcat(f.animal{1,2}));

filename = strcat(animals{1,1},'plinfo',num2str(6));

save(filename,'f');


% plotting

animaldir='CML21_';
datadir = '/data14/jai/';





daylistmain=unique(f.epochs{1,1}(:,1));
for d=1:length(daylistmain);
    
    % find of epochs per day
    [drow, dcol]=find(f.epochs{1,1}(:,1)==daylistmain(d)); % get all epochs for the day
    celltettable=[];
    
    for dataind=1:length(drow);
        dataindsize(dataind)=size(f.data{1,1}{1,drow(dataind)},1); % for each day, concatenate the no. of cells/tetrodes
        dayepoch=repmat(f.epochs{1,1}(dataind,:), length(f.data{1,1}{1,drow(dataind)}),1);
        celltettable=[celltettable;[f.data{1,1}{1,drow(dataind)} dayepoch]];

    end
    
    
    
    %     % get unique celltet pairs
    %     uniquecelltet=unique(celltettable,'rows');
    %     uniquetet=unique(uniquecelltet(:,1));
    %
    %     maxcellind=size(uniquecelltet(:,1),1); % maximum no. of cells per day used as max. column number
    %
    %     for uniqueind=1:length(uniquetet);
    %         [urow ucol]=find(uniquecelltet(:,1)==uniquetet(uniqueind));
    %         % first column indicates tetrode, second column specifies which
    %         % column on the subplot the data should be plotted
    %         uniqueref(uniqueind,1)=uniquetet(uniqueind);
    %         uniqueref(uniqueind,2)=length(urow);
    %     end
    
    %     figure;
    %
    %     set(gcf,'position',[0 0 2000 500]);
    %     set(gcf,'PaperPositionMode','auto');
    %
    %     dist=0;
    %     nrows=length(dcol);
    %     columns=maxcellind;
    %     gap_h=0.0000000001;
    %     %gap_w=0.00000001;
    %     %gap_h=-0.2;
    %     gap_w=0.000005;
    %
    %     marg_h=[0.01 0.1];
    %     marg_w=[0.01 0.1];
    %     ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
    %     ax=1;
    %
    
    for epochind=1:length(celltettable);
        
        % get all data for each day
        % find unique days
        
        i=epochind;
        
        % epoch defined in f.epoch
        day=celltettable(i,3);
        epoch=celltettable(i,4);
        
        % See if day number needs 0
        dsz = '';
        if (day < 10)
            dsz = '0';
        end
        dayt = num2str(day);
        
        
        % Loads the position data
        posfilename = strcat(f.animal{1,2},animals{1,1},'data',dsz,dayt,'.mat');
        load(posfilename);
        
        posfilename = strcat(f.animal{1,2},animals{1,1},'linpos',dsz,dayt,'.mat');
        load(posfilename);
        
        % get rewarded wells
%         Wells=data{1,day}{1,epoch}.Wellinfo.rewardedwells;
%         % get reward well positions
%         if strmatch(epochtype,'Run');
%             Wellpos=linpos{1,day}{1,epoch}.wellSegmentInfo.wellCoord;
%         end
        
        % get segment intersections
        
        intersections=unique([linpos{1,day}{1,epoch}.segmentInfo.segmentCoords(:,[1 2]);linpos{1,day}{1,epoch}.segmentInfo.segmentCoords(:,[3 4])],'rows');
        
        % get epoch index for the current day
        %         daylist=f.epochs{1,1}(f.epochs{1,1}(:,1)==day,:);
        %         [row col]=find(daylist(:,2)==epoch);
        %
        %tetlast=f.data{1,1}{1,i}(1,1);
        
        
        %         % cells defined in f.data
        %         for j=1:size(f.data{1,1}{1,i},1);
        tet=celltettable(i,1);
        cell=celltettable(i,2);
        
        % excluded time defined in f.excluded times
        %excludetime=f.excludetime{1,1}{1,i};
        
        %[imagedatas,output,pospossize,Wellpos]=PlotopenfieldrateRerouteNoSave(animaldir,tet,cell,day,epoch, 2,1,0,0,1,0,excludetime);
        
        % get current axis depending on tet
        % find which subplot the current data should be on according to uniqueref tet
        
        figure;
        
        set(gcf,'position',[0 0 1000 1000]);
        set(gcf,'PaperPositionMode','auto');
        
        dist=0;
        nrows=6;
        columns=4;
        gap_h=0.0000000001;
        %gap_w=0.00000001;
        %gap_h=-0.2;
        gap_w=0.000005;
        
        marg_h=[0.01 0.1];
        marg_w=[0.01 0.1];
        ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
        ax=1;
        
        
        
        for jj=1:length(f.output{1,1}{1,i}.data);
            
            
            
            
            
            %             plotcol=uniqueref(uniqueref(:,1)==tet,2);
            %             if tet~=tetlast;
            %                 [uniquerowref dump]=find(uniqueref(:,1)==tetlast); % get the row reference of previous tetrode
            %                 % get actual column ref
            %                 currcol=ax-(row-1)*columns;
            %                 % test if current column is less than what it should be
            %                 % done by comparing the column the current figure (currcol) should be
            %                 % on with the uniqueref table
            %                 shouldbecol=cumsum(uniqueref(:,2))-uniqueref(:,2)+1; % list of columns each tetrode should be plotted on
            %                 if currcol<  shouldbecol(uniquerowref+1);
            %                     ax=ax+1;
            %                 end
            %                 %                 if j<cumsum(uniqueref(1:uniquerowref,2));
            %                 %                     cumsum(uniqueref(1:uniquerowref+1,2));
            %                 %                 end
            %             end
            
            axes(ha(ax));
            
            
            
            
            
            currax=ax;
            
            %             bmax = max(output.smoothedspikerate(:));
            %             bmin = min(output.smoothedspikerate(:));
            tett = num2str(tet);
            
            % start and end intersections
            startintersection=f.output{1,1}{1,i}.directionlist(jj,1);
            endintersection=f.output{1,1}{1,i}.directionlist(jj,2);
            
            % data is in outputdata arranged linearly, need to generate a
            % lookup inded from cell/tet table in f.data
            
            %dataind=(i-1)*size(f.data{1,1}{1,i},1)+j;
            
            binx=f.output{1,1}{1,i}.data{1,jj}.twoDxticks;
            biny=f.output{1,1}{1,i}.data{1,jj}.twoDyticks;
            
          
            
            imagedatas=f.output{1,1}{1,i}.data{1,jj}.twoDsmoothedspikerate;
            %imagedatas=f.output{1,1}{1,i}.data{1,jj}.occupancy;
            
            spikecount=sum(f.output{1,1}{1,i}.data{1,jj}.twoDspikes(:));
            
            if isempty(imagedatas)
                imagedatas=zeros(size(f.output{1,1}{1,i}.data{1,jj}.twoDoccupancy));
                spikecount=0;
            end
            
            %imagedatas(imagedatas>3)=3;
            
            
            
            nc = 1024;
            cmap = jet(nc);
            cmap(1,:) = 1;
            colormap(cmap);
            
            % Generates a title from information provided
            %title(sprintf('Smoothed occupancy normalised spike rate for %s: Day %s Epoch %s Tetrode %s cell %s \n %s spikes',...
            %    animals{1,1}, day, epoch, tet,cell,pospossize),'FontSize',5);
            
            maxfiringrate=max(imagedatas(:));
            
            %caxis([-1 maxfiringrate]);
            % Firing rate
            imagesc(flipud(imagedatas),[-2 10]);
            %imagesc(fliplr(imagedatas));
            %imagesc(imagedatas);
            
            %set(gca,'YDir','normal')
            hold on;
            %             if strmatch(epochtype,'Run');
            %             [output.welloccupancy output.xticks output.yticks]=hist2(Wellpos(:,1),Wellpos(:,2), binx, biny);
            %             output.welloccupancy = flipud(output.welloccupancy);
            %             [wella, wellb]=find(output.welloccupancy==1);
            %             Wellposs=[wella,wellb];
            %             for wind=1:size(Wellpos,1);
            %                 plot(Wellposs(wind,2),Wellposs(wind,1),'+w','MarkerSize',3);
            %             end
            
            % plot intersections
            [output.intersections output.xticks output.yticks]=hist2(intersections(:,1),intersections(:,2), binx, biny);
            output.intersections = flipud(output.intersections);
            [wella, wellb]=find(output.intersections==1);
            intersectionspos=[wella,wellb];
            for wind=1:size(intersectionspos,1);
                plot(intersectionspos(wind,2),intersectionspos(wind,1),'+r','MarkerSize',3);
            end
           
        
        % Keeps the proportions of the data
        axis image;
        % no ticks
        set(gca,'xtick',[],'ytick',[]);
        
        
        xlimval=[0 size(imagedatas,2)];
        ylimval=[0 size(imagedatas,1)];
        
        distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
        
        text(distlabel(1),distlabel(2),sprintf('from %s to %s, %s spks %s Hz',num2str(startintersection),num2str(endintersection),num2str(spikecount),num2str(maxfiringrate,2)),'FontSize',12,'Color','k');
        
        %             tetlast=tet;
        %             if j== size(f.data{1,1}{1,i},1) && ax==(row-1)*columns+columns;
        %                 ax=ax+1;
        %
        %             elseif j== size(f.data{1,1}{1,i},1) && ax<=(row-1)*columns+columns;
        %                 ax=(row-1)*columns+columns+1;
        %             elseif j< size(f.data{1,1}{1,i},1);
        %                 ax=ax+1;
        %
        %             end
        
         ax=ax+1;
        
        
        
        
        
        
        
        
    end
  ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
        1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 0.95,sprintf('Linear place fields of %s for day %s epoch %s tetrode %s cell %s ',...
        animals{1,1}, num2str(daylistmain(d)), num2str(epoch), num2str(celltettable(i,1)),num2str(celltettable(i,2))),'HorizontalAlignment','center','VerticalAlignment', 'top');
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
    
    
    
end
% Print title for whole figure




% ----Saving----
% Saves figure as pdf
% First checks if a folder called Plot exists in the processed data folder,
% if not, creates it.

%     cd(f.animal{1,2});
%     plotdir = dir('Plot');
%     if (isempty(plotdir))
%         %an a plot folder needs to be created
%         !mkdir Plot
%     end
%
%     % change to that directory and saves the figure with file name
%     % animal_day_epoch
%     cd(strcat(f.animal{1,2},'Plot/'));
%     figurename = strcat(animals{1,1},'_sspr_d_',num2str(daylistmain(d)));
%
%     saveas(gcf, figurename, 'pdf');
%
%     %Closes the figure
%     close;

end
%clear all;
