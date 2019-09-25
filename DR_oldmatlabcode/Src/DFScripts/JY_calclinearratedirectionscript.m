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

epochfilter{1} = ['isequal($epochtype, ''Run'')'];
%epochfilter{1} = ['isequal($epoch,  4)'];

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

cd(strcat(f.animal{1,2}));

day=str2num(days);
dsz = '';
if (day < 10)
    dsz = '0';
end

dayt = num2str(day);

filename = strcat(animals{1,1},'plinfo',dsz,dayt);

save(filename,'f');

%
%write out place field summary into summary table
%{day}{epoch}{tetrode}{cell}[list of 12 segments, 0=no field,1=field in direction one, 2=field in direction two, 3=field in both directions]
% field defined as firing rate about 3Hz

placefieldlookup={};

% plotting

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
        
        
        % get segment intersections
        
        intersections=unique([linpos{1,day}{1,epoch}.segmentInfo.segmentCoords(:,[1 2]);linpos{1,day}{1,epoch}.segmentInfo.segmentCoords(:,[3 4])],'rows');
        
        tet=celltettable(i,1);
        cell=celltettable(i,2);
        
        % find the maxfiring rate of the cell over all directions
       maxfiringrate=max(cellfun(@(x) max(x.linearsmoothedrate), f.output{1,1}{1,i}.data));
       
        
        
        % get current axis depending on tet
        % find which subplot the current data should be on according to uniqueref tet
        
        figure;
        
        set(gcf,'position',[0 0 500 500]);
        set(gcf,'PaperPositionMode','auto');
        
        dist=0;
        nrows=6;
        columns=2;
        gap_h=0.05;
        gap_w=0.05;
        
        marg_h=[0.01 0.1];
        marg_w=[0.01 0.1];
        ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
        ax=1;
        
        jj=1;
        
        while jj<length(f.output{1,1}{1,i}.directionlist);
            
            axes(ha(ax));
            
            currax=ax;
            
            tett = num2str(tet);
            
            directionlist=f.output{1,1}{1,i}.directionlist;
            % plot the smoothed spike rate for each direction
            
            directiondata1=f.output{1,1}{1,i}.data{1,jj}.linearsmoothedrate;
            directiondata2=f.output{1,1}{1,i}.data{1,jj+1}.linearsmoothedrate;
            linearspacestep=f.output{1,1}{1,i}.data{1,jj}.linearspacestep;
            xspacing=[0:0.02:1].*linearspacestep;
            
            % find out if placefield exists on current segment (rate>3Hz)
            
            % max firing rate on track
            
            combinedirections=[ directiondata1'; directiondata2'];
            %maxfiringrate=max(combinedirections,[],1);
            
            if maxfiringrate>=5;
                
                if ~isempty(find(directiondata1>0.1*maxfiringrate));
                    direction1placefield=1;
                else
                    direction1placefield=0;
                end
                
                if ~isempty(find(directiondata2>0.1*maxfiringrate));
                    
                    direction2placefield=2;
                else
                    
                    direction2placefield=0;
                end
                
                placefieldcode=direction1placefield+direction2placefield;
                
                placefieldlookup{day}{epoch}{tet}{cell}(1,f.output{1,1}{1,i}.segmentlist(jj))=placefieldcode;
                
            end
            
            %plot(xspacing,directiondata1,'r',xspacing,directiondata2,'b');
            ylim([0 30]);
            %axis equal;
            
            xlimval=get(gca,'Xlim');
            ylimval=get(gca,'Ylim');
            
            distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
            
            startintersection=f.output{1,1}{1,i}.directionlist(jj,1);
            endintersendtion=f.output{1,1}{1,i}.directionlist(jj,2);
            
            text(distlabel(1),distlabel(2),sprintf('Segment %s %s to %s',num2str(ceil(jj/2)),num2str(startintersection),num2str(endintersendtion)),'FontSize',12,'Color','r');
            
            jj=jj+2;
            
            ax=ax+1;
            
        end
        
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
            1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        text(0.5, 0.95,sprintf('Linear place fields of %s for day %s epoch %s tetrode %s cell %s ',...
            animals{1,1}, num2str(daylistmain(d)),num2str(epoch), num2str(celltettable(i,1)),num2str(celltettable(i,2))),'HorizontalAlignment','center','VerticalAlignment', 'top');
        
        set(gcf,'PaperPositionMode','auto');
        
        papersize = get(gcf, 'PaperSize');
        
        set(gcf,'PaperSize',[25 15]);
        
    end
    % Print title for whole figure
    
%     filename = strcat(animals{1,1},'placefield',dsz,num2str(daylistmain(d)),'.mat');
%     
%     placefield=placefieldlookup;
%     
%     save(filename,'placefield');
    
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
    %
    %
    % clear all;
    
end
