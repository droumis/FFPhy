% plots the combined normalised place fields of all places cells active in
% a rippple



Veqn = '>=0'
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

days='[11]';%,'1:10';



%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

%epochfilter{1} = ['isequal($epochtype, ''Run'')'];
epochfilter{1} = ['isequal($epoch, 2)'];
cellfilter = '(isequal($area, ''CA1'') && ($meanrate <7))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%timefilter = { {'JY_getlinvelocity', '$velocity <2'} };
timefilter = { {'JY_getriptimes','($nripples ==0)', [], 3,'cellfilter', '(isequal($area, ''CA1''))'}};
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
iterator = 'singlecellanal';

f = setfilteriterator(f,iterator);

out=setfilterfunction(f, 'JY_calcopenfieldoccupancy', {'spikes','data'});

out=runfilter(out);

% place fields are stored here
outdata=out.output{1,1};


% Plotting

daylistmain=unique(f.epochs{1,1}(:,1));
animaldir='I1_';
datadir = '/data14/jai/';
animalname = animals{1,1};
daylistmain=unique(f.epochs{1,1}(:,[1 2]),'rows');

for d=1:size(daylistmain,1);
    
    day=daylistmain(d,1)
    epoch=daylistmain(d,2);
    
    dsz = '';
    if (day < 10)
        dsz = '0';
    end
    
    dayt = num2str(day);
    
    % load place field file (plf)
    %     filename=strcat(datadir,animaldir,'/',animalname,'plf',dsz,dayt,'.mat');
    %     load(filename);
    
    % Load the mat file, needed for barrier information
    sfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'data',dsz,dayt,'.mat');
    load(sfilename);
    
    posfilename = strcat(datadir,animaldir,'/',animals{1,1},'linpos',dsz,dayt,'.mat');
    load(posfilename);
    
    % get rewarded wells
    Wells=data{1,day}{1,epoch}.Wellinfo.rewardedwells;
    % get reward well positions
    Wellpos=linpos{1,day}{1,epoch}.wellSegmentInfo.wellCoord;
    
    
    % barrier label
    barrierevents=data{1,day}{1,epoch}.Events.Barrier;
    
    
    
    % load spike summary file
    filename=strcat(f.animal{2},f.animal{1},'ripple_cell_count',dsz,num2str(day),'.mat');
    load(filename);
    % count of active cells in each ripple
    ripplecellcount=ripple_cell_count{day}{epoch};
    % find all ripples with more than one cell active
    rippleind=find(ripplecellcount>0);
    
    % setup figure
    figure;
    
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperSize',[2500 1000]);
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'position',[10 10 2300 800]);
    
    
    % find how many plots are needed
    tmpdata={};
    rippcounter=0;
    dataref={};
    
    for plotind=1:length(rippleind);
        
        % get a list of cells active in each ripple
        % spikes are found in ripple_spikeslist
        uniquecelllist=unique(ripple_spikeslist{day}{epoch}{2,rippleind(plotind)}(:,[8 9]),'rows');
        cellsplacefield={};
        cellcounter=0;
        
        % get the corresponding place field of each cell
        
        for cellind=1:size(uniquecelllist,1);
            
            % get place field information
            %             if max(celldata{1,day}{1,epoch}{1,uniquecelllist(cellind,1)}{1,uniquecelllist(cellind,2)}.plf(:))<100;
            
            % find epoch
            % find tet
            currtetind=uniquecelllist(cellind,1);
            % find cell
            currcellind=uniquecelllist(cellind,2);
            % find corresponding epoch index in f.data
            epochind=find(f.epochs{1,1}(:,2)==epoch);
            % find tet/cell reference in f.data
            currdataind=find(f.data{1,1}{1,epochind}(:,1)==currtetind & f.data{1,1}{1,epochind}(:,2)==currcellind);
            if ~isempty(currdataind)
                
                % find corresponding placefield in outdata
                outdataind=(epochind-1)*size(f.data{1,1}{1,epochind},1)+currdataind;
                
                cellsplacefield{1,cellind}=outdata(1,outdataind).smoothedspikerate;
                binx=outdata(1,outdataind).xticks;
                biny=outdata(1,outdataind).yticks;
                %cellsplacefield{1,cellind}(isnan(cellsplacefield{1,cellind}))=0;
                
                cellcounter=cellcounter+1;
            else
                cellind=cellind+1;
            end
            
            %else cellsplacefield{1,cellind}=zeros(size(celldata{1,day}{1,epoch}{1,uniquecelllist(cellind,1)}{1,uniquecelllist(cellind,2)}.plf_norm));
            
            %end
            %cellind=cellind+1;
            
        end
        
        imagedata=max(cat(3,cellsplacefield{:}),[],3);
        imagedatamin=min(cat(3,cellsplacefield{:}),[],3);
        imagedata(imagedatamin==-1)=-1;
        
        if ~sum(imagedata(:))==0;
            tmpdata{rippcounter+1}=imagedata; % store current image in data
            dataref{rippcounter+1}(1,1)=rippleind(plotind);
            dataref{rippcounter+1}(1,2)=cellcounter;
            rippcounter=rippcounter+1;
            plotind=plotind+1;
            
        else
            plotind=plotind+1;
            
        end
    end
    
    % barrier label
    barrierevents=data{1,day}{1,epoch}.Events.Barrier;
    %plot
    
    dist=0;
    columns=50;
    nrows=ceil(rippcounter/columns);
    %gap_h=0.0000000001;
    %gap_w=0.00000001;
    gap_h=-0.00005;
    gap_w=0.000005;
    
    marg_h=[0.01 0.1];
    marg_w=[0.01 0.1];
    ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
    ax=1;
    
    % go through each ripple with more than one cell active and plot place
    % fields of active cells
    plotind=1;
    nc = 1024;
        cmap = jet(nc);
        cmap(1,:) = 1;
        colormap(cmap);
        
    for plotind=1:length(tmpdata);
        
        axes(ha(ax));
        
        
        imagedata=tmpdata{plotind};
        % max color scale value
        amax = max(imagedata(:));
        
        caxis([-1 1]);
        imagesc(flipud(imagedata));
        hold on;
        
        
        % plot where the ripple happend
        TS=min(ripple_spikeslist{day}{epoch}{2,dataref{plotind}(1,1)}(:,1))*10000;
        ripx =mean(ripple_spikeslist{day}{epoch}{2,dataref{plotind}(1,1)}(:,2));
        ripy=mean(ripple_spikeslist{day}{epoch}{2,dataref{plotind}(1,1)}(:,3));
        [output.rippleoccupancy output.xticks output.yticks]=HIST2(ripx,ripy, binx, biny);
        output.rippleoccupancy = flipud(output.rippleoccupancy);
        [ripa, ripb]=find(output.rippleoccupancy==1);
        plot(ripb,ripa,'xm','MarkerSize',1);
        
        [output.welloccupancy output.xticks output.yticks]=HIST2(Wellpos(:,1),Wellpos(:,2), binx, biny);
        output.welloccupancy = flipud(output.welloccupancy);
        [wella, wellb]=find(output.welloccupancy==1);
        Wellposs=[wella,wellb];
        for wind=1:size(Wellpos,1);
            plot(Wellposs(wind,2),Wellposs(wind,1),'+w','MarkerSize',1);
        end
        
        axis image;
        % no ticks
        set(gca,'xtick',[],'ytick',[]);
        
        
        xlimval=[0 size(imagedata,2)];
        ylimval=[0 size(imagedata,1)];
        
        distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
        %titlelabel=[0.8*(xlimval(2)-xlimval(1))+xlimval(1) 0.3*(ylimval(2)-ylimval(1))+ylimval(1)];
        rippnumber=dataref{plotind}(1,1);
        cellcount=dataref{plotind}(1,2);
        clear imagedata;
        %find which trial the ripple happend in
        test=[];
        trialno=0;
        TS=min(ripple_spikeslist{day}{epoch}{2,dataref{plotind}(1,1)}(:,1))*10000;
            test(:,1)=data{1,day}{1,epoch}.Run(:,3)<=TS;
            test(:,2)=data{1,day}{1,epoch}.Run(:,4)>=TS;
            find(test(:,1)==1 & test(:,2)==1);
            trialno=find(test(:,1)==1 & test(:,2)==1);
        
        if ~isempty(barrierevents);
            
            % get run start and run end times
            TS=min(ripple_spikeslist{day}{epoch}{2,dataref{plotind}(1,1)}(:,1))*10000;
            TE=max(ripple_spikeslist{day}{epoch}{2,dataref{plotind}(1,1)}(:,1))*10000;
            BS=barrierevents(:,1);
            BE=barrierevents(:,1)+barrierevents(:,4);
            barriertest=(TS >= BS & TS < BE) | (TE <= BE & TE > BS) | (TS <= BS & TE >= BE);
            if sum(barriertest)>0;
                text(distlabel(1),distlabel(2),sprintf('ep %s r %s, %s cells tr %s  ',...
                    num2str(epoch),num2str(rippnumber),num2str(cellcount),num2str(trialno)),'FontSize',3,'Color','r');
            else text(distlabel(1),distlabel(2),sprintf('ep %s r %s, %s cells tr %s',...
                    num2str(epoch),num2str(rippnumber),num2str(cellcount),num2str(trialno)),'FontSize',3,'Color','k');
                
            end
        else text(distlabel(1),distlabel(2),sprintf('ep %s r %s, %s cells tr %s ',...
                num2str(epoch),num2str(rippnumber),num2str(cellcount),num2str(trialno)),'FontSize',3,'Color','k');
        end
        
        
        ax=ax+1;
        
    end
    % Print title for whole figure
    
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
        1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 0.9,sprintf('Combined place fields of reactivated cells during ripples of %s for day %s epoch %s',...
        animals{1,1}, num2str(daylistmain(d,1)), num2str(daylistmain(d,2))),'HorizontalAlignment','center','VerticalAlignment', 'top');
    
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperSize',[30 15]);
    
    %set(gcf,'position',[10 10 2300 800]);
    
    
    
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
    
    %set(gcf, 'renderer', 'painters');
    cd(strcat(f.animal{1,2},'Plot/'));
    figurename = strcat(animals{1,1},'ripples','_',num2str(daylistmain(d,1)),'_',num2str(daylistmain(d,2)));
    
    %print('-dpdf', figurename);
    saveas(gcf, figurename, 'pdf');
    % Closes the figure
    close;
    
    
    
    
end

