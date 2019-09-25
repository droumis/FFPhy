function [ outputdata ] = PlotReroutetraj_epochsdays( animaldir,animalname,daylist )
% Plots the spatial occupancy of each epoch for days in daylist
% datadirectory - directory where the output from batchDIOextract is
% located, eg '/data14/jai/F12_/'
% animalname -  name of the animal, eg F1
% days - vector of day or days to be analysed, eg [1], [1:3]
%


% get the current directory to go back to after the function is done
currentdir = pwd;

% ---- File loading ----

% Specify data directory and load the file
%animalname = animaldir(1:end-1);
datadir = '/data14/jai/';

% specify paramters for tight_subplot
figure;
dist=0;
%nrows=size(daylist,2);

nrows=3;
%columns=3; %7 to include sleep sessions

columns=size(daylist,2);
gap_h=0.005;
gap_w=0.005;
marg_h=[0.01 0.1];
marg_w=[0.01 0.2];
ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
ax=1;

% generate colormap

% cmap = colormap(grey(nc));
% % adjusts scaling of colormap
% factor=1;
% %factor=5*log(cmap+0.9);
% cmap=cmap.*factor;
% cmap(cmap<=0)=0;
% %cmap(cmap>=1)=1;
% cmap(1,:) = 0;
cmap=colormap(hot(128));
%cmap=colormap(jet(128));

% loop for each day
while ax<(nrows*columns);
    
    
    
    for d=1:size(daylist,2);
        % See if day number needs 0
        dsz = '';
        if (daylist(d) < 10)
            dsz = '0';
        end
        
        outputdata={};
        
        % Converts the day and epoch into text
        dayt = num2str(daylist(1,d));
        % Load the mat file
        sfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'data',dsz,dayt,'.mat');
        load(sfilename);
        % load posfile
        posfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'linpos',dsz,dayt,'.mat');
        load(posfilename);
        dayepochlist=1;
        %for k=1:size(data{1,daylist(1,d)},2);
        for k=1:6;
            if isempty(data{1,daylist(1,d)}{1,k});
                %ax=ax+1;
                k=k+1;
                dayepochlist=dayepochlist+1;
                
                
            else
                if strcmp(data{1,daylist(1,d)}{1,k}.Stats.animal,animalname)==1;
                    if strcmp(data{1,daylist(1,d)}{1,k}.Stats.epochtype,'Run')==1;
                        
                        epoch = k;
                        epocht = num2str(epoch);
                        % get all positions for epoch
                        cmperpix=data{1,daylist(1,d)}{1,k}.Pos.cmperpixel;
                        traj=data{1,daylist(1,d)}{1,k}.Pos.correcteddata(:,2:3).*cmperpix;
                        
                        wellfactor=1;
                        %                 if max(traj(:,1)>=200)
                        %                     traj=traj./2;
                        %                     wellfactor=2;
                        %                 end
                        
                        % scale factor for cutting off upper occupancy
                        scf=0.0003;
                        % generate smoothed occupancy
                        binsize=2; % in cm
                        edge=15;
                        %                 minx = floor(min(traj(:,1)))-edge;
                        %                 maxx = ceil(max(traj(:,1)))+edge;
                        %
                        minx = 90
                        maxx = 300;
                        
                        
                        
                        binx = (minx:binsize:maxx);
                        
                        %miny = floor(min(traj(:,2)))-edge;
                        %maxy = ceil(max(traj(:,2)))+edge;
                        miny = 5;
                        maxy = 220;
                        biny = (miny:binsize:maxy);
                        
                        
                        
                        [output.occupancy output.xticks output.yticks] = hist2(traj(:,1),traj(:,2), binx, biny);
                        maxv=sum(sum(output.occupancy));
                        output.occupancy=output.occupancy*100/maxv;
                        % get quantile
                        linlist=output.occupancy(output.occupancy~=0);
                        linlist=sort(linlist);
                        qscf=ceil(scf*size(linlist,1));
                        
                        outputdata{daylist(1,d)}{k}=output;
                        
                        
                        % set gaussian parameters
                        std=2;
                        g = gaussian2(std,(6*std));
                        output.occupancy = filter2(g,(output.occupancy));
                        linlist2=output.occupancy(output.occupancy~=0);
                        linlist2=sort(linlist2);
                        qscf2=ceil(scf*size(linlist2,1));
                        output.occupancy(output.occupancy>=scf*sum(sum(output.occupancy)))=scf*sum(sum(output.occupancy));
                        
                        % get current axis
                        
                        ax=(k/2)*columns-columns+d;
                        
                        axes(ha(ax));
                        
                        % plot image
                        imagedata = flipud(output.occupancy);
                        %image(imagedata,'CDataMapping', 'scaled');
                        imagesc(imagedata);
                        %image(imagedata);
                        
                        if strcmp(data{1,daylist(1,d)}{1,k}.Stats.epochtype,'Sleep')~=1;
                            % plot circle round rewarded wells
                            % get rewarded wells
                            Wells=data{1,daylist(1,d)}{1,k}.Wellinfo.rewardedwells;
                            
                            % get reward well positions
                            if ~isempty(linpos{1,daylist(1,d)}{1,k});
                                Wellpos=(linpos{1,daylist(1,d)}{1,k}.wellSegmentInfo.wellCoord./wellfactor).*cmperpix;
                                [output.occupancy output.xticks output.yticks]=hist2(Wellpos(:,1),Wellpos(:,2), binx, biny);
                                
                                % plot circles around wells
                                hold on;
                                output.occupancy = flipud(output.occupancy);
                                [wella, wellb]=find(output.occupancy==1);
                                Wellpos=[wella,wellb];
                                
                                for wind=1:size(Wellpos,1);
                                    plot(Wellpos(wind,2),Wellpos(wind,1),'og','MarkerSize',15);
                                    
                                end
                            end
                            
                        end
                        
                        xlimval=[0 size(output.occupancy,2)];
                        ylimval=[0 size(output.occupancy,1)];
                        
                        distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.95*(ylimval(2)-ylimval(1))+ylimval(1)];
                        %titlelabel=[0.8*(xlimval(2)-xlimval(1))+xlimval(1) 0.3*(ylimval(2)-ylimval(1))+ylimval(1)];
                        
                        %text(distlabel(1),distlabel(2),sprintf('day %s epoch %s %s s to %s s',num2str(daylist(1,d)),num2str(k),num2str(scf*0.033*maxv,2),num2str(0.033*maxv,3)),'FontSize',8,'Color','w');
                        %text(titlelabel(1),titlelabel(2),sprintf('day %s',num2str(daylist(1,d))),'FontSize',8,'Color','w');
                        
                        set(gca,'xtick',[],'ytick',[]);
                        set(gca,'PlotBoxAspectRatio',[1 1 1]);
                        
                        if dayepochlist==1 && k==size(data{1,daylist(1,d)},2)
                            ax=ax+2;
                        end
                        
                        if dayepochlist==2 && k==size(data{1,daylist(1,d)},2)
                            ax=ax+1;
                        end
                        
                        dayepochlist=dayepochlist+1;
                        k=k+1;
                        ax=ax+1;
                        
                    else display(sprintf('No %s data in day %s epoch %s ',animalname,num2str(daylist(1,d)),num2str(k)));
                        if dayepochlist==2 && k==size(data{1,daylist(1,d)},2)
                            ax=ax+2;
                        end
                        if dayepochlist==3 && k==size(data{1,daylist(1,d)},2)
                            ax=ax+1;
                        end
                        k=k+1;
                        
                    end
                end
            end
        end
        
        occupancyoutput=outputdata;
        
        directoryname=sprintf('/data14/jai/%s_/',animalname');
        
        cd(directoryname);

        save(sprintf('%soccupancyoutput%s%s.mat',animalname,dsz,dayt),'occupancyoutput');

    end
end

% Print title for whole figure

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 0.98,sprintf('Occupancy of %s for days %s to %s',animalname, num2str(min(daylist)),num2str(max(daylist))),'HorizontalAlignment'...
    ,'center','VerticalAlignment', 'top');
step=scf/6;
colorbarlabel=[0:step:scf]*100;
hcol=colorbar('YTickLabel', num2str(colorbarlabel(2:end)',2),'Location','EastOutside');
%hcol=colorbar('Location','EastOutside');
cpos=get(hcol,'Position');
cpos(4)=cpos(4)/2;      % Halve the thickness
cpos(1)=cpos(1)-0.05;
cpos(2)=cpos(2)+0.3;   %  Move it down outside the plot
set(hcol,...
    'Position',cpos)

set(gcf,'PaperPositionMode','auto')
set(gcf,'Position',[1145 457 1375 450])
set(gcf,'PaperSize',[20 8.5])


% % saving
%directoryname='/home/jai/Documents/Presentations/Labmeetings/20130121/';
directoryname='/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/';

print(sprintf('%s%s_occupancyafterinactivation',directoryname,animalname),'-depsc');

close;




end