day=[1];
epoch=[1:7];
tet=[1:14];
minrip=7;

% plot subplots
% specify paramters for tight_subplot
figure;
dist=0;
nrows=size(day,2);
columns=size(epoch,2);
gap_h=0.01;
gap_w=0.000001;
marg_h=[0.01 0.1];
marg_w=[0.01 0.2];
ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
ax=1;


while ax<=(nrows*columns);
    
    for d=day;
        % open posdata file
            % Load the mat file
            % See if day number needs 0
            dsz = '';
            if (d < 10)
                dsz = '0';
            end
            sfilename = strcat('/data14/jai/H2_/','H2_',dsz,num2str(d),'.mat');
        load(sfilename);
        % load posfile
        posfilename = strcat('/data14/jai/H2_/','H2_',dsz,num2str(d),'linpos.mat');
        load(posfilename);
        for e=epoch;
            
            % gets time when more than minrip tetrodes have ripples
            riptimes =JY_getriptimes('/data14/jai/H2_/','H2',[d e],tet,minrip);
            
            
            
            
            % find times with multiple ripples
            tind=lookup(riptimes(:,1),Data{1,d}{1,e}.Pos.correcteddata(:,1));
            tind=unique(tind);
            % plots where ripples happen
            posx=Data{1,d}{1,e}.Pos.correcteddata(tind,2);
            posy=Data{1,d}{1,e}.Pos.correcteddata(tind,3);
            
            % plots all trajectory
            posallx=Data{1,d}{1,e}.Pos.correcteddata(:,2);
            posally=Data{1,d}{1,e}.Pos.correcteddata(:,3);
            
            
            % get current axis
            axes(ha(ax));
            
            
            plot(posallx,posally,'.',...
                'MarkerEdgeColor',[0.5 0.5 0.5],...
                'MarkerFaceColor',[0.5 0.5 0.5],...
                'MarkerSize',10);
            hold on;
            plot(posx,posy,'.r','MarkerSize',10);
            hold on;
            
            
            if strcmp(Data{1,d}{1,e}.Stats.epochtype,'Sleep')~=1;
                % plot circle round rewarded wells
                % get rewarded wells
                Wells=Data{1,d}{1,e}.Wellinfo.rewardedwells;
                % get reward well positions
                Wellpos=linpos{1,d}{1,e}.wellSegmentInfo.wellCoord(Wells,:);

                % plot circles around wells
                
                for wind=1:size(Wellpos,1);
                    plot(Wellpos(wind,1),Wellpos(wind,2),'ok','MarkerSize',20);
                    
                end
            end
            
            %set boundries
            xmin=min(posallx);
            ymin=min(posally);
            xmax=max(posallx);
            ymax=max(posally);
            fac=0.2;
            xmarg=fac*(xmax-xmin);
            ymarg=fac*(ymax-ymin);
            
            xlim([xmin-xmarg xmax+xmarg]);
            ylim([ymin-ymarg ymax+ymarg]);
            set(gca,'xtick',[],'ytick',[]);
            set(gca,'PlotBoxAspectRatio',[1 1 1]);
            
            % print text
            xlimval=xlim;
            ylimval=ylim;
            % distlabel=[0.35*(xlimval(2)-xlimval(1))+xlimval(1) 0.05*(ylimval(2)-ylimval(1))+ylimval(1)];
            % text(distlabel(1),distlabel(2),sprintf('day %s epoch %s',num2str(d),num2str(e)),'FontSize',12,'Color','k');
            % print number of ripples
            distlabel2=[0.10*(xlimval(2)-xlimval(1))+xlimval(1) 0.85*(ylimval(2)-ylimval(1))+ylimval(1)];
            text(distlabel2(1),distlabel2(2),sprintf('%s ripples',num2str(size(tind,1))),'FontSize',12,'Color','k');
            ax=ax+1;
        end
    end
end

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
    1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.5, 0.98,sprintf('Ripple locations of H2 for days %s to %s \n ripples on at least %s tetrodes',...
  %  num2str(min(day)),num2str(max(day)),num2str(minrip)),'HorizontalAlignment'...
   % ,'center','VerticalAlignment', 'top');




