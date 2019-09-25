day=[6];
epoch=[2];
tet=[1:14];
minrip=7;
twindow=100; % time window in ms to look at ripples

% plot subplots
% specify paramters for tight_subplot
%figure;

rippletable={};

%while ax<=(nrows*columns);

for d=day;
    for e=epoch;
        
        % gets time when more than minrip tetrodes have ripples
        riptimes =getriptimes('/data14/jai/H2_/','H2',[d e],tet,minrip);
        
        % open posdata file
        % Load the mat file
        % See if day number needs 0
        dsz = '';
        if (d < 10)
            dsz = '0';
        end
        
        sfilename = strcat('/data14/jai/H2_/','H2_',dsz,num2str(d),'.mat');
        load(sfilename);
        % find times with multiple ripples
        tind=lookup(riptimes(:,1),Data{1,d}{1,e}.Pos.correcteddata(:,1));
        tind=unique(tind);
        % translate tind into time
        rpt=Data{1,d}{1,e}.Pos.correcteddata(tind,1);
        % lookup ripples
        
        for tetind=7;
            tsz = '';
            if (tetind < 10)
                tsz = '0';
            end
            %open  eeg file
            tetfile=strcat('/data14/jai/H2_/EEG/','H2eeg',dsz,num2str(d),'-',num2str(e),'-',tsz,num2str(tetind),'.mat');
            load(tetfile);
            
            
            
            %get index of ripple start time in EEG file
            fieldref=eeg{1,d}{1,e}{1,tetind};
            rpind=(rpt-fieldref.starttime)./(1/fieldref.samprate);
            rpind=round(rpind);
            
            dist=0;
            nrows=round(size(rpind,1)/10)+1;
            columns=10; %size(epoch,2);
            gap_h=0.01;
            gap_w=0.000001;
            marg_h=[0.01 0.1];
            marg_w=[0.01 0.2];
            ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
            ax=1;
            
            
            
            
            % translate time window into number of cells
            tcells=round((twindow/1000)*fieldref.samprate);
            
            % make table
            ript=[];
            
            % get eeg in window after start time
            
            for rind=1:size(rpind,1);
                currript=fieldref.data(rpind(rind,1):rpind(rind,1)+tcells,1)';
                ript=[ript;currript];
                rind=rind+1;
                
                axes(ha(ax));
                plot(currript,'k');
                ax=ax+1;

            end
            rippletable{1,tetind}=ript;
            
            % plot eeg
            
            
            
        end
    end
end






%
%
%             % plots where ripples happen
%             posx=Data{1,d}{1,e}.Pos.correcteddata(tind,2);
%             posy=Data{1,d}{1,e}.Pos.correcteddata(tind,3);
%
%             % plots all trajectory
%             posallx=Data{1,d}{1,e}.Pos.correcteddata(:,2);
%             posally=Data{1,d}{1,e}.Pos.correcteddata(:,3);
%
%
%             % get current axis
%             axes(ha(ax));
%
%
%             plot(posallx,posally,'.',...
%                 'MarkerEdgeColor',[0.5 0.5 0.5],...
%                 'MarkerFaceColor',[0.5 0.5 0.5],...
%                 'MarkerSize',10);
%             hold on;
%             plot(posx,posy,'.r','MarkerSize',10);
%
%
%
%
%             %set boundries
%             xmin=min(posallx);
%             ymin=min(posally);
%             xmax=max(posallx);
%             ymax=max(posally);
%             fac=0.2;
%             xmarg=fac*(xmax-xmin);
%             ymarg=fac*(ymax-ymin);
%
%             xlim([xmin-xmarg xmax+xmarg]);
%             ylim([ymin-ymarg ymax+ymarg]);
%             set(gca,'xtick',[],'ytick',[]);
%             set(gca,'PlotBoxAspectRatio',[1 1 1]);
%
%             % print text
%             xlimval=xlim;
%             ylimval=ylim;
%             distlabel=[0.35*(xlimval(2)-xlimval(1))+xlimval(1) 0.05*(ylimval(2)-ylimval(1))+ylimval(1)];
%             text(distlabel(1),distlabel(2),sprintf('day %s epoch %s',num2str(d),num2str(e)),'FontSize',12,'Color','k');
%             % print number of ripples
%             distlabel2=[0.10*(xlimval(2)-xlimval(1))+xlimval(1) 0.85*(ylimval(2)-ylimval(1))+ylimval(1)];
%             text(distlabel2(1),distlabel2(2),sprintf('%s ripples',num2str(size(tind,1))),'FontSize',12,'Color','k');
%             ax=ax+1;
%         end
%     end
% %end
%
% ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
%     1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
% text(0.5, 0.98,sprintf('Ripple locations of H2 for days %s to %s \n ripples on at least %s tetrodes',...
%     num2str(min(day)),num2str(max(day)),num2str(minrip)),'HorizontalAlignment'...
%     ,'center','VerticalAlignment', 'top');
%



