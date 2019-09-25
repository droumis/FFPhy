global subplot_count;
subplot_count = 1;

Veqn = '>=0';
minV =  str2num(Veqn(end));
maxstage = 3; % [1 2 3]
minVPF = 0; %cm/sec
minPeakPF = 3;
lessthan=0;
includestates = 6;

%Animal selection
%-----------------------------------------------------
%animals = {'CML21','K3'};
%animals = {'M2','M1','M3','K3','L2','L3','N1'};
%animals = {'N1'};
animals = {'L3'};


%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filterionno

days='[1:8]';%,'1:10';
%days = '[1:1]';
%days = '[9:9]';

%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

%epochfilter{1} = [''];

epochtype='Run';

epochfilter{1} = ['isequal($epochtype, ''Run'')'];
%epochfilter{1} = ['isequal($epoch,  4)'];

%cellfilter = '(isequal($area, ''ACC'') && ($meanrate >0 ) )'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%cellfilter = '(isequal($area, ''CA1'') && ($meanrate >0 ) )'  ;

%timefilter = { {'JY_getlinvelocity', '($velocity <0.07) & ($velocity >-0.07)'} };

%timefilter = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter = { {'JY_getriptimes','($nripples == 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
%    {'JY_getlinvelocity', strcat('$velocity > ',num2str(minVPF))},{'JY_getbarrier','($barrier== 0)'}};
%timefilter = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}, {'JY_getlinvelocity', '$velocity <0.01'}};
%timefilter2 = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter2 = { {'getriptimes','($nripples > 0)',[],'minthresh',2,'cellfilter',cellfilter,'tetfilter',[1 2 3 4]}};

%timefilter = { 'JY_getlinvelocity', '$velocity > 0' };

%f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter);
%if minPeakPF>0
%    includecells = calcincludecells(minVPF, minPeakPF, animals, epochfPF, lessthan);
%5    f = excludecellsf(f, includecells);
%end
%-----------------------------------------------------------


%run function- single cells
%--------------------------------------------
iterator = 'singleepochanal';

f = setfilteriterator(f,iterator);

outall=setfilterfunction(f, 'JY_trialtrajectory', {'data','linpos'});

outall=runfilter(outall);




%plot



for mm=1:size(outall,2)
out=outall(1,mm);
uniquedays=unique(out.epochs{1,1}(:,1));
outdata=out.output{1,1};
uniquedays=unique(out.epochs{1,1}(:,1));
for ii = 1:size(uniquedays)
    trials_day = [];
    
    jind=find(out.epochs{1,1}(:,1)==uniquedays(ii,1));
        daytext=num2str(uniquedays(ii));

       for jj=1:size(jind,1)
           epochtext=num2str(out.epochs{1,1}(jind(jj),2));
        figure;

   
        
        
        data=outdata(1,jind(jj)).trajectory;
        allpos=outdata(1,jind(jj)).allpos;
        % plot subplots
        % specify paramters for tight_subplot

        dist=0;
        nrows=ceil(size(data,2)/5);
        columns=5;
        gap_h=0.01;
        gap_w=0.000001;
        marg_h=[0.01 0.1];
        marg_w=[0.01 0.2];
        ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
        ax=1;
        set(gca,'xtick',[],'ytick',[]);
            %while ax<=(nrows*columns); 
            
        for kk=1:size(data,2)
                    % get current axis

                    axes(ha(ax));
                    
                    
                    hold on;
                    
                    if ~isempty(data{1,kk}.pos)
                    
                    % plot all the positions
                    plot(allpos(:,1),allpos(:,2),'Color',[0.9 0.9 0.9]);
                    hold on;

                    %plot the start
                    plot(data{1,kk}.pos(1,2),data{1,kk}.pos(1,3),'+r');
                    
                    % plot the trajectory
                    hold on;
                    plot(data{1,kk}.pos(:,2),data{1,kk}.pos(:,3),'Color',[0 0 0]);

                    end
                  

%                         for wind=1:size(Wellpos,1);
%                             plot(Wellpos(wind,1),Wellpos(wind,2),'ok','MarkerSize',20);
%                             
%       
                    
                    %set boundries
                    xmin=min(0);
                    ymin=min(0);
                    xmax=max(allpos(:,1))*1.1;
                    ymax=max(allpos(:,2))*1.1;
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
                    %text(distlabel2(1),distlabel2(2),sprintf('%s spikes',num2str(size(data{1,kk}.spikes,1))),'FontSize',12,'Color','k');
                   
                    
%                     if outdata{1,jind(jj)}.trajectorybarrier(kk)==1
%                      % barrier label
%                      barrierlabel=[0.6*(xlimval(2)-xlimval(1))+xlimval(1) 0.10*(ylimval(2)-ylimval(1))+ylimval(1)];
%                     
%                      text(barrierlabel(1),barrierlabel(2),'barrier','FontSize',10);
%                     end
                    
                    ax=ax+1;
                end
                 
            %end
    
     
        ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
            1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
        % text(0.5, 0.98,sprintf('Ripple locations of H2 for days %s to %s \n ripples on at least %s tetrodes',...
        %  num2str(min(day)),num2str(max(day)),num2str(minrip)),'HorizontalAlignment'...
        % ,'center','VerticalAlignment', 'top');
        
      text(0.5, 0.98,sprintf('Trajectory for %s Day %s Epoch %s + start', out.animal{1,3},...
        daytext, epochtext),'FontSize',12,'HorizontalAlignment',...
         'center','VerticalAlignment', 'top');   
    
    
    % ----Saving----
    % Saves figure as pdf
    % First checks if a folder called Plot exists in the processed data folder,
    % if not, creates it.
    
    cd(out.animal{1,2});
    plotdir = dir('Plot');
    if (isempty(plotdir))
        %an a plot folder needs to be created
        !mkdir Plot
    end
    
    % change to that directory and saves the figure with file name
    % animal_day_epoch
    cd(strcat(out.animal{1,2},'Plot/'));
    figurename = strcat(out.animal{1,3},'_trialtrajectory_',daytext,'_',epochtext);
    
    saveas(gcf, figurename, 'fig');
    saveas(gcf, figurename, 'pdf');
    
    %Closes the figure
    close;
     
     
     
     
    end
end
end
       
  %% firing between trials

% for ii=1:size(uniquecells,1);
%     
%     jind=find(rowfind(daydataindex(:,3:4),uniquecells(ii,:))==1);
%    
%     hold on;
%     epochlegend=[];
%     for jj=1:size(jind,1)
%         
%         figure;
%         col= colormap(lines(size(jind,1)));
%         dayindex=find(rowfind(out.epochs{1,1},daydataindex(jind(jj),1:2))==1);
%         
%         data=outdata{1,jind(jj)}.intertrial;
%         allpos=outdata{1,jind(jj)}.allpos;
%         % plot subplots
%         % specify paramters for tight_subplot
% 
%         dist=0;
%         nrows=ceil(size(data,2)/5);
%         columns=5;
%         gap_h=0.01;
%         gap_w=0.000001;
%         marg_h=[0.01 0.1];
%         marg_w=[0.01 0.2];
%         ha = tight_subplot(nrows, columns, [gap_h gap_w], marg_h,marg_w);
%         ax=1;
% 
%             %while ax<=(nrows*columns); 
%             
%         for kk=1:size(data,2)
%                     % get current axis
% 
%                     axes(ha(ax));
%                     
%                     
%                     hold on;
%      
%                     
%                     % plot all the positions
%                     plot(allpos(:,1),allpos(:,2),'Color',[0.9 0.9 0.9]);
%                     hold on;
%                     if ~isempty(data{1,kk}.pos)
%                     %plot the start
%                     plot(data{1,kk}.pos(1,2),data{1,kk}.pos(1,3),'+g');
%                     
%                     % plot the trajectory
%                     hold on;
%                     plot(data{1,kk}.pos(:,2),data{1,kk}.pos(:,3),'Color',[0.5 0.5 0.5]);
% 
%                     hold on;
%                     % plot the spikes
%                     plot(data{1,kk}.spikes(:,2),data{1,kk}.spikes(:,3),'.r');
%                     
%                     end
% 
% %                         for wind=1:size(Wellpos,1);
% %                             plot(Wellpos(wind,1),Wellpos(wind,2),'ok','MarkerSize',20);
% %                             
% %       
%                     
%                     %set boundries
%                     xmin=min(0);
%                     ymin=min(0);
%                     xmax=160;
%                     ymax=120;
%                     fac=0.2;
%                     xmarg=fac*(xmax-xmin);
%                     ymarg=fac*(ymax-ymin);
%                     
%                     xlim([xmin-xmarg xmax+xmarg]);
%                     ylim([ymin-ymarg ymax+ymarg]);
%                     set(gca,'xtick',[],'ytick',[]);
%                     set(gca,'PlotBoxAspectRatio',[1 1 1]);
%                     
%                     % print text
%                     xlimval=xlim;
%                     ylimval=ylim;
%                     % distlabel=[0.35*(xlimval(2)-xlimval(1))+xlimval(1) 0.05*(ylimval(2)-ylimval(1))+ylimval(1)];
%                     % text(distlabel(1),distlabel(2),sprintf('day %s epoch %s',num2str(d),num2str(e)),'FontSize',12,'Color','k');
%                     % print number of ripples
%                     distlabel2=[0.10*(xlimval(2)-xlimval(1))+xlimval(1) 0.85*(ylimval(2)-ylimval(1))+ylimval(1)];
%                     text(distlabel2(1),distlabel2(2),sprintf('%s spikes',num2str(size(data{1,kk}.spikes,1))),'FontSize',12,'Color','k');
%                    
%                     
%                     if outdata{1,jind(jj)}.trajectorybarrier(kk)==1
%                      % barrier label
%                      barrierlabel=[0.6*(xlimval(2)-xlimval(1))+xlimval(1) 0.10*(ylimval(2)-ylimval(1))+ylimval(1)];
%                     
%                      text(barrierlabel(1),barrierlabel(2),'barrier','FontSize',10);
%                     end
%                     
%                     ax=ax+1;
%                 end
%                  
%             %end
%     
%      
%         ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
%             1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%         % text(0.5, 0.98,sprintf('Ripple locations of H2 for days %s to %s \n ripples on at least %s tetrodes',...
%         %  num2str(min(day)),num2str(max(day)),num2str(minrip)),'HorizontalAlignment'...
%         % ,'center','VerticalAlignment', 'top');
%         
%       text(0.5, 0.98,sprintf('Spike location between trials for %s \n Day %s Epoch %s Tetrode %s cell %s', animals{1,1},...
%         num2str(unique(daydataindex(:,1))), num2str(dayindex), num2str(uniquecells(ii,1)),num2str(uniquecells(ii,2))),'FontSize',12,'HorizontalAlignment',...
%          'center','VerticalAlignment', 'top');   
%     
%     
%     % ----Saving----
%     % Saves figure as pdf
%     % First checks if a folder called Plot exists in the processed data folder,
%     % if not, creates it.
%     
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
%     figurename = strcat(animals{1,1},'_intertrialspikelocation_d',num2str(daytext),'_e',num2str(dayindex),'_t',num2str(uniquecells(ii,1)),'_c',num2str(uniquecells(ii,2)));
%     
%     saveas(gcf, figurename, 'pdf');
%     
%     %Closes the figure
%     close;
%      
%      
%      
%      
%     end
% end  
    
