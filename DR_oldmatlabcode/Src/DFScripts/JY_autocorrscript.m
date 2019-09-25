% plots cross correlograms for co-firing of cell parts under different
% conditions

% compare ripples during stillness and runs
% with barrier and without

% r_b run with barrier
% r_nb run no barrier
% s_b stillness with barrier
% s_nb stillness without barrier



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
animals = {'CML21'};
%-----------------------------------------------------

%Filter creation
%--------------------------------------------------------
% day filterionno


days = '[2]';

epochfilter = ['isequal($epochtype, ''Run'')'];
%epochfilter{1} = ['isequal($epoch, 4)'];

%cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 10) && ($numspikes>100))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%timefilter = { {'JY_getlinvelocity', '($velocity <0.07) & ($velocity >-0.07)'} };

cellfilter = '(isequal($area, ''ACC'') )'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%timefilter = { {'JY_getlinvelocity', '($velocity <0.07) & ($velocity >-0.07)'} };




%timefilter = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
 %   {'JY_getlinvelocity', strcat('abs($velocity > ',num2str(minVPF),')')}};



%------------------------------------------------------

 f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter);

%-------------------------------------------------------
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
 f = setfilterfunction(f, 'JY_spkautocorr', {'spikes'}, 0.5, includestates, minV, f);
 f = runfilter(f);


%day=f_s_nb.epochs{1,1}(1,1);
day=f.epochs{1,1}(1,1);


% normalise


% generate index to show day epoch and tetrode channel
daydataindex=[];



for ii=1:size(f.epochs{1,1},1)
    tempind=f.data{1,1}{1,ii};
    tempdayind=repmat(f.epochs{1,1}(ii,:),size(tempind,1),1);
    daydataindex=[daydataindex;tempdayind tempind];
end

daytext=unique(daydataindex(:,1));

% find unique cells

uniquecells=unique(daydataindex(:,3:4),'rows');

% plot unique cells

% compare firing rate to trajectory distance

for ii=1:size(uniquecells,1);
    
    uniquecellindex=find(rowfind(daydataindex(:,3:4),uniquecells(ii,:))==1);
    figure;
    hold on;
    epochlegend=[];
    col= colormap(lines(size(uniquecellindex,1)));
    for jj=1:size(uniquecellindex,1)
        
        
        

        plot(f.output{1,1}{1,uniquecellindex(jj)}.time, f.output{1,1}{1,uniquecellindex(jj)}.convautocorr,'Color',col(jj,:));
        set(gca,'YScale','log')
        %set(gca,'YScale','linear')
        

    end
    
      title(sprintf('Autocorrelogram for %s day %s \n tet %d cell %d', animals{1,1},num2str(day), ...
        uniquecells(ii,1),uniquecells(ii,2)));
    
    legend('2','4','6','Location','EastOutside');
    
    set(gca,'XTick',-0.5:0.1:0.5)
    
      % change to that directory and saves the figure with file name
    % animal_day_epoch
    cd(strcat(f.animal{1,2},'Plot/'));
    figurename = sprintf('autocorr_%s_d%d_t%d_c%d',animals{1,1},day,...
        uniquecells(ii,1),uniquecells(ii,2));
    
    saveas(gcf, figurename, 'pdf');
    close;
    
end

        
        
        
%         
%         
%         
%         
%         
%         
% % gaussian smoothing
% sigma=1;
% width = round((6*sigma - 1)/2);
% support = (-width:width);
% gaussFilter = exp( -(support).^2 ./ (2*sigma^2) );
% gaussFilter = gaussFilter/ sum(gaussFilter);
% 
% win=gausswin(64,1);
% win=win./sum(win);
% %convc1vsc2{ii,jj}=conv(spike_xcorr.c1vsc2,win);
% %convc1vsc2=[convc1vsc2;conv(spike_xcorr.c1vsc2,gaussFilter,'same')];
% 
% 
% for i=1:size(f.output{3,1},1)
%     figure;
%     hold on;
%     
%     plot(f.output{1,1}(1,2).time(1,:),conv(f_s_nb.output{4,1}(i,:),gaussFilter,'same'),'-r'); % second epoch has time references
%   
%     set(gca,'XTick',-0.5:0.1:0.5);
%     title(sprintf('Autocorrelogram for %s day %s \n tet %d cell %d vs tet %d vs %d', animals{1,1},num2str(day), ...
%         commoncells(i,1),commoncells(i,2),commoncells(i,3),commoncells(i,4)));
%     legend(sprintf('< %s -barrier',num2str(minVPF)),sprintf('< %s +barrier',num2str(minVPF)),...
%         sprintf('> %s -barrier',num2str(minVPF)),sprintf('> %s +barrier',num2str(minVPF)));
%     
%     % change to that directory and saves the figure with file name
%     % animal_day_epoch
%     cd(strcat(f_r_b.animal{1,2},'Plot/'));
%     figurename = sprintf('CCG_%s_%d_%d_%dv%d_%d',animals{1,1},day,...
%         commoncells(i,1),commoncells(i,2),commoncells(i,3),commoncells(i,4));
%     
%     saveas(gcf, figurename, 'pdf');
%     close;
%     
%     
%     
%     
% end
% 
% 
% 
% figure;
% boxplot([f_s_b.output{3,1}(:,1);f_s_nb.output{3,1}(:,1);f_r_b.output{3,1}(:,1);f_r_nb.output{3,1}(:,1)],...
%     [ones(size(commoncells,1),1).*1;ones(size(commoncells,1),1).*2;ones(size(commoncells,1),1).*3;ones(size(commoncells,1),1).*4],...
%     'labels',{'<2cm/s +barrier','<2cm/s -barrier','>2cm/s +barrier','>2cm/s -barrier'},'labelorientation','inline');
% title(sprintf('Spike number for each group for %s day %s', animals{1,1},num2str(day)));
% 
% 
% cd(strcat(f_r_b.animal{1,2},'Plot/'));
% figurename = sprintf('spikesummary_%s_%d',animals{1,1},day);
% 
% saveas(gcf, figurename, 'pdf');
% close;
% 
% 












%  day=f_r_b.epochs{1,1}(1,1);
%  epoch=f_r_b.epochs{1,1}(1,2);
% for ii=1:size(f_s_nb.output{1,1}.convc1vsc2,1)
%     for jj=1:size(f_s_nb.output{1,1}.convc1vsc2,2);
%         if ii ~= jj
%             figure;
%             hold on;
%             plot(f_s_nb.output{1,1}.time{ii,jj},f_s_nb.output{1,1}.convc1vsc2{ii,jj},'-r');
%             plot(f_s_b.output{1,1}.time{ii,jj},f_s_b.output{1,1}.convc1vsc2{ii,jj},'--r');
%             plot(f_r_nb.output{1,1}.time{ii,jj},f_r_nb.output{1,1}.convc1vsc2{ii,jj},'-b');
%             plot(f_r_b.output{1,1}.time{ii,jj},f_r_b.output{1,1}.convc1vsc2{ii,jj},'--b');
%
%             set(gca,'XTick',-0.5:0.1:0.5);
%
%             title(sprintf('Cross correlogram for %s day %s epoch %s \n tet %s cell %s vs tet %s vs %s', animals{1,1},num2str(day),num2str(epoch), ...
%                 num2str(f_r_b.data{1,1}{1,1}(ii,1)), num2str(f_r_b.data{1,1}{1,1}(ii,2)),...
%                 num2str(f_r_b.data{1,1}{1,1}(jj,1)), num2str(f_r_b.data{1,1}{1,1}(jj,2))));
%             legend(sprintf('< %s -barrier',num2str(minVPF)),sprintf('< %s +barrier',num2str(minVPF)),...
%                 sprintf('> %s -barrier',num2str(minVPF)),sprintf('> %s +barrier',num2str(minVPF)));
%
%
%
%                 % change to that directory and saves the figure with file name
%                 % animal_day_epoch
%                 cd(strcat(f_r_b.animal{1,2},'Plot/'));
%                 figurename = sprintf('CCG_%s_%d_%d_%d_%dv%d_%d',animals{1,1},day,epoch,...
%                     f_r_b.data{1,1}{1,1}(ii,1), f_r_b.data{1,1}{1,1}(ii,2),...
%                 f_r_b.data{1,1}{1,1}(jj,1), f_r_b.data{1,1}{1,1}(jj,2));
%
%                 saveas(gcf, figurename, 'pdf');
%                 close;
%
%
%
%
%         end
%     end
% end
%
% % figure;
% % imagesc(f.output{1,1});
