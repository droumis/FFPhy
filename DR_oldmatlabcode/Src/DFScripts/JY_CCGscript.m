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


days = '[4]';

epochfilter = ['isequal($epochtype, ''Run'')'];
%epochfilter{1} = ['isequal($epoch, 4)'];

%cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 10) && ($numspikes>100))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%timefilter = { {'JY_getlinvelocity', '($velocity <0.07) & ($velocity >-0.07)'} };

cellfilter = '(isequal($area, ''ACC'') )'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%timefilter = { {'JY_getlinvelocity', '($velocity <0.07) & ($velocity >-0.07)'} };

% timefilter_r_b = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
%     {'JY_getlinvelocity', strcat('abs($velocity > ',num2str(minVPF),')')},...
%     {'JY_getbarrier','($barrier== 1)'}}; % barrier=0 means get all non-barrier events
% 
% timefilter_r_nb = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
%     {'JY_getlinvelocity', strcat('abs($velocity > ',num2str(minVPF),')')},...
%     {'JY_getbarrier','($barrier== 0)'}}; % barrier=0 means get all non-barrier events
% 
% timefilter_s_b = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
%     {'JY_getlinvelocity', strcat('abs($velocity < ',num2str(minVPF),')')},...
%     {'JY_getbarrier','($barrier== 1)'}}; % barrier=0 means get all non-barrier events
% 
% timefilter_s_nb = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
%     {'JY_getlinvelocity', strcat('abs($velocity < ',num2str(minVPF),')')},...
%     {'JY_getbarrier','($barrier== 0)'}}; % barrier=0 means get all non-barrier events


%timefilter = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
 %   {'JY_getlinvelocity', strcat('abs($velocity > ',num2str(minVPF),')')}};



%------------------------------------------------------
% f_r_b = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter_r_b);
% f_r_nb = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter_r_nb);
% f_s_b = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter_s_b);
% f_s_nb = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter_s_nb);
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
iterator = 'multicellanal';

% f_r_b = setfilteriterator(f_r_b,iterator);
% f_r_nb = setfilteriterator(f_r_nb,iterator);
% f_s_b = setfilteriterator(f_s_b,iterator);
% f_s_nb = setfilteriterator(f_s_nb,iterator);
% 
% 
% f_r_b = setfilterfunction(f_r_b, 'JY_spkxcorr', {'spikes'}, 0.5, includestates, minV, f_r_b);
% f_r_nb = setfilterfunction(f_r_nb, 'JY_spkxcorr', {'spikes'}, 0.5, includestates, minV, f_r_nb);
% f_s_b = setfilterfunction(f_s_b, 'JY_spkxcorr', {'spikes'}, 0.5, includestates, minV, f_s_b);
% f_s_nb = setfilterfunction(f_s_nb, 'JY_spkxcorr', {'spikes'}, 0.5, includestates, minV, f_s_nb);
% 
% f_r_b = runfilter(f_r_b);
% f_r_nb = runfilter(f_r_nb);
% f_s_b = runfilter(f_s_b);
% f_s_nb = runfilter(f_s_nb);
% 
% f_r_b.output{2,1}={};
% f_r_nb.output{2,1}={};
% f_s_b.output{2,1}={};
% f_s_nb.output{2,1}={};

 f = setfilteriterator(f,iterator);
 f = setfilterfunction(f, 'JY_spkautocorr', {'spikes'}, 0.2, includestates, minV, f);
 f = runfilter(f);


%day=f_s_nb.epochs{1,1}(1,1);
day=f{1,1}(1,1);



plot(f.output{1,1}(1,1).time(3,:), f.output{1,1}(1,1).convautocorr(3,:))
set(gca,'YScale','log')
set(gca,'YScale','linear')

% for i=1:1:size(f_r_b.epochs{1,1},1);
%
%
%
%
%     sum=arrayfun(@(x)sum(x),
%
%
% % get all epochs together
% for ii=1:size(f_s_nb.output{1,1}.(i,i)convc1vsc2,1)
%     for jj=1:size(f_s_nb.output{1,1}.convc1vsc2,2);
%         if ii ~= jj
%             currcell=
% f_r_b.output{2,1}=

% 
% commoncells=mintersect(f_r_nb.output{1,1}(1,1).cellpairindex,...
%     f_r_nb.output{1,1}(1,2).cellpairindex,...
%     f_r_nb.output{1,1}(1,3).cellpairindex,...
%     f_r_b.output{1,1}(1,2).cellpairindex,...
%     f_r_b.output{1,1}(1,3).cellpairindex,...
%     f_s_nb.output{1,1}(1,1).cellpairindex,...
%     f_s_nb.output{1,1}(1,2).cellpairindex,...
%     f_s_nb.output{1,1}(1,3).cellpairindex,...
%     f_s_b.output{1,1}(1,2).cellpairindex,...
%     f_s_b.output{1,1}(1,3).cellpairindex,...
%     'rows');



%                         f_s_nb.output{1,1}(1,1).cellpairindex,...fgsfg
%                         f_s_nb.output{1,1}(1,2).cellpairindex,...
%                         f_s_nb.output{1,1}(1,3).cellpairindex,...
%                         f_s_b.output{1,1}(1,1).cellpairindex,...
%                         f_s_b.output{1,1}(1,2).cellpairindex,...
%                         f_s_b.output{1,1}(1,3).cellpairindex,...
%f_r_b.output{1,1}(1,1).cellpairindex,...

output=zeros(size(commoncells,1), size(f_s_nb.output{1,1}(1,1).time,2));

f_s_nb.output{2,1}=output;
f_s_b.output{2,1}=output;
f_r_nb.output{2,1}=output;
f_r_b.output{2,1}=output;
f_s_nb.output{3,1}=output(:,1);
f_s_b.output{3,1}=output(:,1);
f_r_nb.output{3,1}=output(:,1);
f_r_b.output{3,1}=output(:,1);
for i=1:1:size(f_r_b.epochs{1,1},1);
    
    %     %[c,diffindex]=setdiff(f_r_nb.output{1,1}(1,1).cellpairindex,commoncells ,'rows');
    %
    %     if ~isempty(diffindex)
    
    goodcellindex1=ismember(f_s_nb.output{1,1}(1,i).cellpairindex,commoncells,'rows');
    goodcellindex2=ismember(f_s_b.output{1,1}(1,i).cellpairindex,commoncells,'rows');
    goodcellindex3=ismember(f_r_nb.output{1,1}(1,i).cellpairindex,commoncells,'rows');
    goodcellindex4=ismember(f_r_b.output{1,1}(1,i).cellpairindex,commoncells,'rows');
    if ~isempty(goodcellindex1)
        
        f_s_nb.output{2,1}=f_s_nb.output{2,1}+f_s_nb.output{1,1}(1,i).unnorm(goodcellindex1,:);
        f_s_nb.output{3,1}=f_s_nb.output{3,1}+f_s_nb.output{1,1}(1,i).normalisenumber(goodcellindex1,:);
    end
    if ~isempty(goodcellindex2)
        
        f_s_b.output{2,1}=f_s_b.output{2,1}+f_s_b.output{1,1}(1,i).unnorm(goodcellindex2,:);
        f_s_b.output{3,1}=f_s_b.output{3,1}+f_s_b.output{1,1}(1,i).normalisenumber(goodcellindex2,:);
        
    end
    if ~isempty(goodcellindex3)
        f_r_nb.output{2,1}=f_r_nb.output{2,1}+f_r_nb.output{1,1}(1,i).unnorm(goodcellindex3,:);
        f_r_nb.output{3,1}=f_r_nb.output{3,1}+f_r_nb.output{1,1}(1,i).normalisenumber(goodcellindex3,:);
        
        
    end
    if ~isempty(goodcellindex4)
        f_r_b.output{2,1}=f_r_b.output{2,1}+f_r_b.output{1,1}(1,i).unnorm(goodcellindex4,:);
        f_r_b.output{3,1}=f_r_b.output{3,1}+f_r_b.output{1,1}(1,i).normalisenumber(goodcellindex4,:);
        
        
    end
    
end

% normalise



f_s_nb.output{3,1}=repmat(f_s_nb.output{3,1},1,size(f_s_nb.output{2,1},2));
f_s_b.output{3,1}=repmat(f_s_b.output{3,1},1,size(f_s_b.output{2,1},2));
f_r_nb.output{3,1}=repmat(f_r_nb.output{3,1},1,size(f_r_nb.output{2,1},2));
f_r_b.output{3,1}=repmat(f_r_b.output{3,1},1,size(f_r_b.output{2,1},2));




f_s_nb.output{4,1}=f_s_nb.output{2,1}./f_s_nb.output{3,1};

f_s_b.output{4,1}=f_s_b.output{2,1}./f_s_b.output{3,1};

f_r_nb.output{4,1}=f_r_nb.output{2,1}./f_r_nb.output{3,1};

f_r_b.output{4,1}=f_r_b.output{2,1}./f_r_b.output{3,1};
% gaussian smoothing
sigma=1;
width = round((6*sigma - 1)/2);
support = (-width:width);
gaussFilter = exp( -(support).^2 ./ (2*sigma^2) );
gaussFilter = gaussFilter/ sum(gaussFilter);

win=gausswin(64,1);
win=win./sum(win);
%convc1vsc2{ii,jj}=conv(spike_xcorr.c1vsc2,win);
%convc1vsc2=[convc1vsc2;conv(spike_xcorr.c1vsc2,gaussFilter,'same')];

for i=1:size(f_s_nb.output{3,1},1)
    figure;
    hold on;
    
    plot(f_s_nb.output{1,1}(1,2).time(1,:),conv(f_s_nb.output{4,1}(i,:),gaussFilter,'same'),'-r'); % second epoch has time references
    
    
    
    
    
    plot(f_s_b.output{1,1}(1,2).time(1,:),conv(f_s_b.output{4,1}(i,:),gaussFilter,'same'),'--r');
    plot(f_r_nb.output{1,1}(1,2).time(1,:),conv(f_r_nb.output{4,1}(i,:),gaussFilter,'same'),'-b');
    plot(f_r_b.output{1,1}(1,2).time(1,:),conv(f_r_b.output{4,1}(i,:),gaussFilter,'same'),'--b');
    set(gca,'XTick',-0.5:0.1:0.5);
    title(sprintf('Cross correlogram for %s day %s \n tet %d cell %d vs tet %d vs %d', animals{1,1},num2str(day), ...
        commoncells(i,1),commoncells(i,2),commoncells(i,3),commoncells(i,4)));
    legend(sprintf('< %s -barrier',num2str(minVPF)),sprintf('< %s +barrier',num2str(minVPF)),...
        sprintf('> %s -barrier',num2str(minVPF)),sprintf('> %s +barrier',num2str(minVPF)));
    
    % change to that directory and saves the figure with file name
    % animal_day_epoch
    cd(strcat(f_r_b.animal{1,2},'Plot/'));
    figurename = sprintf('CCG_%s_%d_%d_%dv%d_%d',animals{1,1},day,...
        commoncells(i,1),commoncells(i,2),commoncells(i,3),commoncells(i,4));
    
    saveas(gcf, figurename, 'pdf');
    close;
    
    
    
    
end


figure;
boxplot([f_s_b.output{3,1}(:,1);f_s_nb.output{3,1}(:,1);f_r_b.output{3,1}(:,1);f_r_nb.output{3,1}(:,1)],...
    [ones(size(commoncells,1),1).*1;ones(size(commoncells,1),1).*2;ones(size(commoncells,1),1).*3;ones(size(commoncells,1),1).*4],...
    'labels',{'<2cm/s +barrier','<2cm/s -barrier','>2cm/s +barrier','>2cm/s -barrier'},'labelorientation','inline');
title(sprintf('Spike number for each group for %s day %s', animals{1,1},num2str(day)));


cd(strcat(f_r_b.animal{1,2},'Plot/'));
figurename = sprintf('spikesummary_%s_%d',animals{1,1},day);

saveas(gcf, figurename, 'pdf');
close;














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
