% plots cross correlograms for co-firing of cell parts under different
% conditions
% also shows place fields

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
minVPF = 0.5; %cm/sec
minPeakPF = 3;
lessthan=0;
includestates = 6;

%Animal selection
%-----------------------------------------------------
animals = {'I1'};
%-----------------------------------------------------

%Filter creation
%--------------------------------------------------------
% day filterionno


days = '[10]';

%% Get data for CGG
epochfilter = ['isequal($epochtype, ''Run'')'];
%epochfilter{1} = ['isequal($epoch, 4)'];

cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 10) && ($numspikes>100))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%timefilter = { {'JY_getlinvelocity', '($velocity <0.07) & ($velocity >-0.07)'} };

timefilter_r_b = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
    {'JY_getlinvelocity', strcat('abs($velocity > ',num2str(minVPF),')')},...
    {'JY_getbarrier','($barrier== 1)'}}; % barrier=0 means get all non-barrier events

timefilter_r_nb = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
    {'JY_getlinvelocity', strcat('abs($velocity > ',num2str(minVPF),')')},...
    {'JY_getbarrier','($barrier== 0)'}}; % barrier=0 means get all non-barrier events

timefilter_s_b = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
    {'JY_getlinvelocity', strcat('abs($velocity < ',num2str(minVPF),')')},...
    {'JY_getbarrier','($barrier== 1)'}}; % barrier=0 means get all non-barrier events

timefilter_s_nb = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
    {'JY_getlinvelocity', strcat('abs($velocity < ',num2str(minVPF),')')},...
    {'JY_getbarrier','($barrier== 0)'}}; % barrier=0 means get all non-barrier events


%------------------------------------------------------
f_r_b = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter_r_b);
f_r_nb = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter_r_nb);
f_s_b = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter_s_b);
f_s_nb = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter_s_nb);


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

f_r_b = setfilteriterator(f_r_b,iterator);
f_r_nb = setfilteriterator(f_r_nb,iterator);
f_s_b = setfilteriterator(f_s_b,iterator);
f_s_nb = setfilteriterator(f_s_nb,iterator);


f_r_b = setfilterfunction(f_r_b, 'JY_spkxcorr', {'spikes'}, 0.5, includestates, minV, f_r_b);
f_r_nb = setfilterfunction(f_r_nb, 'JY_spkxcorr', {'spikes'}, 0.5, includestates, minV, f_r_nb);
f_s_b = setfilterfunction(f_s_b, 'JY_spkxcorr', {'spikes'}, 0.5, includestates, minV, f_s_b);
f_s_nb = setfilterfunction(f_s_nb, 'JY_spkxcorr', {'spikes'}, 0.5, includestates, minV, f_s_nb);

f_r_b = runfilter(f_r_b);
f_r_nb = runfilter(f_r_nb);
f_s_b = runfilter(f_s_b);
f_s_nb = runfilter(f_s_nb);

%% Get data for place fields


epochtype='Run';
epochfilter = ['isequal($epochtype, ''Run'')'];
%epochfilter{1} = ['isequal($epoch, 2)'];
cellfilter = '(isequal($area, ''CA1'') && ($meanrate >0 ))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
timefilter = { {'JY_getriptimes','($nripples == 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},{'JY_getlinvelocity', strcat('$velocity > ',num2str(minVPF))}};
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);
iterator = 'singlecellanal';
f = setfilteriterator(f,iterator);
out=setfilterfunction(f, 'JY_calcopenfieldoccupancy', {'spikes','data'});
out=runfilter(out);
outdata=out.output{1,1};


%% load well information
% posfilename = strcat(datadir,animaldir,'/',animals{1,1},'data',dsz,dayt,'.mat');
% load(posfilename);
% 
% posfilename = strcat(datadir,animaldir,'/',animals{1,1},'linpos',dsz,dayt,'.mat');
% load(posfilename);
% 
% % get rewarded wells
% Wells=data{1,day}{1,epoch}.Wellinfo.rewardedwells;
% % get reward well positions
% if strmatch(epochtype,'Run');
%     Wellpos=linpos{1,day}{1,epoch}.wellSegmentInfo.wellCoord;
% end


%% put data together, placefields of the two cells for each epoch of the
% day along side eachother


f_r_b.output{2,1}={};
f_r_nb.output{2,1}={};
f_s_b.output{2,1}={};
f_s_nb.output{2,1}={};


day=f_s_nb.epochs{1,1}(1,1);

if ismember(day,[10 11 12 13])
    % comment out for non barrier days
    commoncells=mintersect(f_r_nb.output{1,1}(1,1).cellpairindex,...
        f_r_nb.output{1,1}(1,2).cellpairindex,...
        f_r_nb.output{1,1}(1,3).cellpairindex,...
        f_r_b.output{1,1}(1,2).cellpairindex,...
        f_r_b.output{1,1}(1,3).cellpairindex,...
        f_s_nb.output{1,1}(1,1).cellpairindex,...
        f_s_nb.output{1,1}(1,2).cellpairindex,...
        f_s_nb.output{1,1}(1,3).cellpairindex,...
        f_s_b.output{1,1}(1,2).cellpairindex,...
        f_s_b.output{1,1}(1,3).cellpairindex,...
        'rows');
else
    commoncells=mintersect(f_r_nb.output{1,1}(1,1).cellpairindex,...
        f_r_nb.output{1,1}(1,2).cellpairindex,...
        f_r_nb.output{1,1}(1,3).cellpairindex,...
        f_s_nb.output{1,1}(1,1).cellpairindex,...
        f_s_nb.output{1,1}(1,2).cellpairindex,...
        f_s_nb.output{1,1}(1,3).cellpairindex,...
        'rows');
    
end


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

% collect epoch data
daydata=[];
epochindex=[];
for i=1:size(out.data{1,1},2)
    daydata=[daydata;out.data{1,1}{1,i}];
    epochindex=[epochindex;repmat(out.epochs{1,1}(i,:),size(out.data{1,1}{1,i},1),1)];
end

%% plot placefield and CCG on same image

% generate colour map
nc = 1024;
cmap = jet(nc);
cmap(1,:) = 1;
colormap(cmap);
plfoverlap=[];
% for i=1:size(f_s_nb.output{3,1},1)
%     %figure;
%     
%     set(gcf,'position',[0 0 2000 800]);
%     set(gcf,'PaperPositionMode','auto');
%     
%     
%     cell1=[commoncells(i,1),commoncells(i,2)];
%     cell2=[commoncells(i,3),commoncells(i,4)];
%     
%     % find place fields matching the cells
%     
%     cell1index=find(daydata(:,1)==cell1(1,1) & daydata(:,2)==cell1(1,2));
%     cell1epochindex=epochindex(cell1index,:);
%     cell2index=find(daydata(:,1)==cell2(1,1) & daydata(:,2)==cell2(1,2));
%     cell2epochindex=epochindex(cell2index,:);
%     subwidth=6;
%     
%     % calculate placefield overlap
%     
%     
%     curroverlap=[];
%     
%     for k=1:size(cell1index,1)
%         
%         %%% plot placefield 1st cell
%         
%         subplot(3,subwidth,(k*subwidth-5))
%         imagedatas1=outdata{1,cell1index(k)}.smoothedspikerate;
%         imagesc(flipud(imagedatas1));
%         
%         %plot rewarded wells
%         if strmatch(epochtype,'Run');
%             [output.welloccupancy output.xticks output.yticks]=hist2(Wellpos(:,1),Wellpos(:,2), binx, biny);
%             output.welloccupancy = flipud(output.welloccupancy);
%             
%             Wellpos=linpos{1,cell1epochindex(k,1)}{1,cell1epochindex(k,2)}.wellSegmentInfo.wellCoord;
%             
%             
%             [wella, wellb]=find(output.welloccupancy==1);
%             Wellposs=[wella,wellb];
%             for wind=1:size(Wellpos,1);
%                 plot(Wellposs(wind,2),Wellposs(wind,1),'+w','MarkerSize',3);
%             end
%         end
%         
%         xlimval=[0 size(imagedatas1,2)];
%         ylimval=[0 size(imagedatas1,1)];
%         
%         distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
%         
%         text(distlabel(1),distlabel(2),sprintf('max %s Hz',num2str(max(imagedatas1(:)),2),'FontSize',3,'Color','k'));
%         
%         %%% plot placefield 2nd cell
%         
%         axis image;
%         set(gca,'xtick',[],'ytick',[]);
%         
%         subplot(3,subwidth,(k*subwidth-4))
%         imagedatas2=outdata{1,cell2index(k)}.smoothedspikerate;
%         imagesc(flipud(imagedatas2));
%         
%         if strmatch(epochtype,'Run');
%             [output.welloccupancy output.xticks output.yticks]=hist2(Wellpos(:,1),Wellpos(:,2), binx, biny);
%             output.welloccupancy = flipud(output.welloccupancy);
%             [wella, wellb]=find(output.welloccupancy==1);
%             Wellposs=[wella,wellb];
%             for wind=1:size(Wellpos,1);
%                 plot(Wellposs(wind,2),Wellposs(wind,1),'+w','MarkerSize',3);
%             end
%         end
%         
%         xlimval=[0 size(imagedatas2,2)];
%         ylimval=[0 size(imagedatas2,1)];
%         
%         distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
%         
%         text(distlabel(1),distlabel(2),sprintf('max %s Hz',num2str(max(imagedatas2(:)),2),'FontSize',3,'Color','k'));
%        
%         axis image;
%         set(gca,'xtick',[],'ytick',[]);
%         
%         imagedatas1=imagedatas1./max(imagedatas1(:));
%         
%         imagedatas2=imagedatas2./max(imagedatas2(:));
%         imagedatas1(imagedatas1<3)=0;
%         imagedatas1(imagedatas1>=3)=1;
%         imagedatas2(imagedatas2<3)=0;
%         imagedatas2(imagedatas2>=3)=1;
%         
%         overlap=corr2(imagedatas1,imagedatas2);
%         curroverlap=[curroverlap;overlap];
%         
%     end
%     
%     plfoverlap=[plfoverlap;mean(curroverlap)];
%     
%     %%% plot CCG
%     
%     subplot(3,subwidth,[3:subwidth (subwidth+3):2*subwidth (2*subwidth+3):3*subwidth])
%     
%     hold on;
% 
%     plot(f_s_nb.output{1,1}(1,2).time(1,:),conv(f_s_nb.output{4,1}(i,:),gaussFilter,'same'),'-r'); % second epoch has time references
%     %plot(f_s_b.output{1,1}(1,2).time(1,:),conv(f_s_b.output{4,1}(i,:),gaussFilter,'same'),'--r');
%     plot(f_r_nb.output{1,1}(1,2).time(1,:),conv(f_r_nb.output{4,1}(i,:),gaussFilter,'same'),'-b');
%     %plot(f_r_b.output{1,1}(1,2).time(1,:),conv(f_r_b.output{4,1}(i,:),gaussFilter,'same'),'--b');
%     set(gca,'XTick',-0.5:0.1:0.5);
%     title(sprintf('Cross correlogram for %s day %s \n tet %d cell %d vs tet %d cell %d Placefield corr %s', animals{1,1},num2str(day), ...
%         commoncells(i,1),commoncells(i,2),commoncells(i,3),commoncells(i,4), num2str(plfoverlap(i,1),4)));
%     
%     legend(sprintf('< %s -barrier %.0f spikes',num2str(minVPF),f_s_nb.output{3,1}(i,1)),sprintf('< %s +barrier %.0f spikes',num2str(minVPF),f_s_b.output{3,1}(i,1)),...
%         sprintf('> %s -barrier %.0f spikes',num2str(minVPF),f_r_nb.output{3,1}(i,1)),sprintf('> %s +barrier %.0f spikes',num2str(minVPF),f_r_b.output{3,1}(i,1)));
%     
%     
%     %   legend(sprintf('< %s -barrier %.0f spikes',num2str(minVPF),f_s_nb.output{3,1}(i,1)),...
%     %       sprintf('> %s -barrier %.0f spikes',num2str(minVPF),f_r_nb.output{3,1}(i,1)));
%     
%     % change to that directory and saves the figure with file name
%     % animal_day_epoch
%     cd(strcat(f_r_b.animal{1,2},'Plot/'));
%     figurename = sprintf('CCG_plf_%s_%d_%d_%dv%d_%d',animals{1,1},day,...
%         commoncells(i,1),commoncells(i,2),commoncells(i,3),commoncells(i,4));
%     
%     
%     
%     set(gcf,'PaperPositionMode','auto');
%     
%     papersize = get(gcf, 'PaperSize');
%     
%     set(gcf,'PaperSize',[25 15]);
%     
%     
%     saveas(gcf, figurename, 'pdf');
%     close;
% end


% figure;
% boxplot([f_s_b.output{3,1}(:,1);f_s_nb.output{3,1}(:,1);f_r_b.output{3,1}(:,1);f_r_nb.output{3,1}(:,1)],...
%     [ones(size(commoncells,1),1).*1;ones(size(commoncells,1),1).*2;ones(size(commoncells,1),1).*3;ones(size(commoncells,1),1).*4],...
%     'labels',{'<2cm/s +barrier','<2cm/s -barrier','>2cm/s +barrier','>2cm/s -barrier'},'labelorientation','inline');
% title(sprintf('Spike number for each group for %s day %s', animals{1,1},num2str(day)));
%
% cd(strcat(f_r_b.animal{1,2},'Plot/'));
% figurename = sprintf('spikesummary_%s_%d',animals{1,1},day);
%
% saveas(gcf, figurename, 'pdf');
% close;


%% plot max normalised CCG


f_s_nb_max=max(f_s_nb.output{4,1},[],2);
f_s_b_max=max(f_s_b.output{4,1},[],2);
f_r_nb_max=max(f_r_nb.output{4,1},[],2);
f_r_b_max=max(f_r_b.output{4,1},[],2);

% find if point is above or below line of equivalence
f_s_test=f_s_b_max>f_s_nb_max;
f_r_test=f_r_b_max>f_r_nb_max;

f_s_test_p=chisquarecont([sum(f_s_test==1) sum(f_s_test==0);size(f_s_test,1)/2 size(f_s_test,1)/2]);
f_r_test_p=chisquarecont([sum(f_r_test==1) sum(f_r_test==0);size(f_r_test,1)/2 size(f_r_test,1)/2]);

%figure;
hold on;

f_s=polyfit(f_s_nb_max,f_s_b_max,1);
f_r=polyfit(f_r_nb_max,f_r_b_max,1);

fs=polyval(f_s,f_s_nb_max);
fr=polyval(f_r,f_r_nb_max);

% draw equivalence line
l=line([0 1],[0 1],'LineStyle','--','Color','k');
set(get(get(l,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off');

plot(f_s_nb_max,f_s_b_max,'or'); %, f_s_nb_max,fs,'--r');
plot(f_r_nb_max,f_r_b_max,'ob'); %,f_r_nb_max,fr,'--b');

title(sprintf('Maximum normalised cross correlation for %s day %s \n <2cm/s p= %.4f >2cm/s p= %.4f',...
    animals{1,1},num2str(day),f_s_test_p,f_r_test_p));
xlabel('Without barrier');
ylabel('With barrier');
legend('<2cm/s ripples','>2cm/s ripples');

cd(strcat(f_r_b.animal{1,2},'Plot/'));
figurename = sprintf('MaxCCG_%s_%d',animals{1,1},day);

saveas(gcf, figurename, 'pdf');
close;


%% plot placefield similarity vc correlation
%
% figure; hold on;
% plot(plfoverlap,f_s_nb_max,'.r');
% plot(plfoverlap,f_s_b_max,'.b');
%
% title(sprintf('Maximum normalised cross correlation vs placefield correlation \n <2cm/s for %s day %s', animals{1,1},num2str(day)));
% xlabel('Placefield 2D correlation');
% ylabel('Normalised cross correlation');
% legend('Without barrier','With barrier');
%
% cd(strcat(f_r_b.animal{1,2},'Plot/'));
% figurename = sprintf('MaxCCG_plf_still_%s_%d',animals{1,1},day);
%
% saveas(gcf, figurename, 'pdf');
% close;

%% Maximum normalised cross correlation vs placefield correlation
%
% figure; hold on;
% plot(plfoverlap,f_r_nb_max,'.r');
% plot(plfoverlap,f_r_b_max,'.b');
%
% title(sprintf('Maximum normalised cross correlation vs placefield correlation \n >2cm/s for %s day %s', animals{1,1},num2str(day)));
% xlabel('Placefield 2D correlation');
% ylabel('Normalised cross correlation');
% legend('Without barrier','With barrier');
%
% cd(strcat(f_r_b.animal{1,2},'Plot/'));
% figurename = sprintf('MaxCCG_plf_move_%s_%d',animals{1,1},day);
%
% saveas(gcf, figurename, 'pdf');
% close;
%
%
%
% figure; hold on;
% plot(plfoverlap,f_r_nb_max,'.r');
% plot(plfoverlap,f_s_nb_max,'.b');
%
% title(sprintf('Maximum normalised cross correlation vs placefield correlation \n for %s day %s', animals{1,1},num2str(day)));
% xlabel('Placefield 2D correlation');
% ylabel('Normalised cross correlation');
% legend('>2cm/s','<2cm/s');
%
% cd(strcat(f_r_b.animal{1,2},'Plot/'));
% figurename = sprintf('MaxCCG_plf_nobarrier_%s_%d',animals{1,1},day);
%
% saveas(gcf, figurename, 'pdf');
% close;
%

%% show distribution of max CCG
%
% Rank sum test

% [p,h,stats]=ranksum(f_r_nb_max,f_s_nb_max);
% [pk,hk]=kstest2(f_r_nb_max,f_s_nb_max);
% 
% [ps,hs,stats]=ranksum(f_r_nb_max,f_r_b_max);
% 
% figure;
% boxplot([f_r_nb_max;f_r_b_max],[zeros(size(f_r_b_max,1)); ones(size(f_s_b_max,1))],'labels',{sprintf('> %d cm/s',minVPF),sprintf('< %d cm/s',minVPF)});
% title(sprintf('Maximum normalised cross correlation for all cell pairs \n for %s day %s p= %.2d', animals{1,1},num2str(day), ps));
% 
% 
% figure;
% boxplot([f_r_nb_max;f_s_nb_max],[zeros(size(f_r_nb_max,1)); ones(size(f_s_nb_max,1))],'labels',{sprintf('> %d cm/s',minVPF),sprintf('< %d cm/s',minVPF)});
% title(sprintf('Maximum normalised cross correlation for all cell pairs \n for %s day %s p= %.2d', animals{1,1},num2str(day), p));
% 
% ylabel('Normalised cross correlation');
% cd(strcat(f_r_b.animal{1,2},'Plot/'));
% figurename = sprintf('MaxCCG_nobarrier_%d_%s_%d',minVPF,animals{1,1},day);
% saveas(gcf, figurename, 'pdf');
% close;



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
