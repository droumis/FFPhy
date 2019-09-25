Veqn = '>=0';
minV =  str2num(Veqn(end));
maxstage = 3; % [1 2 3]
minVPF = 2; %cm/sec
minPeakPF = 3;
lessthan=0;
includestates = 6;

%Animal selection
%-----------------------------------------------------
animals = {'I1'};
%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filterionno

%days='[1:10]';%,'1:10';
%days = '[1:1]';
days = '[11]';

%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

%epochfilter{1} = ['isequal($epochtype, ''Run'')'];

epochfilter{1} = ['isequal($epoch,  4)'];


cellfilter = '(isequal($area, ''CA1'') && ($meanrate < 200))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%timefilter = { {'JY_getlinvelocity', '($velocity <0.07) & ($velocity >-0.07)'} };

%timefilter = { {'JY_getriptimes','($nripples == 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter = { {'JY_getriptimes','($nripples == 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
%    {'JY_getlinvelocity', strcat('$velocity > ',num2str(minVPF))},{'JY_getbarrier','($barrier== 0)'}};
timefilter = { {'JY_getriptimes','($nripples == 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
    {'JY_getlinvelocity', strcat('$velocity > ',num2str(minVPF))}};

%timefilter = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}, {'JY_getlinvelocity', '$velocity <0.01'}};
%timefilter2 = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter2 = { {'getriptimes','($nripples > 0)',[],'minthresh',2,'cellfilter',cellfilter,'tetfilter',[1 2 3 4]}};

%timefilter = { {'JY_getlinvelocity', '(($velocity) >= 0))', 6} };
f = DFL_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);

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

f = setfilterfunction(f, 'JY_DFL_linplacefield', {'spikes', 'linpos','data'}, minV,f);
% out = plottrajdata(index, excludetimes, spikes, linpos, includestates, minV, varargin)

f = runfilter(f);

% find cells with placefields on segment with barrier

segmentwithbarrier=f.output{1,1}{1,1}.barriersegment;

cellindex=cellfun(@(x) x.segplace(1,segmentwithbarrier),f.output{1,1});
plfonseg=f.data{1,1}{1,1}(cellindex==1,:);

%% calculate data for differen conditions
epochfilter{1} = ['isequal($epochtype, ''Run'')'];
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

% calculate area under the normalised curve

% x values
xvals=f_s_nb.output{1,1}(1,1).time(1,:);
%f_s_nb.output{5,1}=trapz(repmat(xvals,size(f_s_nb.output{4,1},1),1),f_s_nb.output{4,1},1);

f_s_nb.output{5,1}=trapz(xvals,f_s_nb.output{4,1},2);
f_s_b.output{5,1}=trapz(xvals,f_s_b.output{4,1},2);
f_r_nb.output{5,1}=trapz(xvals,f_r_nb.output{4,1},2);
f_r_b.output{5,1}=trapz(xvals,f_r_b.output{4,1},2);

%% classify cell pair comparisons
% either as one cell is on barrier segment, both are or none are

cell1=ismember(commoncells(:,1:2),plfonseg,'rows');
cell2=ismember(commoncells(:,3:4),plfonseg,'rows');
comparisontype=cell1+cell2;








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

% figure;
% hold on;
%
% f_s=polyfit(f_s_nb_max,f_s_b_max,1);
% f_r=polyfit(f_r_nb_max,f_r_b_max,1);
%
% fs=polyval(f_s,f_s_nb_max);
% fr=polyval(f_r,f_r_nb_max);
%
% % draw equivalence line
% l=line([0 1],[0 1],'LineStyle','--','Color','k');
% set(get(get(l,'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off');
%
% plot(f_s_nb_max,f_s_b_max,'or'); %, f_s_nb_max,fs,'--r');
% plot(f_r_nb_max,f_r_b_max,'ob'); %,f_r_nb_max,fr,'--b');
%
% title(sprintf('Maximum normalised cross correlation for %s day %s \n <2cm/s p= %.4f >2cm/s p= %.4f',...
%     animals{1,1},num2str(day),f_s_test_p,f_r_test_p));
% xlabel('Without barrier');
% ylabel('With barrier');
% legend('<2cm/s ripples','>2cm/s ripples');
%
% cd(strcat(f_r_b.animal{1,2},'Plot/'));
% figurename = sprintf('MaxCCG_%s_%d',animals{1,1},day);

%saveas(gcf, figurename, 'pdf');
%close;

%% plot max normalised CCG for different cell pair comparisons

% figure;
%      set(gcf,'position',[0 0 2000 500]);
%      set(gcf,'PaperPositionMode','auto');
% subplot(1,3,1);
%  hold on;
%
%  f_s_test=f_s_b_max(comparisontype==0)>f_s_nb_max(comparisontype==0);
% f_r_test=f_r_b_max(comparisontype==0)>f_r_nb_max(comparisontype==0);
%
% f_s_test_p=chisquarecont([sum(f_s_test==1) sum(f_s_test==0);size(f_s_test,1)/2 size(f_s_test,1)/2]);
% f_r_test_p=chisquarecont([sum(f_r_test==1) sum(f_r_test==0);size(f_r_test,1)/2 size(f_r_test,1)/2]);
%
%
% % draw equivalence line
% l=line([0 1],[0 1],'LineStyle','--','Color','k');
% set(get(get(l,'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off');
%
% plot(f_s_nb_max(comparisontype==0),f_s_b_max(comparisontype==0),'or'); %, f_s_nb_max,fs,'--r');
% plot(f_r_nb_max(comparisontype==0),f_r_b_max(comparisontype==0),'ob'); %,f_r_nb_max,fr,'--b');
%
% title(sprintf('None barrier vs none barrier for %s day %s \n <2cm/s p= %.4f >2cm/s p= %.4f',...
%     animals{1,1},num2str(day),f_s_test_p,f_r_test_p));
% xlabel('Without barrier');
% ylabel('With barrier');
%
%
% subplot(1,3,2);
%
%
% f_s_test=f_s_b_max(comparisontype==1)>f_s_nb_max(comparisontype==1);
% f_r_test=f_r_b_max(comparisontype==1)>f_r_nb_max(comparisontype==1);
%
% f_s_test_p=chisquarecont([sum(f_s_test==1) sum(f_s_test==0);size(f_s_test,1)/2 size(f_s_test,1)/2]);
% f_r_test_p=chisquarecont([sum(f_r_test==1) sum(f_r_test==0);size(f_r_test,1)/2 size(f_r_test,1)/2]);
%
%
%
% hold on;
% % draw equivalence line
% l=line([0 1],[0 1],'LineStyle','--','Color','k');
% set(get(get(l,'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off');
%
% plot(f_s_nb_max(comparisontype==1),f_s_b_max(comparisontype==1),'or'); %, f_s_nb_max,fs,'--r');
% plot(f_r_nb_max(comparisontype==1),f_r_b_max(comparisontype==1),'ob'); %,f_r_nb_max,fr,'--b');
%
% title(sprintf('Barrier vs non-barrier for %s day %s \n <2cm/s p= %.4f >2cm/s p= %.4f',...
%     animals{1,1},num2str(day),f_s_test_p,f_r_test_p));
% xlabel('Without barrier');
% ylabel('With barrier');
%
% subplot(1,3,3);
% hold on;
%
% f_s_test=f_s_b_max(comparisontype==2)>f_s_nb_max(comparisontype==2);
% f_r_test=f_r_b_max(comparisontype==2)>f_r_nb_max(comparisontype==2);
%
% f_s_test_p=chisquarecont([sum(f_s_test==1) sum(f_s_test==0);size(f_s_test,1)/2 size(f_s_test,1)/2]);
% f_r_test_p=chisquarecont([sum(f_r_test==1) sum(f_r_test==0);size(f_r_test,1)/2 size(f_r_test,1)/2]);
%
% % draw equivalence line
% l=line([0 1],[0 1],'LineStyle','--','Color','k');
% set(get(get(l,'Annotation'),'LegendInformation'),...
%     'IconDisplayStyle','off');
%
% plot(f_s_nb_max(comparisontype==2),f_s_b_max(comparisontype==2),'or'); %, f_s_nb_max,fs,'--r');
% plot(f_r_nb_max(comparisontype==2),f_r_b_max(comparisontype==2),'ob'); %,f_r_nb_max,fr,'--b');
%
% title(sprintf('Barrier vs Barrier for %s day %s \n <2cm/s p= %.4f >2cm/s p= %.4f',...
%     animals{1,1},num2str(day),f_s_test_p,f_r_test_p));
% xlabel('Without barrier');
% ylabel('With barrier');
%
%  set(gcf,'PaperSize',[20 8]);
%
% cd(strcat(f_r_b.animal{1,2},'Plot/'));
%  figurename = sprintf('MaxCCG_comparecondition_%s_%d',animals{1,1},day);
%  saveas(gcf, figurename, 'pdf');
%
%  close;


%% plot max normalised CCG for different cell pair comparisons

figure;
set(gcf,'position',[0 0 2000 500]);
set(gcf,'PaperPositionMode','auto');
subplot(1,2,1);

diff0=(f_s_nb_max(comparisontype==0)-f_s_b_max(comparisontype==0))./(f_s_nb_max(comparisontype==0)+f_s_b_max(comparisontype==0));
diff1=(f_s_nb_max(comparisontype==1)-f_s_b_max(comparisontype==1))./(f_s_nb_max(comparisontype==1)+f_s_b_max(comparisontype==1));
diff2=(f_s_nb_max(comparisontype==2)-f_s_b_max(comparisontype==2))./(f_s_nb_max(comparisontype==2)+f_s_b_max(comparisontype==2));

diff0_s=diff0;
diff1_s=diff1;
diff2_s=diff2;

if isempty(diff0)
    diff0=0;
end
if isempty(diff1)
    diff1=0;
end
if isempty(diff2)
    diff2=0;
end

% boxplot([diff0;diff1;diff2],[zeros(size(diff0));ones(size(diff1));ones(size(diff2))*2],'labels',...
%     {sprintf('nb vs nb %d',length(f_s_nb_max(comparisontype==0))),sprintf('nb vs b %d',length(f_s_nb_max(comparisontype==1))),...
%     sprintf('b vs b %d',length(f_s_nb_max(comparisontype==2)))},'labelorientation','inline');
% 
% title(sprintf('Change in normalised CCG <2cm/s for %s day %s', animals{1,1},num2str(day)));
% subplot(1,2,2);

diff0=(f_r_nb_max(comparisontype==0)-f_r_b_max(comparisontype==0))./(f_r_nb_max(comparisontype==0)+f_r_b_max(comparisontype==0));
diff1=(f_r_nb_max(comparisontype==1)-f_r_b_max(comparisontype==1))./(f_r_nb_max(comparisontype==1)+f_r_b_max(comparisontype==1));
diff2=(f_r_nb_max(comparisontype==2)-f_r_b_max(comparisontype==2))./(f_r_nb_max(comparisontype==2)+f_r_b_max(comparisontype==2));

diff0_r=diff0;
diff1_r=diff1;
diff2_r=diff2;


if isempty(diff0)
    diff0=0;
end
if isempty(diff1)
    diff1=0;
end
if isempty(diff2)
    diff2=0;
end

% boxplot([diff0;diff1;diff2],[zeros(size(diff0));ones(size(diff1));ones(size(diff2))*2],'labels',...
%     {sprintf('nb vs nb %d',length(f_r_nb_max(comparisontype==0))),sprintf('nb vs b %d',length(f_r_nb_max(comparisontype==1))),...
%     sprintf('b vs b %d',length(f_r_nb_max(comparisontype==2)))},'labelorientation','inline');
% 
% title(sprintf('Change in normalised CCG >2cm/s for %s day %s', animals{1,1},num2str(day)));
% cd(strcat(f_r_b.animal{1,2},'Plot/'));
% figurename = sprintf('MaxCCG_change_%s_%d',animals{1,1},day);
% 
% set(gcf,'PaperSize',[20 8]);
% 
% saveas(gcf, figurename, 'pdf');
% close;


%% plot differences in area under the curve of CCG for different cell pair comparisons
% 
% figure;
% set(gcf,'position',[0 0 2000 500]);
% set(gcf,'PaperPositionMode','auto');
% subplot(1,2,1);
% 
% diff0=(f_s_nb.output{5,1}(comparisontype==0)-f_s_b.output{5,1}(comparisontype==0))./(f_s_nb.output{5,1}(comparisontype==0)+f_s_b.output{5,1}(comparisontype==0));
% diff1=(f_s_nb.output{5,1}(comparisontype==1)-f_s_b.output{5,1}(comparisontype==1))./(f_s_nb.output{5,1}(comparisontype==1)+f_s_b.output{5,1}(comparisontype==1));
% diff2=(f_s_nb.output{5,1}(comparisontype==2)-f_s_b.output{5,1}(comparisontype==2))./(f_s_nb.output{5,1}(comparisontype==2)+f_s_b.output{5,1}(comparisontype==2));
% 
% if isempty(diff0)
%     diff0=0;
% end
% if isempty(diff1)
%     diff1=0;
% end
% if isempty(diff2)
%     diff2=0;
% end
% 
% 
% 
% 
% boxplot([diff0;diff1;diff2],[zeros(size(diff0));ones(size(diff1));ones(size(diff2))*2],'labels',...
%     {sprintf('nb vs nb %d',length(f_s_nb.output{5,1}(comparisontype==0))),sprintf('nb vs b %d',length(f_s_nb.output{5,1}(comparisontype==1))),...
%     sprintf('b vs b %d',length(f_s_nb.output{5,1}(comparisontype==2)))},'labelorientation','inline');
% 
% title(sprintf('Change in normalised CCG during ripples barrier vs no barrier \n <2cm/s for %s day %s', animals{1,1},num2str(day)));
% subplot(1,2,2);
% 
% diff0=(f_r_nb.output{5,1}(comparisontype==0)-f_r_b.output{5,1}(comparisontype==0))./(f_r_nb.output{5,1}(comparisontype==0)+f_r_b.output{5,1}(comparisontype==0));
% diff1=(f_r_nb.output{5,1}(comparisontype==1)-f_r_b.output{5,1}(comparisontype==1))./(f_r_nb.output{5,1}(comparisontype==1)+f_r_b.output{5,1}(comparisontype==1));
% diff2=(f_r_nb.output{5,1}(comparisontype==2)-f_r_b.output{5,1}(comparisontype==2))./(f_r_nb.output{5,1}(comparisontype==2)+f_r_b.output{5,1}(comparisontype==2));
% 
% if isempty(diff0)
%     diff0=0;
% end
% if isempty(diff1)
%     diff1=0;
% end
% if isempty(diff2)
%     diff2=0;
% end
% boxplot([diff0;diff1;diff2],[zeros(size(diff0));ones(size(diff1));ones(size(diff2))*2],'labels',...
%     {sprintf('nb vs nb %d',length(f_r_nb.output{5,1}(comparisontype==0))),sprintf('nb vs b %d',length(f_r_nb.output{5,1}(comparisontype==1))),...
%     sprintf('b vs b %d',length(f_r_nb.output{5,1}(comparisontype==2)))},'labelorientation','inline');
% 
% title(sprintf('Change in normalised CCG during ripples barrier vs no barrier \n >2cm/s for %s day %s', animals{1,1},num2str(day)));
% cd(strcat(f_r_b.animal{1,2},'Plot/'));
% figurename = sprintf('AreaCCG_change_%s_%d',animals{1,1},day);
% 
% 
% set(gcf,'PaperSize',[20 8]);
% 
% saveas(gcf, figurename, 'pdf');
% close;

%% sort cells according to change in max CCG

% sort indices

[diff0_s_s diff0_s_i]=sort(diff0_s);
[diff1_s_s diff1_s_i]=sort(diff1_s);
[diff2_s_s diff2_s_i]=sort(diff2_s);

[diff0_r_s diff0_r_i]=sort(diff0_r);
[diff1_r_s diff1_r_i]=sort(diff1_r);
[diff2_r_s diff2_r_i]=sort(diff2_r);

% find cell pairs

goodcellpairs0s=commoncells(comparisontype==0,:);
goodcellpairs1s=commoncells(comparisontype==1,:);
goodcellpairs2s=commoncells(comparisontype==2,:);

goodcellpairs0s=goodcellpairs0s(diff0_s_i,:);
goodcellpairs0s=[goodcellpairs0s diff0_s_s];

goodcellpairs1s=goodcellpairs1s(diff1_s_i,:);
goodcellpairs1s=[goodcellpairs1s diff1_s_s];

goodcellpairs2s=goodcellpairs2s(diff2_s_i,:);
goodcellpairs2s=[goodcellpairs2s diff2_s_s];


goodcellpairs0r=commoncells(comparisontype==0,:);
goodcellpairs1r=commoncells(comparisontype==1,:);
goodcellpairs2r=commoncells(comparisontype==2,:);

% show differences in max CCG and also area under the CCG

goodcellpairs0r=goodcellpairs0r(diff0_r_i,:);
goodcellpairs0r=[goodcellpairs0r diff0_r_s];

goodcellpairs1r=goodcellpairs1r(diff1_r_i,:);
goodcellpairs1r=[goodcellpairs1r diff1_r_s];

goodcellpairs2r=goodcellpairs2r(diff2_r_i,:);
goodcellpairs2r=[goodcellpairs2r diff2_r_s];


% plot placefields of cell pairs



%% get placefield data
epochtype='Run';
epochfilter = ['isequal($epochtype, ''Run'')'];
%epochfilter{1} = ['isequal($epoch, 2)'];
cellfilter = '(isequal($area, ''CA1'') && ($meanrate <10 ))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
timefilter = { {'JY_getriptimes','($nripples == 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},{'JY_getlinvelocity', strcat('$velocity > ',num2str(minVPF))}};
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);
iterator = 'singlecellanal';
f = setfilteriterator(f,iterator);
out=setfilterfunction(f, 'JY_calcopenfieldoccupancy', {'spikes','data'});
out=runfilter(out);
outdata=out.output{1,1};




% collect epoch data
daydata=[];
epochindex=[];
for i=1:size(out.data{1,1},2)
    daydata=[daydata;out.data{1,1}{1,i}];
    epochindex=[epochindex;repmat(out.epochs{1,1}(i,:),size(out.data{1,1}{1,i},1),1)];
end

% plot combined place fields for all cells

for i=1:size(epochindex,1)
    totalimage(:,:,i)=outdata{1,i}.smoothedspikerate;
    end
    
    imagedata=mean(totalimage,3);
    imagesc(flipud(imagedata));
    title(sprintf('Combined place field for all cells \n %s for day %s ',...
        animals{1,1}, num2str(day)));
    cd(strcat(f_r_b.animal{1,2},'Plot/'));
figurename = sprintf('Combined_plf_%s_%s',animals{1,1},days);


saveas(gcf, figurename, 'pdf');
close;



% plot cell pairs both without placefields on barrier segment
figure;
 set(gcf,'position',[0 0 500 2000]);
 set(gcf,'PaperPositionMode','auto');
% number of rows
nrow0r=sum(isnan(goodcellpairs0r(:,5)));
n=1
  while n<=nrow0r
    
    % lookup indices of placefield for each cell
    cell1=goodcellpairs0r(n,1:2);
    cell2=goodcellpairs0r(n,3:4);
    change=goodcellpairs1r(n,5);
    % find place fields matching the cells
    cell1index=find(daydata(:,1)==cell1(1,1) & daydata(:,2)==cell1(1,2));
    cell1epochindex=epochindex(cell1index,:);
    cell2index=find(daydata(:,1)==cell2(1,1) & daydata(:,2)==cell2(1,2));
    cell2epochindex=epochindex(cell2index,:);
    subwidth=6;
    subplot(nrow0r/2,2,n-1+1)
     imagedatas1=[];
     
    for k=1:size(cell1index,1)
        imagedatas1(:,:,k)=outdata{1,cell1index(k)}.smoothedspikerate;
    end
    
    imagedata=mean(imagedatas1,3);
    imagesc(flipud(imagedata));
    
    % Keeps the proportions of the data
            axis image;
            % no ticks
            set(gca,'xtick',[],'ytick',[]);
      xlimval=[0 size(imagedata,2)];
            ylimval=[0 size(imagedata,1)];
            
            distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
            
            text(distlabel(1),distlabel(2),sprintf('t %s c %s %s',num2str(cell1(1)),num2str(cell1(2)),num2str(change)),'FontSize',10,'Color','w');
            
    subplot(nrow0r/2,2,n-1+2)
    imagedatas2=[];
    
    for k=1:size(cell2index,1)
        imagedatas2(:,:,k)=outdata{1,cell2index(k)}.smoothedspikerate;
    end
    
    imagedata=mean(imagedatas2,3);
    imagesc(flipud(imagedata));
    
    % Keeps the proportions of the data
            axis image;
            % no ticks
            set(gca,'xtick',[],'ytick',[]);
      xlimval=[0 size(imagedata,2)];
            ylimval=[0 size(imagedata,1)];
            
            distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
            
            text(distlabel(1),distlabel(2),sprintf('t %s c %s',num2str(cell2(1)),num2str(cell2(2)),num2str(change)),'FontSize',8,'Color','w');
          
            
            n=n+2
    end
    ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
        1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 0.99,sprintf('Placefield of cell pairs with no cell on barrier segment ripples >2cm/s \n %s for day %s ',...
        animals{1,1}, num2str(day)),'HorizontalAlignment','center','VerticalAlignment', 'top'); 

    cd(strcat(f_r_b.animal{1,2},'Plot/'));
figurename = sprintf('Cellpair_plf_0r_%s_%d',animals{1,1},day);

set(gcf,'PaperSize',[8 20]);

saveas(gcf, figurename, 'pdf');
close;
% plot cell pairs with one placefield on barrier segment
figure;

 set(gcf,'position',[0 0 500 2000]);
 set(gcf,'PaperPositionMode','auto');

 
 % number of rows
nrow1r=sum(isnan(goodcellpairs1r(:,5)));
n=1
  while n<=nrow1r
    
    % lookup indices of placefield for each cell
    cell1=goodcellpairs1r(n,1:2);
    cell2=goodcellpairs1r(n,3:4);
    change=goodcellpairs1r(n,5);
    % find place fields matching the cells
    cell1index=find(daydata(:,1)==cell1(1,1) & daydata(:,2)==cell1(1,2));
    cell1epochindex=epochindex(cell1index,:);
    cell2index=find(daydata(:,1)==cell2(1,1) & daydata(:,2)==cell2(1,2));
    cell2epochindex=epochindex(cell2index,:);
    subwidth=6;
    subplot(nrow1r/2,2,n-1+1)
     imagedatas1=[];
     
    for k=1:size(cell1index,1)
        imagedatas1(:,:,k)=outdata{1,cell1index(k)}.smoothedspikerate;
    end
    
    imagedata=mean(imagedatas1,3);
    imagesc(flipud(imagedata));
    
    % Keeps the proportions of the data
            axis image;
            % no ticks
            set(gca,'xtick',[],'ytick',[]);
      xlimval=[0 size(imagedata,2)];
            ylimval=[0 size(imagedata,1)];
            
            distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
            
            text(distlabel(1),distlabel(2),sprintf('t %s c %s %s',num2str(cell1(1)),num2str(cell1(2)),num2str(change)),'FontSize',8,'Color','w');
            
    subplot(nrow1r/2,2,n-1+2)
    imagedatas2=[];
    
    for k=1:size(cell2index,1)
        imagedatas2(:,:,k)=outdata{1,cell2index(k)}.smoothedspikerate;
    end
    
    imagedata=mean(imagedatas2,3);
    imagesc(flipud(imagedata));
    
    % Keeps the proportions of the data
            axis image;
            % no ticks
            set(gca,'xtick',[],'ytick',[]);
      xlimval=[0 size(imagedata,2)];
            ylimval=[0 size(imagedata,1)];
            
            distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
            
            text(distlabel(1),distlabel(2),sprintf('t %s c %s %s',num2str(cell2(1)),num2str(cell2(2)),num2str(change)),'FontSize',10,'Color','w');
          
            
            n=n+2
    end
    
  ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
        1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 0.99,sprintf('Placefield of cell pairs with 1 cell on barrier segment ripples >2cm/s \n %s for day %s ',...
        animals{1,1}, num2str(day)),'HorizontalAlignment','center','VerticalAlignment', 'top');

cd(strcat(f_r_b.animal{1,2},'Plot/'));
figurename = sprintf('Cellpair_plf_1r_%s_%d',animals{1,1},day);

set(gcf,'PaperSize',[8 20]);

saveas(gcf, figurename, 'pdf');
close;

% plot cell pairs with both placefield on barrier segment
figure;

 set(gcf,'position',[0 0 500 2000]);
 set(gcf,'PaperPositionMode','auto');
% number of rows
nrow2r=sum(isnan(goodcellpairs2r(:,5)));
n=1
  while n<=nrow2r
    
    % lookup indices of placefield for each cell
    cell1=goodcellpairs2r(n,1:2);
    cell2=goodcellpairs2r(n,3:4);
    change=goodcellpairs2r(n,5);
    % find place fields matching the cells
    cell1index=find(daydata(:,1)==cell1(1,1) & daydata(:,2)==cell1(1,2));
    cell1epochindex=epochindex(cell1index,:);
    cell2index=find(daydata(:,1)==cell2(1,1) & daydata(:,2)==cell2(1,2));
    cell2epochindex=epochindex(cell2index,:);
    subwidth=6;
    subplot(nrow2r/2,2,n-1+1)
     imagedatas1=[];
     
    for k=1:size(cell1index,1)
        imagedatas1(:,:,k)=outdata{1,cell1index(k)}.smoothedspikerate;
    end
    
    imagedata=mean(imagedatas1,3);
    imagesc(flipud(imagedata));
    
    % Keeps the proportions of the data
            axis image;
            % no ticks
            set(gca,'xtick',[],'ytick',[]);
      xlimval=[0 size(imagedata,2)];
            ylimval=[0 size(imagedata,1)];
            
            distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
            
            text(distlabel(1),distlabel(2),sprintf('t %s c %s %s',num2str(cell1(1)),num2str(cell1(2)),num2str(change)),'FontSize',8,'Color','w');
            
    subplot(nrow2r/2,2,n-1+2)
    imagedatas2=[];
    
    for k=1:size(cell2index,1)
        imagedatas2(:,:,k)=outdata{1,cell2index(k)}.smoothedspikerate;
    end
    
    imagedata=mean(imagedatas2,3);
    imagesc(flipud(imagedata));
    
    % Keeps the proportions of the data
            axis image;
            % no ticks
            set(gca,'xtick',[],'ytick',[]);
      xlimval=[0 size(imagedata,2)];
            ylimval=[0 size(imagedata,1)];
            
            distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
            
            text(distlabel(1),distlabel(2),sprintf('t %s c %s %s',num2str(cell2(1)),num2str(cell2(2)),num2str(change)),'FontSize',10,'Color','w');
          
            
            n=n+2
    end
    
  ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
        1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 0.99,sprintf('Placefield of cell pairs with both cells on barrier segment ripples >2cm/s \n %s for day %s ',...
        animals{1,1}, num2str(day)),'HorizontalAlignment','center','VerticalAlignment', 'top');
cd(strcat(f_r_b.animal{1,2},'Plot/'));
figurename = sprintf('Cellpair_plf_2r_%s_%d',animals{1,1},day);

set(gcf,'PaperSize',[8 20]);

saveas(gcf, figurename, 'pdf');
close;
% plot cell pairs with both placefield on barrier segment <2cm/s
figure;
 set(gcf,'position',[0 0 500 2000]);
 set(gcf,'PaperPositionMode','auto');

% number of rows
nrow0s=sum(isnan(goodcellpairs0s(:,5)));
n=1
  while n<=nrow0s
    
    % lookup indices of placefield for each cell
    cell1=goodcellpairs0s(n,1:2);
    cell2=goodcellpairs0s(n,3:4);
    change=goodcellpairs0s(n,5);
    % find place fields matching the cells
    cell1index=find(daydata(:,1)==cell1(1,1) & daydata(:,2)==cell1(1,2));
    cell1epochindex=epochindex(cell1index,:);
    cell2index=find(daydata(:,1)==cell2(1,1) & daydata(:,2)==cell2(1,2));
    cell2epochindex=epochindex(cell2index,:);
    subwidth=6;
    subplot(nrow0s/2,2,n-1+1)
     imagedatas1=[];
     
    for k=1:size(cell1index,1)
        imagedatas1(:,:,k)=outdata{1,cell1index(k)}.smoothedspikerate;
    end
    
    imagedata=mean(imagedatas1,3);
    imagesc(flipud(imagedata));
    
    % Keeps the proportions of the data
            axis image;
            % no ticks
            set(gca,'xtick',[],'ytick',[]);
      xlimval=[0 size(imagedata,2)];
            ylimval=[0 size(imagedata,1)];
            
            distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
            
            text(distlabel(1),distlabel(2),sprintf('t %s c %s %s',num2str(cell1(1)),num2str(cell1(2)),num2str(change)),'FontSize',8,'Color','w');
            
    subplot(nrow0s/2,2,n-1+2)
    imagedatas2=[];
    
    for k=1:size(cell2index,1)
        imagedatas2(:,:,k)=outdata{1,cell2index(k)}.smoothedspikerate;
    end
    
    imagedata=mean(imagedatas2,3);
    imagesc(flipud(imagedata));
    
    % Keeps the proportions of the data
            axis image;
            % no ticks
            set(gca,'xtick',[],'ytick',[]);
      xlimval=[0 size(imagedata,2)];
            ylimval=[0 size(imagedata,1)];
            
            distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
            
            text(distlabel(1),distlabel(2),sprintf('t %s c %s %s',num2str(cell2(1)),num2str(cell2(2)),num2str(change)),'FontSize',10,'Color','w');
          
            
            n=n+2
    end
    
  ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
        1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 0.99,sprintf('Placefield of cell pairs with no cells on barrier segment ripples <2cm/s \n %s for day %s ',...
        animals{1,1}, num2str(day)),'HorizontalAlignment','center','VerticalAlignment', 'top');

    cd(strcat(f_r_b.animal{1,2},'Plot/'));
figurename = sprintf('Cellpair_plf_0s_%s_%d',animals{1,1},day);

set(gcf,'PaperSize',[8 20]);

saveas(gcf, figurename, 'pdf');
close;
    
  % plot cell pairs with one placefield on barrier segment <2cm/s
figure;

 set(gcf,'position',[0 0 500 2000]);
 set(gcf,'PaperPositionMode','auto');
 
% number of rows
nrow1s=sum(isnan(goodcellpairs1s(:,5)));
n=1
  while n<=nrow1s
    
    % lookup indices of placefield for each cell
    cell1=goodcellpairs1s(n,1:2);
    cell2=goodcellpairs1s(n,3:4);
    change=goodcellpairs1s(n,5);
    % find place fields matching the cells
    cell1index=find(daydata(:,1)==cell1(1,1) & daydata(:,2)==cell1(1,2));
    cell1epochindex=epochindex(cell1index,:);
    cell2index=find(daydata(:,1)==cell2(1,1) & daydata(:,2)==cell2(1,2));
    cell2epochindex=epochindex(cell2index,:);
    subwidth=6;
    subplot(nrow1s/2,2,n-1+1)
     imagedatas1=[];
     
    for k=1:size(cell1index,1)
        imagedatas1(:,:,k)=outdata{1,cell1index(k)}.smoothedspikerate;
    end
    
    imagedata=mean(imagedatas1,3);
    imagesc(flipud(imagedata));
    
    % Keeps the proportions of the data
            axis image;
            % no ticks
            set(gca,'xtick',[],'ytick',[]);
      xlimval=[0 size(imagedata,2)];
            ylimval=[0 size(imagedata,1)];
            
            distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
            
            text(distlabel(1),distlabel(2),sprintf('t %s c %s %s',num2str(cell1(1)),num2str(cell1(2)),num2str(change)),'FontSize',8,'Color','w');
            
    subplot(nrow1s/2,2,n-1+2)
    imagedatas2=[];
    
    for k=1:size(cell2index,1)
        imagedatas2(:,:,k)=outdata{1,cell2index(k)}.smoothedspikerate;
    end
    
    imagedata=mean(imagedatas2,3);
    imagesc(flipud(imagedata));
    
    % Keeps the proportions of the data
            axis image;
            % no ticks
            set(gca,'xtick',[],'ytick',[]);
      xlimval=[0 size(imagedata,2)];
            ylimval=[0 size(imagedata,1)];
            
            distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
            
            text(distlabel(1),distlabel(2),sprintf('t %s c %s %s',num2str(cell2(1)),num2str(cell2(2)),num2str(change)),'FontSize',10,'Color','w');
          
            
            n=n+2
    end
    
  ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
        1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 0.99,sprintf('Placefield of cell pairs with 1 cell on barrier segment ripples <2cm/s \n %s for day %s ',...
        animals{1,1}, num2str(day)),'HorizontalAlignment','center','VerticalAlignment', 'top');
  
    cd(strcat(f_r_b.animal{1,2},'Plot/'));
figurename = sprintf('Cellpair_plf_1s_%s_%d',animals{1,1},day);

set(gcf,'PaperSize',[8 20]);

saveas(gcf, figurename, 'pdf');
close;
  % plot cell pairs with one placefield on barrier segment <2cm/s
figure;
 set(gcf,'position',[0 0 500 2000]);
 set(gcf,'PaperPositionMode','auto');
% number of rows
nrow2s=sum(isnan(goodcellpairs2s(:,5)));
n=1
  while n<=nrow2s
    
    % lookup indices of placefield for each cell
    cell1=goodcellpairs2s(n,1:2);
    cell2=goodcellpairs2s(n,3:4);
    change=goodcellpairs2s(n,5);
    % find place fields matching the cells
    cell1index=find(daydata(:,1)==cell1(1,1) & daydata(:,2)==cell1(1,2));
    cell1epochindex=epochindex(cell1index,:);
    cell2index=find(daydata(:,1)==cell2(1,1) & daydata(:,2)==cell2(1,2));
    cell2epochindex=epochindex(cell2index,:);
    subwidth=6;
    subplot(nrow2s/2,2,n-1+1)
     imagedatas1=[];
     
    for k=1:size(cell1index,1)
        imagedatas1(:,:,k)=outdata{1,cell1index(k)}.smoothedspikerate;
    end
    
    imagedata=mean(imagedatas1,3);
    imagesc(flipud(imagedata));
    
    % Keeps the proportions of the data
            axis image;
            % no ticks
            set(gca,'xtick',[],'ytick',[]);
      xlimval=[0 size(imagedata,2)];
            ylimval=[0 size(imagedata,1)];
            
            distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
            
            text(distlabel(1),distlabel(2),sprintf('t %s c %s %s',num2str(cell1(1)),num2str(cell1(2)),num2str(change)),'FontSize',8,'Color','w');
            
    subplot(nrow2s/2,2,n-1+2)
    imagedatas2=[];
    
    for k=1:size(cell2index,1)
        imagedatas2(:,:,k)=outdata{1,cell2index(k)}.smoothedspikerate;
    end
    
    imagedata=mean(imagedatas2,3);
    imagesc(flipud(imagedata));
    
    % Keeps the proportions of the data
            axis image;
            % no ticks
            set(gca,'xtick',[],'ytick',[]);
      xlimval=[0 size(imagedata,2)];
            ylimval=[0 size(imagedata,1)];
            
            distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.1*(ylimval(2)-ylimval(1))+ylimval(1)];
            
            text(distlabel(1),distlabel(2),sprintf('t %s c %s %s',num2str(cell2(1)),num2str(cell2(2)),num2str(change)),'FontSize',10,'Color','w');
          
            
            n=n+2
    end
    
  ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 ...
        1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
    text(0.5, 0.99,sprintf('Placefield of cell pairs with 1 cell on barrier segment ripples <2cm/s \n %s for day %s ',...
        animals{1,1}, num2str(day)),'HorizontalAlignment','center','VerticalAlignment', 'top');    
    
    cd(strcat(f_r_b.animal{1,2},'Plot/'));
figurename = sprintf('Cellpair_plf_2s_%s_%d',animals{1,1},day);

set(gcf,'PaperSize',[8 20]);

saveas(gcf, figurename, 'pdf');
close;
    
goodcellindex2=ismember(f_s_b.output{1,1}(1,i).cellpairindex,commoncells,'rows');
goodcellindex3=ismember(f_r_nb.output{1,1}(1,i).cellpairindex,commoncells,'rows');
goodcellindex4=ismember(f_r_b.output{1,1}(1,i).cellpairindex,commoncells,'rows');