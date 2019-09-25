%% RUN FILTER FOR ALL CELLS
%Animal selection
animals = {'Bond','Frank','Ten'};

%Define epochs
epochfilter = [];
%epochfilter{1} = 'isequal($type, ''run'')';
epochfilter{1} = 'isequal($type, ''sleep'')';

cellfilter = '($meanrate<7)';

%Define iterator
iterator = 'multicellanal';

%timefilter = {{'get2dstate', '$velocity<4'}};
timefilter = {{'get2dstate', '$velocity<4 & $immobilitytime>60'}};

%Define iterator
iterator = 'multicellanal';

%create training data by calulating the linearized rates of all cells
f3 = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);
f3 = setfilterfunction(f3, 'calcswrrate', {'spikes','ripples','cellinfo'},'min_cells',5);
f3 = runfilter(f3);

%% PLOT
% out is a vector with the following fields:
% 1:    index
% 2:    number of ca1 cells
% 3:    SWR firing rate for ca1 place cells
% 4:    # of SWR spikes fired by ca1 place cells
% 5:    SWR firing rate for ca1 non-place cells
% 6:    # of SWR spikes fired by ca1 non-place cells
% 7:    pre-SWR firing rate for ca1 place cells
% 8:    # of pre-SWR spikes fired by ca1 place cells
% 9:    pre-SWR firing rate for ca1 non-place cells
% 10:    # of pre-SWR spikes fired by ca1 non-place cells
% 11:    number of ca3 cells
% 12:    SWR firing rate for ca3 place cells
% 13:    # of SWR spikes fired by ca3 place cells
% 14:    SWR firing rate for ca3 non-place cells
% 15:    # of SWR spikes fired by ca3 non-place cells
% 16:    pre-SWR firing rate for ca3 place cells
% 17:    # of pre-SWR spikes fired by ca3 place cells
% 18:    pre-SWR firing rate for ca3 non-place cells
% 19:    # of pre-SWR spikes fired by ca3 non-place cells


run = numericgroupcombine(f3);
sleep = numericgroupcombine(f3);

% Compare SWR firing rates for sleep and run, place cells and non place
% cells, ca1 and ca3
sleep1 = sleep(:,3)./sleep(:,2); invalid = isinf(sleep1)|isnan(sleep1); sleep1(invalid) = [];
sleep3 = sleep(:,12)./sleep(:,11); invalid = isinf(sleep3)|isnan(sleep3); sleep3(invalid) = [];
sleep1_non = sleep(:,5)./sleep(:,2); invalid = isinf(sleep1_non)|isnan(sleep1_non); sleep1_non(invalid) = [];
sleep3_non = sleep(:,14)./sleep(:,11); invalid = isinf(sleep3_non)|isnan(sleep3_non); sleep3_non(invalid) = [];
run1 = run(:,3)./run(:,2); invalid = isinf(run1)|isnan(run1); run1(invalid) = [];
run3 = run(:,12)./run(:,11); invalid = isinf(run3)|isnan(run3); run3(invalid) = [];
run1_non = run(:,5)./run(:,2); invalid = isinf(run1_non)|isnan(run1_non); run1_non(invalid) = [];
run3_non = run(:,14)./run(:,11); invalid = isinf(run3_non)|isnan(run3_non); run3_non(invalid) = [];

figure
hold on
bar(1:4,[mean(run1) mean(sleep1) mean(run1_non) mean(sleep1_non)],'b')
errorbar2(1:4,[mean(run1) mean(sleep1) mean(run1_non) mean(sleep1_non) ],...
    [std(run1)./sqrt(length(run1)-1) std(sleep1)./sqrt(length(sleep1)-1) std(run1_non)./sqrt(length(run1_non)-1) std(sleep1_non)./sqrt(length(sleep1_non)-1)],'k')
set(gca,'xtick',1:4,'xticklabel',[{'Awake pc'},{'Quiescent pc'},{'Awake non pc'},{'Quiescent non pc'}],'xlim',[0.5 4.5],'ylim',[0 3])

figure
hold on
bar(1:4,[mean(run3) mean(sleep3) mean(run3_non) mean(sleep3_non)],'b')
errorbar2(1:4,[mean(run3) mean(sleep3) mean(run3_non) mean(sleep3_non) ],...
    [std(run3)./sqrt(length(run3)-1) std(sleep3)./sqrt(length(sleep3)-1) std(run3_non)./sqrt(length(run3_non)-1) std(sleep3_non)./sqrt(length(sleep1_non)-1)],'k')
set(gca,'xticklabel',[{'Awake pc'},{'Quiescent pc'},{'Awake non pc'},{'Quiescent non pc'}],'xlim',[0.5 4.5],'ylim',[0 3])

%Compare pre-SWR firing rates for sleep and run, place cells and non place
%cells, ca1 and ca3
% Compare SWR firing rates for sleep and run, place cells and non place
% cells, ca1 and ca3
sleep1 = sleep(:,7)./sleep(:,2); invalid = isinf(sleep1)|isnan(sleep1); sleep1(invalid) = [];
sleep3 = sleep(:,16)./sleep(:,11); invalid = isinf(sleep3)|isnan(sleep3); sleep3(invalid) = [];
sleep1_non = sleep(:,9)./sleep(:,2); invalid = isinf(sleep1_non)|isnan(sleep1_non); sleep1_non(invalid) = [];
sleep3_non = sleep(:,18)./sleep(:,11); invalid = isinf(sleep3_non)|isnan(sleep3_non); sleep3_non(invalid) = [];
run1 = run(:,7)./run(:,2); invalid = isinf(run1)|isnan(run1); run1(invalid) = [];
run3 = run(:,16)./run(:,11); invalid = isinf(run3)|isnan(run3); run3(invalid) = [];
run1_non = run(:,9)./run(:,2); invalid = isinf(run1_non)|isnan(run1_non); run1_non(invalid) = [];
run3_non = run(:,18)./run(:,11); invalid = isinf(run3_non)|isnan(run3_non); run3_non(invalid) = [];

figure
hold on
bar(1:4,[mean(run1) mean(sleep1) mean(run1_non) mean(sleep1_non)],'b')
errorbar2(1:4,[mean(run1) mean(sleep1) mean(run1_non) mean(sleep1_non) ],...
    [std(run1)./sqrt(length(run1)-1) std(sleep1)./sqrt(length(sleep1)-1) std(run1_non)./sqrt(length(run1_non)-1) std(sleep1_non)./sqrt(length(sleep1_non)-1)],'k')
set(gca,'xtick',1:4,'xticklabel',[{'Awake pc'},{'Quiescent pc'},{'Awake non pc'},{'Quiescent non pc'}],'xlim',[0.5 4.5],'ylim',[0 0.55])
ylabel('Mean firing rate before SWR')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_preswr_firingrate_ca1.pdf', m, d, y);
print('-dpdf', savestring)

figure
hold on
bar(1:4,[mean(run3) mean(sleep3) mean(run3_non) mean(sleep3_non)],'b')
errorbar2(1:4,[mean(run3) mean(sleep3) mean(run3_non) mean(sleep3_non) ],...
    [std(run3)./sqrt(length(run3)-1) std(sleep3)./sqrt(length(sleep3)-1) std(run3_non)./sqrt(length(run3_non)-1) std(sleep3_non)./sqrt(length(sleep1_non)-1)],'k')
set(gca,'xticklabel',[{'Awake pc'},{'Quiescent pc'},{'Awake non pc'},{'Quiescent non pc'}],'xlim',[0.5 4.5],'ylim',[0 0.55])
ylabel('Mean firing rate before SWR')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_preswr_firingrate_ca3.pdf', m, d, y);
print('-dpdf', savestring)

%Test significance using anova
group_run = [ones(size(run1)); ones(size(run1_non)); 2*ones(size(sleep1)); 2*ones(size(sleep1_non))];
group_place = [ones(size(run1)); 2*ones(size(run1_non)); ones(size(sleep1)); 2*ones(size(sleep1_non))];

[p table stats] = anovan([run1; run1_non; sleep1;sleep1_non],{group_run group_place},'model','full');
c = multcompare(stats,'dimension',[1 2]);
%Main effect of place cells, p<1e-5, interaction between place cells and
%run/slee p<1e-5, no main effect of run/sleep