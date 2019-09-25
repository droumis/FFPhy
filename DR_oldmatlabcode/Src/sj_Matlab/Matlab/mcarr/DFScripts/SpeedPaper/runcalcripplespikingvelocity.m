%% DEFINE EPOCHS AND ANIMALS

%animal selection
animals = {'Conley','Corriander','Eight','Five','Miles','Ten'};

% Epoch selection
epochfilter = [];
for i = 1:10
    epochfilter{i} = ['($exposure == ',num2str(i),') & isequal($description,''TrackA'')'];
end

% Time selection
timefilter = {};

%% RUN FILTER TO LOOK AT PROPORTION ACTIVE AND RIPPLE SIZE FOR CA1
cellfilter = '($meanrate<7) & isequal($area,''CA1'')';

%Define iterator
iterator = 'multicellanal';

timefilter = [];

%create training data by calulating the linearized rates of all cells
f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);

f = setfilterfunction(f, 'calcripplespikingvelocity', {'spikes','ripples','pos','cellinfo'});
f = runfilter(f);

spikeCA1 = f;
save('/data13/mcarr/VelocityPaper/spikeCA1.mat','spikeCA1')

%% RUN FILTER TO LOOK AT PROPORTION ACTIVE AND RIPPLE SIZE FOR CA3
cellfilter = '($meanrate<7) & isequal($area,''CA3'')';

%Define iterator
iterator = 'multicellanal';

timefilter = [];

%create training data by calulating the linearized rates of all cells
f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);

f = setfilterfunction(f, 'calcripplespikingvelocity', {'spikes','ripples','pos','cellinfo'});
f = runfilter(f);

spikeCA3 = f;
save('/data13/mcarr/VelocityPaper/spikeCA3.mat','spikeCA3')

%% LOAD FILTERS
%CA1
load '/data13/mcarr/VelocityPaper/spikeCA1.mat'

%CA3
load '/data13/mcarr/VelocityPaper/spikeCA3.mat'

%% PLOT CA1 AND CA3 PROPORTION ACTIVE FOR FIRST EXPOSURE
bin = [1/2 1 2 4 8 16 32];

ca1 = []; ca3 = []; ca1total = []; ca3total = [];
day = 1;
for an = 1:length(spikeCA1)
    if ~isempty(spikeCA1(an).output)
        if sum(sum(isnan(spikeCA1(an).output{day}.active_prob))>0)<2
            ca1 = [ca1; spikeCA1(an).output{day}(1).active_prob];
            ca1total = [ca1total; spikeCA1(an).output{day}(1).total_proportion_active];
        end
    end
    if ~isempty(spikeCA3(an).output)
        if sum(sum(isnan(spikeCA3(an).output{day}.active_prob))>0)<2
            ca3 = [ca3; spikeCA3(an).output{day}(1).active_prob];
            ca3total = [ca3total; spikeCA3(an).output{day}(1).total_proportion_active];
        end
    end
end

subs1 = lookup(repmat(bin,1,size(ca1,1)),bin);
tmp1 = reshape(ca1',size(ca1,1)*size(ca1,2),1);

a = accumarray(subs1,tmp1,[length(bin) 1],@(x) nanmean(x),NaN);
Lstd = accumarray(subs1,tmp1,[length(bin) 1],@(x) nanstd(x),NaN);
Lsqrt = accumarray(subs1,tmp1,[length(bin) 1],@(x) sqrt(sum(~isnan(x))),NaN);
seL = -Lstd./Lsqrt;
seU = +Lstd./Lsqrt;
subs1 = repmat(bin,1,size(ca1,1));
[L U slope intercept] = regress_fill(log(subs1),tmp1,log(bin));

figure
fill(log([bin fliplr(bin)]), [L fliplr(U)],'r','FaceAlpha',1,'EdgeColor','none')
hold on
plot(log(bin),log(bin)*slope+intercept,'r')    
plot(log(bin),a,'ro','MarkerFace','r')
errorbar2(log(bin),a,[seL seU]',0.001,'r','plottype','semilogx')
set(gca,'xtick',log(bin),'xtickLabel',bin)
set(gca,'yLim',[-0.01 0.16],'ytick',0:0.05:.2,'FontSize',20)
xlabel('Speed (cm/sec)','FontSize',22)
ylabel('Activation Probability','FontSize',22)
box off

subs3 = lookup(repmat(bin,1,size(ca3,1)),bin);
tmp3 = reshape(ca3',size(ca3,1)*size(ca3,2),1);

a = accumarray(subs3,tmp3,[length(bin) 1],@(x) nanmean(x),NaN);
Lstd = accumarray(subs3,tmp3,[length(bin) 1],@(x) nanstd(x),NaN);
Lsqrt = accumarray(subs3,tmp3,[length(bin) 1],@(x) sqrt(sum(~isnan(x))),NaN);
seL = -Lstd./Lsqrt;
seU = +Lstd./Lsqrt;
subs3 = repmat(bin,1,size(ca3,1));
b = regress(tmp3,[ones(length(subs3),1) log(subs3)'],0.05);

[L U slope_L slope_U] = regress_fill(log(subs3),tmp3,log(bin));

fill(log([bin fliplr(bin)]), [L fliplr(U)],'k','FaceAlpha',1,'EdgeColor','none')
hold on
plot(log(bin),log(bin)*b(2)+b(1),'k')    
plot(log(bin),a,'ko','MarkerFace','k')
errorbar2(log(bin),a,[seL seU]',0.001,'k','plottype','semilogx')
set(gca,'xtick',log(bin),'xtickLabel',bin)
set(gca,'yLim',[-0.02 0.16],'ytick',0:0.05:0.2,'FontSize',20)
xlabel('Speed (cm/sec)','FontSize',22)
ylabel('Activation Probability','FontSize',22)
box off


% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_velocity_activationprobability_exposure1.pdf', m, d, y);
print('-dpdf', savestring)

%% PLOT CA1 RIPPLE SIZE FOR FIRST EXPOSURE
bin = log([1/4 1/2 1 2 4 8 16 32]);

maxthresh = []; speed = [];
day = 1;
for an = 1:length(spikeCA1)
    if ~isempty(spikeCA1(an).output)
        speed = [speed; spikeCA1(an).output{day}(1).speed];
        maxthresh = [maxthresh; spikeCA1(an).output{day}(1).maxthresh'];
    end
end

speed = log(speed);
maxthresh(isinf(speed) | isnan(speed) | speed < log(1/8)) = [];
maxthresh = maxthresh - 3;
speed(isinf(speed) | isnan(speed) | speed < log(1/8)) = [];
subs = lookup(speed,bin);

a = accumarray(subs, maxthresh,[length(bin) 1],@(x) nanmean(x),NaN);
Lstd = accumarray(subs,maxthresh,[length(bin) 1],@(x) nanstd(x),NaN);
Lsqrt = accumarray(subs,maxthresh,[length(bin) 1],@(x) sqrt(sum(~isnan(x))),NaN);
seL = -Lstd./Lsqrt;
seU = +Lstd./Lsqrt;
b = regress(maxthresh,[ones(length(speed),1) speed],0.05);

[L U] = regress_fill(speed,maxthresh,bin);

figure
fill([bin fliplr(bin)], [L fliplr(U)],'k','FaceAlpha',1,'EdgeColor','none')
hold on
plot(bin,bin*b(2)+b(1),'k')    
plot(bin,a,'ko','MarkerFace','k')
errorbar2(bin,a,[seL seU]',0.001,'k','plottype','semilogx')
set(gca,'xtick',bin,'xtickLabel',exp(bin))
set(gca,'yLim',[0 3],'ytick',0:0.5:3.5,'FontSize',20)
xlabel('Speed (cm/sec)','FontSize',22)
ylabel('Normalized Ripple Power','FontSize',22)
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_velocity_ripplepower_exposure1.pdf', m, d, y);
print('-dpdf', savestring)

%% PLOT CA1 AND CA3 PROPORTION ACTIVE OVER DAYS
bin = [1/2 1 2 4 8 16 32];

ca1 = cell(10,1); ca3 = cell(10,1);
ca1total = cell(10,1); ca3total = cell(10,1);

for an = 1:length(spikeCA1)
    for day = 1:length(spikeCA1(an).output)
        if ~isempty(spikeCA1(an).output)
            if sum(sum(isnan(spikeCA1(an).output{day}.active_prob))>0)<1
                ca1{day} = [ca1{day}; spikeCA1(an).output{day}(1).active_prob];
                ca1total{day} = [ca1total{day}; spikeCA1(an).output{day}(1).total_proportion_active];
            end
        end
        if ~isempty(spikeCA3(an).output)
            if sum(sum(isnan(spikeCA3(an).output{day}.active_prob))>0)<1
                ca3{day} = [ca3{day}; spikeCA3(an).output{day}(1).active_prob];
                ca3total{day} = [ca3total{day}; spikeCA3(an).output{day}(1).total_proportion_active];
            end
        end
    end
end

B1 = nan(10,1); L1 = nan(10,1); U1 = nan(10,1);
B3 = nan(10,1); L3 = nan(10,1); U3 = nan(10,1);
pvalue_1 = nan(10,1); pvalue_3 =nan(10,1);

for day = 1:10
    tmp1 = reshape(ca1{day}',size(ca1{day},1)*size(ca1{day},2),1);
    subs1 = log(repmat(bin,1,size(ca1{day},1)));
    subs1(isnan(tmp1))=[]; tmp1(isnan(tmp1))=[];
    tmp3 = reshape(ca3{day}',size(ca3{day},1)*size(ca3{day},2),1);
    subs3 = log(repmat(bin,1,size(ca3{day},1)));
    subs3(isnan(tmp3))=[]; tmp3(isnan(tmp3))=[];
    
    [b slope_L slope_U] = regress_boot(subs1,tmp1);
    B1(day) = b; L1(day) = slope_L; U1(day) = slope_U;
    pvalue_1(day) = regress_permutation(subs1,tmp1);

    [b slope_L slope_U] = regress_boot(subs3,tmp3);
    pvalue_3(day) = regress_permutation(subs3,tmp3);
    B3(day) = b; L3(day) = slope_L; U3(day) = slope_U;
end

% Plot slope over days
figure
hold on
plot(1:10,B1,'k-o','MarkerFace','k')
errorbar2(1:10,B1,[L1-B1 U1-B1]',0.2,'k')
plot(1:10,B3,'r-o','MarkerFace','r')
errorbar2(1:10,B3,[L3-B3 U3-B3]',0.2,'r')
xlabel('Exposure to Environment','FontSize',24)
ylabel('Slope of Speed vs. Activation Probability','FontSize',24)
set(gca,'xtick',1:10,'xlim',[0 11])
set(gca,'ylim',[-0.03 0.005],'ytick',-0.03:0.005:0.005,'FontSize',16)
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_velocity_activationprobability_overdays.pdf', m, d, y);
print('-dpdf', savestring)


% Plot total activation over days
A1 = nan(10,2); A3 = nan(10,2);
ca1_total = []; ca3_total = [];
for day = 1:10
    A1(day,1) = mean(ca1total{day});
    A1(day,2) = std(ca1total{day})./sqrt(length(ca1total{day})+1);
    ca1_total = [ca1_total; [ones(length(ca1total{day}),1).*day ca1total{day}]];

    A3(day,1) = mean(ca3total{day});
    A3(day,2) = std(ca3total{day})./sqrt(length(ca3total{day})+1);
    ca3_total = [ca3_total; [ones(length(ca3total{day}),1).*day ca3total{day}]];
end

%Determine significance
pvalue1 = regress_permutation(ca1_total(:,1),ca1_total(:,2))
pvalue3 = regress_permutation(ca3_total(:,1),ca3_total(:,2))

figure
hold on
plot(1:10,A1(:,1),'k-o','MarkerFace','k')
plot(1:10,A3(:,1),'r-o','MarkerFace','r')
errorbar2(1:10,A1(:,1),A1(:,2),0.001,'k')
errorbar2(1:10,A3(:,1),A3(:,2),0.001,'r')
xlabel('Exposure to Environment','FontSize',24)
ylabel('Activation Probability during SWRs','FontSize',24)
set(gca,'xtick',1:10,'xlim',[0 11])
set(gca,'ylim',[0 0.15],'ytick',0:0.025:0.15,'FontSize',16)
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_velocity_totalactivationprobability_overdays.pdf', m, d, y);
print('-dpdf', savestring)

%Pvalues for all animals:
%CA1 = 0 CA3 = 0.07
%Pvalues for individual animals:
%CA1 is significant for all animals with sufficient cells
%(i.e. cells recorded on each day)
%Conley:    CA1: 0.27   CA3: 0.77
%Corriander: CA1: 0.024 CA3: 0.16
%Eight:     CA1:  0
%Five:      CA1: 0.003
%Miles:     CA1: 0.331  CA3: 0.736
%Ten:       CA1: 0.049   CA3: 0.079

%% PLOT CA1 RIPPLE SIZE OVER DAYS
bin = log([1/4 1/2 1 2 4 8 16 32]);

maxthresh = cell(10,1); speed = cell(10,1);

for an = 1:length(spikeCA1)
    for day = 1:length(spikeCA1(an).output)
        if ~isempty(spikeCA1(an).output{day})
            speed{day} = [speed{day}; spikeCA1(an).output{day}(1).speed];
            maxthresh{day} = [maxthresh{day}; spikeCA1(an).output{day}(1).maxthresh'];
        end
    end
end

nboot = 5000;
qboot = nan(nboot,10);
for day = 1:10
    speed{day} = log(speed{day});
    maxthresh{day}(isinf(speed{day}) | isnan(speed{day}) | speed{day} < log(1/8)) = [];
    maxthresh{day} = maxthresh{day} - 3;
    speed{day}(isinf(speed{day}) | isnan(speed{day}) | speed{day} < log(1/8)) = [];
    speed{day}(isnan(maxthresh{day})) = [];
    maxthresh{day}(isnan(maxthresh{day})) = [];
    
    for s = 1:nboot
        boot = ceil(length(speed{day})*rand(length(speed{day}),1));
        xboot = speed{day}(boot);
        yboot = maxthresh{day}(boot);
        b = regress(yboot,[ones(length(xboot),1) xboot],0.05);
        qboot(s,day) = b(2);
    end
end

a = mean(qboot);
seL = a - prctile(qboot,2.5);
seU = a - prctile(qboot,97.5);

figure
hold on
plot(1:10,a,'k-o','MarkerFace','k')
errorbar2(1:10,a,[seL; seU],0.001,'k')
set(gca,'xtick',1:10,'xtickLabel',1:10)
set(gca,'yLim',[-0.55 0.05],'ytick',-0.5:0.1:0,'FontSize',20)
xlabel('Exposure to Environment','FontSize',22)
ylabel('Normalized Ripple Power','FontSize',22)
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_velocity_slopevsripplepower_overdays.pdf', m, d, y);
print('-dpdf', savestring)