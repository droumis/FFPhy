%% DEFINE EPOCHS AND ANIMALS

%animal selection
animals = {'Corriander','Eight','Five','Miles','Ten'};

% Epoch selection
epochfilterA = []; epochfilterB = [];
epochfilterA{1} = ['($experimentday == 4) & $dailyexposure == 1 & isequal($description,''TrackA'')'];
epochfilterB{1} = ['($experimentday == 4) & $dailyexposure == 1 & isequal($description,''TrackB'')'];

% Time selection
timefilter = {};

%% RUN FILTER TO LOOK AT PROPORTION ACTIVE AND RIPPLE SIZE FOR CA1
cellfilter = '($meanrate<7) & isequal($area,''CA1'')';

%Define iterator
iterator = 'multicellanal';

timefilter = [];

%create training data by calulating the linearized rates of all cells
f = createfilter('animal',animals,'epochs',epochfilterA,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);

f = setfilterfunction(f, 'calcripplespikingvelocity', {'spikes','ripples','pos','cellinfo'});
f = runfilter(f);

withinday_spikeCA1A = f;
save('/data13/mcarr/VelocityPaper/withinday_spikeCA1A.mat','withinday_spikeCA1A')

f = createfilter('animal',animals,'epochs',epochfilterB,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);

f = setfilterfunction(f, 'calcripplespikingvelocity', {'spikes','ripples','pos','cellinfo'});
f = runfilter(f);

withinday_spikeCA1B = f;
save('/data13/mcarr/VelocityPaper/withinday_spikeCA1B.mat','withinday_spikeCA1B')

%% RUN FILTER TO LOOK AT PROPORTION ACTIVE AND RIPPLE SIZE FOR CA3
cellfilter = '($meanrate<7) & isequal($area,''CA3'')';

%Define iterator
iterator = 'multicellanal';

timefilter = [];

%create training data by calulating the linearized rates of all cells
f = createfilter('animal',animals,'epochs',epochfilterA,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);

f = setfilterfunction(f, 'calcripplespikingvelocity', {'spikes','ripples','pos','cellinfo'});
f = runfilter(f);

withinday_spikeCA3A = f;
save('/data13/mcarr/VelocityPaper/withinday_spikeCA3A.mat','withinday_spikeCA3A')

f = createfilter('animal',animals,'epochs',epochfilterB,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);

f = setfilterfunction(f, 'calcripplespikingvelocity', {'spikes','ripples','pos','cellinfo'});
f = runfilter(f);

withinday_spikeCA3B = f;
save('/data13/mcarr/VelocityPaper/withinday_spikeCA3B.mat','withinday_spikeCA3B')

%% LOAD DATA
load('/data13/mcarr/VelocityPaper/withinday_spikeCA1A.mat')
load('/data13/mcarr/VelocityPaper/withinday_spikeCA1B.mat')
load('/data13/mcarr/VelocityPaper/withinday_spikeCA3A.mat')
load('/data13/mcarr/VelocityPaper/withinday_spikeCA3B.mat')

ca1_A = withinday_spikeCA1A; clear withinday_spikeCA1A
ca1_B = withinday_spikeCA1B; clear withinday_spikeCA1B
ca3_A = withinday_spikeCA3A; clear withinday_spikeCA3A
ca3_B = withinday_spikeCA3B; clear withinday_spikeCA3B

%% PLOT ACTIVATION PROBABILTY
bin = [1/2 1 2 4 8 16 32];
day = 1;
ca1_a = []; ca1_b = []; ca3_a = []; ca3_b = [];    
ca1totalA = []; ca1totalB = []; ca3totalA = []; ca3totalB = [];
for an = 1:length(ca1_A)
    if ~isempty(ca1_A(an).output{day}) && ~isempty(ca1_B(an).output{day})
        if sum(sum(isnan(ca1_A(an).output{day}.active_prob))>0)<2
            ca1_a = [ca1_a; ca1_A(an).output{day}(1).active_prob];
            ca1totalA = [ca1totalA; ca1_A(an).output{day}(1).total_proportion_active];
        end
        if sum(sum(isnan(ca1_B(an).output{day}.active_prob))>0)<2
            
        	ca1_b = [ca1_b; ca1_B(an).output{day}(1).active_prob];
            ca1totalB = [ca1totalB; ca1_B(an).output{day}(1).total_proportion_active];
        end
    
        if any(an == [1 4 5])
            if sum(sum(isnan(ca3_A(an).output{day}.active_prob))>0)<2
            	ca3_a = [ca3_a; ca3_A(an).output{day}(1).active_prob];
                ca3totalA = [ca3totalA; ca3_A(an).output{day}(1).total_proportion_active];
            end
            if sum(sum(isnan(ca3_B(an).output{day}.active_prob))>0)<2
            	ca3_b = [ca3_b; ca3_B(an).output{day}(1).active_prob];
                ca3totalB = [ca3totalB; ca3_B(an).output{day}(1).total_proportion_active];
            end
        end  
    end
end

subsf1 = lookup(repmat(bin,1,size(ca1_a,1)),bin);
f1 = reshape(ca1_a',size(ca1_a,1)*size(ca1_a,2),1);
subsf1 = subsf1 - mean(subsf1);
f1 = f1 - mean(f1);
b_1f = regress(f1,[ones(length(subsf1),1) subsf1]);

subsn1 = lookup(repmat(bin,1,size(ca1_b,1)),bin);
n1 = reshape(ca1_b',size(ca1_b,1)*size(ca1_b,2),1);
subsn1 = subsn1 - mean(subsn1);
n1 = n1 - mean(n1);
b_1n = regress(n1,[ones(length(subsn1),1) subsn1]);
Q = b_1n(2) - b_1f(2);

nperm = 1000;
qperm = nan(nperm,1);
z = [subsf1 f1; subsn1 n1];
Nn1 = length(subsn1); Nf1 = length(subsf1);
for s = 1:nperm
    ind = randperm(length(z));
    zperm = z(ind,:);
    tmpn = regress(zperm(1:Nn1,2),[ones(Nn1,1) zperm(1:Nn1,1)]);
    tmpf = regress(zperm(end-Nf1+1:end,2),[ones(Nf1,1) zperm(end-Nf1+1:end,1)]);
    qperm(s) = tmpn(2)-tmpf(2);
end
pvalue_1 = mean(abs(qperm) > abs(Q));

subsf3 = lookup(repmat(bin,1,size(ca3_a,1)),bin);
f3 = reshape(ca3_a',size(ca3_a,1)*size(ca3_a,2),1);
subsf3 = subsf3 - mean(subsf3);
f3 = f3 - mean(f3);
b_3f = regress(f3,[ones(length(subsf3),1) subsf3]);

subsn3 = lookup(repmat(bin,1,size(ca3_b,1)),bin);
n3 = reshape(ca3_b',size(ca3_b,1)*size(ca3_b,2),1);
subsn3 = subsn3 - mean(subsn3);
n3 = n3 - mean(n3);
b_3n = regress(n3,[ones(length(subsn3),1) subsn3]);
Q = b_3f(2) - b_3n(2);

nperm = 1000;
qperm = nan(nperm,1);
z = [subsf3 f3; subsn3 n3];
Nn3 = length(subsn3); Nf3 = length(subsf3);
for s = 1:nperm
    ind = randperm(length(z));
    zperm = z(ind,:);
    tmpn = regress(zperm(1:Nn3,2),[ones(Nn3,1) zperm(1:Nn3,1)]);
    tmpf = regress(zperm(end-Nf3+1:end,2),[ones(Nf3,1) zperm(end-Nf3+1:end,1)]);
    qperm(s) = tmpf(2)-tmpn(2);
end
pvalue_3 = mean(abs(qperm) > abs(Q));

%% DETERMINE WITHIN DAY RIPPLE DIFFERENCES

maxthreshF = []; maxthreshN = []; speedF = []; speedN = [];
day = 1;

for an = 1:length(ca1_A)
    if ~isempty(ca1_A(an).output{day}) && ~isempty(ca1_B(an).output{day})
        maxthreshF = [maxthreshF; ca1_A(an).output{day}(1).maxthresh'];
        speedF = [speedF; ca1_A(an).output{day}(1).speed];          
        maxthreshN = [maxthreshN; ca1_B(an).output{day}(1).maxthresh'];
    	speedN = [speedN; ca1_B(an).output{day}(1).speed];
	end
end

speedF = log(speedF); speedN = log(speedN);
maxthreshF(isinf(speedF) | isnan(speedF) | speedF < log(1/8)) = [];
speedF(isinf(speedF) | isnan(speedF) | speedF < log(1/8)) = [];
speedF(maxthreshF>12)=[]; maxthreshF(maxthreshF>12)=[];

maxthreshN(isinf(speedN) | isnan(speedN) | speedN < log(1/8)) = [];
speedN(isinf(speedN) | isnan(speedN) | speedN < log(1/8)) = [];
speedN(maxthreshN>12)=[]; maxthreshN(maxthreshN>12)=[];

%Determine Significance
maxthreshF = maxthreshF - mean(maxthreshF);
speedF = speedF - mean(speedF);
b_1f = regress(maxthreshF,[ones(length(speedF),1) speedF]);

maxthreshN = maxthreshN - mean(maxthreshN);
speedN = speedN - mean(speedN);
b_1n = regress(maxthreshN,[ones(length(speedN),1) speedN]);

Q = b_1n(2) - b_1f(2);

nperm = 1000;
qperm = nan(nperm,1);
z = [speedF maxthreshF; speedN maxthreshN];
Nn1 = length(speedN); Nf1 = length(speedF);
for s = 1:nperm
    ind = randperm(length(z));
    zperm = z(ind,:);
    tmpn = regress(zperm(1:Nn1,2),[ones(Nn1,1) zperm(1:Nn1,1)]);
    tmpf = regress(zperm(end-Nf1+1:end,2),[ones(Nf1,1) zperm(end-Nf1+1:end,1)]);
    qperm(s) = tmpn(2)-tmpf(2);
end
pvalue_1 = mean(abs(qperm) > abs(Q));

%Plot
[b_1f Lf Uf] = regress_boot(speedF,maxthreshF);
[b_1n Ln Un] = regress_boot(speedN,maxthreshN);

figure
hold on
plot([1 2],[b_1n b_1f],'k')
errorbar2([1 2],[b_1n b_1f],[Ln-b_1n Lf-b_1f; Un-b_1n Uf-b_1f],0.001,'k')
set(gca,'xtick',[1 2],'xticklabel','Novel|Familiar','xlim',[0.8 2.2])
set(gca,'yLim',[-0.35 0.05],'ytick',-0.3:0.05:0.1,'FontSize',18)
ylabel('Slope','FontSize',18)
title({'Within Day Comparison of Slope' , 'Speed vs. Ripple Amplitude'})
box off

%Novel Slope = -0.261 Familiar Slope = -0.046
%Significantly Different, p = 0
%Individual Animals: p = [0 0.001 0.6 0.07 0.01]
%For animal 3: Novel is significantly different from 0 but Familiar is not

%% DETERMINE WITHIN DAY DIFFERENCES IN RIPPLE AMPLITUDE
bin = log([1/4 1 4 16]);
maxthreshF = []; maxthreshN = []; speedF = []; speedN = [];
day = 1;

for an = 1%1:length(ca1_A)
    if ~isempty(ca1_A(an).output{day}) && ~isempty(ca1_B(an).output{day})
        maxthreshF = [maxthreshF; ca1_A(an).output{day}(1).totalmaxthresh'];
        speedF = [speedF; ca1_A(an).output{day}(1).speed];
        maxthreshN = [maxthreshN; ca1_B(an).output{day}(1).totalmaxthresh'];
        speedN = [speedN; ca1_B(an).output{day}(1).speed];
	end
end

speedF = log(speedF); speedN = log(speedN);
maxthreshF(isinf(speedF) | isnan(speedF) | speedF < log(1/8)) = [];
speedF(isinf(speedF) | isnan(speedF) | speedF < log(1/8)) = [];
speedF(maxthreshF>12) = []; maxthreshF(maxthreshF>12) = [];
subsF = lookup(speedF,bin);

maxthreshN(isinf(speedN) | isnan(speedN) | speedN < log(1/8)) = [];
speedN(isinf(speedN) | isnan(speedN) | speedN < log(1/8)) = [];
speedN(maxthreshN>12) = []; maxthreshN(maxthreshN>12) = [];
subsN = lookup(speedN,bin);

nboot = 1000;
AbootN = nan(nboot,length(bin));
AbootF = nan(nboot,length(bin));
count = 30;

for s = 1:nboot
    subsF_boot = []; subsN_boot = [];
    maxthreshF_boot = []; maxthreshN_boot = [];
    
    boot = ceil(length(subsN)*rand(2*length(subsN),1));
    for q = 1:max(subsN)
        tmpboot = find(subsN(boot)==q);
        tmpboot = boot(tmpboot);
        tmp = subsN(tmpboot);
        tmp = tmp(1:count);
        subsN_boot = [subsN_boot; tmp];
        tmp = maxthreshN(tmpboot);
        tmp = tmp(1:count);
        maxthreshN_boot = [maxthreshN_boot; tmp];
        clear tmpboot tmp
    end
    boot = ceil(length(subsF)*rand(2*length(subsF),1));
    for q = 1:max(subsF)
        tmpboot = find(subsF(boot)==q);
        tmpboot = boot(tmpboot);
        tmp = subsF(tmpboot);
        tmp = tmp(1:count);
        subsF_boot = [subsF_boot; tmp];
        tmp = maxthreshF(tmpboot);
        tmp = tmp(1:count);
        maxthreshF_boot = [maxthreshF_boot; tmp];
        clear tmpboot tmp
    end
    AbootN(s,:) = accumarray(subsN_boot,maxthreshN_boot,[max(subsN) 1],@(x) nanmean(x),NaN);
    AbootF(s,:) = accumarray(subsF_boot,maxthreshF_boot,[max(subsF) 1],@(x) nanmean(x),NaN);
end

Q = nan(length(bin),1);
qperm =nan(nboot,length(bin));
pvalue = nan(length(bin),1);
for q = 1:length(bin)
	F = maxthreshF(subsF==q);
    N = maxthreshN(subsN==q);
    z = [F; N];
    Q(q) = mean(N)-mean(F);
    for s = 1:nboot
        ind = randperm(length(z));
        zperm = z(ind,:);
        qperm(s,q) = mean(zperm(1:length(N)))-mean(zperm(end-length(F)+1:end));
    end
    pvalue(q) = mean(abs(qperm(:,q)) > abs(Q(q)));
end

% Plot
figure
hold on
errorbar2(bin,mean(AbootN),std(AbootN),0.001,'r')
errorbar2(bin,mean(AbootF),std(AbootF),0.001,'k')
legend('Novel','Familiar')
plot(bin,mean(AbootN),'r-o','MarkerFace','r')
plot(bin,mean(AbootF),'k-o','MarkerFace','k')
set(gca,'xtick',bin,'xticklabel',exp(bin),'FontSize',18)
set(gca,'ylim',[0.5 2],'ytick',-0.5:0.5:2)
ylabel('Mean Ripple Power (Standard deviations from mean)')
xlabel('Speed cm/sec')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_ripplepower_withinday.pdf', m, d, y);
print('-dpdf', savestring)

%Significantly Different, p = 0
%Individual Animals: p = [{0 0.06 0.17 0} 0 0 {0.007 0.001 0.001 0} 0]

%% PLOT SLOPES FOR INDIVIDUAL ANIMALS
figure
hold on
title({'Within Day Comparison of Slope' , 'Speed vs. Ripple Amplitude'})

day = 1;
for an = 1:length(ca1_A)
    maxthreshF = []; maxthreshN = []; speedF = []; speedN = [];
    if ~isempty(ca1_A(an).output{day}) && ~isempty(ca1_B(an).output{day})
        maxthreshF = [maxthreshF; ca1_A(an).output{day}(1).maxthresh'];
        speedF = [speedF; ca1_A(an).output{day}(1).speed];          
        maxthreshN = [maxthreshN; ca1_B(an).output{day}(1).maxthresh'];
    	speedN = [speedN; ca1_B(an).output{day}(1).speed];
    end

    speedF = log(speedF); speedN = log(speedN);
    maxthreshF(isinf(speedF) | isnan(speedF) | speedF < log(1/8)) = [];
    speedF(isinf(speedF) | isnan(speedF) | speedF < log(1/8)) = [];
    speedF(maxthreshF>12)=[]; maxthreshF(maxthreshF>12)=[];

    maxthreshN(isinf(speedN) | isnan(speedN) | speedN < log(1/8)) = [];
    speedN(isinf(speedN) | isnan(speedN) | speedN < log(1/8)) = [];
    speedN(maxthreshN>12)=[]; maxthreshN(maxthreshN>12)=[];

    %Determine Significance
    maxthreshF = maxthreshF - mean(maxthreshF);
    speedF = speedF - mean(speedF);
    b_1f = regress(maxthreshF,[ones(length(speedF),1) speedF]);

    maxthreshN = maxthreshN - mean(maxthreshN);
    speedN = speedN - mean(speedN);
    b_1n = regress(maxthreshN,[ones(length(speedN),1) speedN]);

    Q = b_1n(2) - b_1f(2);

    nperm = 1000;
    qperm = nan(nperm,1);
    z = [speedF maxthreshF; speedN maxthreshN];
    Nn1 = length(speedN); Nf1 = length(speedF);
    for s = 1:nperm
        ind = randperm(length(z));
        zperm = z(ind,:);
        tmpn = regress(zperm(1:Nn1,2),[ones(Nn1,1) zperm(1:Nn1,1)]);
        tmpf = regress(zperm(end-Nf1+1:end,2),[ones(Nf1,1) zperm(end-Nf1+1:end,1)]);
        qperm(s) = tmpn(2)-tmpf(2);
    end
    pvalue = mean(abs(qperm) > abs(Q));
    
    [b_1f Lf Uf] = regress_boot(speedF,maxthreshF);
    [b_1n Ln Un] = regress_boot(speedN,maxthreshN);

    %Plot
    subplot(2,3,an)
    hold on
    plot([1 2],[b_1n b_1f],'k-o','MarkerFace','k')
    errorbar2([1 2],[b_1n b_1f],[Ln-b_1n Lf-b_1f; Un-b_1n Uf-b_1f],0.001,'k')
    set(gca,'xtick',[1 2],'xticklabel','Novel|Familiar','xlim',[0.8 2.2])
    set(gca,'yLim',[-0.55 0.2],'ytick',-0.6:0.1:0.1,'FontSize',18)
    ylabel('Slope','FontSize',18)
    title(num2str(pvalue))
    box off
end

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_individualanimal_withinday_rippleslope.pdf', m, d, y);
print('-dpdf', savestring)

