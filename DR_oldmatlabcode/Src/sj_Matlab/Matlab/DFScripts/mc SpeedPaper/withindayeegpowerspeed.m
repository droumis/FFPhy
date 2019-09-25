%% DEFINE COMMON FILTER CRITERIA

%animal selection
animals = {'Eight','Five','Miles','Ten'};

% Epoch selection
epochfilterA = []; epochfilterB = [];
epochfilterA{1} = ['($experimentday == 4) & $dailyexposure == 1 & isequal($description,''TrackA'')'];
epochfilterB{1} = ['($experimentday == 4) & $dailyexposure == 1 & isequal($description,''TrackB'')'];

% Time selection
timefilter = {};

% Tetrode selection
ca1tetfilter =  '(isequal($area, ''CA1'') & $numcells>2 )';

%Select iterator
iterator = 'epocheegnonreferenceanal';

%% RUN FILTER TO LOOK AT DIFFERENCE IN SLOPE, NOVEL VS. FAMILIAR

% Create and Run Filter A
f = createfilter('animal',animals,'epochs',epochfilterA,'excludetime',timefilter, 'eegtetrodes',ca1tetfilter,'iterator', iterator);
f = setfilterfunction(f,'calceegpowerspeed',{'eeg','pos'});
f = runfilter(f);
withinday_eegA = f;

save('/data13/mcarr/VelocityPaper/withinday_eegA.mat','withinday_eegA')

% Create and Run Filter B
f = extendfilter(f,'epochs',epochfilterB,'iterator','epocheegnonreferenceoutputanal','keepoutput','input');
f = runfilter(f);
withinday_eegB = f;

save('/data13/mcarr/VelocityPaper/withinday_eegB.mat','withinday_eegB')

%clear all
%% LOAD FILTER
load '/data13/mcarr/VelocityPaper/withinday_eegA.mat'
load '/data13/mcarr/VelocityPaper/withinday_eegB.mat'

eegA = withinday_eegA; clear withinday_eegA
eegB = withinday_eegB; clear withinday_eegB

%% COMPARE NOVEL AND FAMILIAR WITHIN DAY
bin = [1/4 1 4 16];

lowN = []; lowF = []; highN = []; highF = []; speedN = []; speedF = [];

for an =1:length(eegA)
        tmp_lowN = []; tmp_highN = []; tmp_lowF = []; tmp_highF = [];
        speedF = [speedF; eegA(an).output{1}(1).speed];
        speedN = [speedN; eegB(an).output{1}(1).speed];

        for tet = 1:length(eegA(an).output{1})
            tmp_lowF = [tmp_lowF eegA(an).output{1}(tet).lowgamma_power];
            tmp_highF = [tmp_highF eegA(an).output{1}(tet).highgamma_power];
        end
        tmp_lowF = mean(tmp_lowF,2);
        tmp_highF = mean(tmp_highF,2);
        
        for tet = 1:length(eegB(an).output{1})
            tmp_lowN = [tmp_lowN eegB(an).output{1}(tet).lowgamma_power];
            tmp_highN = [tmp_highN eegB(an).output{1}(tet).highgamma_power];
        end
        tmp_lowN = mean(tmp_lowN,2);
        tmp_highN = mean(tmp_highN,2);
  
        lowF = [lowF; (tmp_lowF - mean(tmp_lowF))./std(tmp_lowF)];
        highF = [highF; (tmp_highF - mean(tmp_highF))./std(tmp_highF)];
          
        lowN = [lowN; (tmp_lowN - mean(tmp_lowF))./std(tmp_lowF)];
        highN = [highN; (tmp_highN - mean(tmp_highF))./std(tmp_highF)];
end

invalid = speedN<1/8 | lowN>12 | highN >12;
speedN(invalid)=[]; highN(invalid)=[]; lowN(invalid)=[];
speedN = log(speedN);
subsN = lookup(speedN,log(bin));

invalid = speedF<1/8 | lowF>10 | highF >10;
speedF(invalid)=[]; highF(invalid)=[]; lowF(invalid)=[];
speedF = log(speedF);
subsF = lookup(speedF,log(bin));

% Plot Bar graph and determine significance
AF_low = accumarray(subsF,lowF,[length(bin) 1],@(x) nanmean(x),NaN);
AN_low = accumarray(subsN,lowN,[length(bin) 1],@(x) nanmean(x),NaN);
EF_low = accumarray(subsF,lowF,[length(bin) 1],@(x) nanstd(x)./sqrt(sum(~isnan(x))-1),NaN);
EN_low = accumarray(subsN,lowN,[length(bin) 1],@(x) nanstd(x)./sqrt(sum(~isnan(x))-1),NaN);

AF_high = accumarray(subsF,highF,[length(bin) 1],@(x) nanmean(x),NaN);
AN_high = accumarray(subsN,highN,[length(bin) 1],@(x) nanmean(x),NaN);
EF_high = accumarray(subsF,highF,[length(bin) 1],@(x) nanstd(x)./sqrt(sum(~isnan(x))-1),NaN);
EN_high = accumarray(subsN,highN,[length(bin) 1],@(x) nanstd(x)./sqrt(sum(~isnan(x))-1),NaN);

figure
barerror([1 3 5 7],AN_low,EN_low,0.5,0.5,'r');
barerror([2 4 6 8],AF_low,EF_low,0.5,0.5);
set(gca,'xlim',[0 9],'ylim',[-0.5 2],'xtick',1.5:2:7.5,'xticklabel',[1/4 1 4 16])
xlabel('Speed')
ylabel('Normalized Slow Gamma Power')

figure
barerror([1 3 5 7],AN_high,EN_high,0.5,0.5,'r');
barerror([2 4 6 8],AF_high,EF_high,0.5,0.5);
set(gca,'xlim',[0 9],'ylim',[-0.75 1.5],'xtick',1.5:2:7.5,'xticklabel',[1/4 1 4 16])
xlabel('Speed')
ylabel('Normalized Fast Gamma Power')

pvalue = nan(length(bin),2);
for i = 1:length(bin)
    pvalue(i,1) = ranksum(lowF(subsF==i),lowN(subsN==i));
    pvalue(i,2) = ranksum(highF(subsF==i),highN(subsN==i));
end

%Save Slow Gamma Figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_withinday_slowgamma_power.pdf', m, d, y);
print('-dpdf', savestring)

%Save Fast Gamma Figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_withinday_fastgamma_power.pdf', m, d, y);
print('-dpdf', savestring)

%% Compare Slopes

nboot = 1000;
count = 100;
BN_low = nan(nboot,2);  BF_low = nan(nboot,2);
BN_high = nan(nboot,2); BF_high = nan(nboot,2);

for s = 1:nboot
    speedF_boot = []; lowF_boot = []; highF_boot = [];
    speedN_boot = []; lowN_boot = []; highN_boot = [];

    bootF = ceil(length(subsF)*rand(2*length(subsF),1));
    bootN = ceil(length(subsN)*rand(2*length(subsN),1));
    
    for q = 1:max(subsN)
        tmpboot = bootF(subsF(bootF)==q);
        tmp = speedF(tmpboot);
        tmp = tmp(1:count);
        speedF_boot = [speedF_boot; tmp];
        tmp = lowF(tmpboot);
        tmp = tmp(1:count);
        lowF_boot = [lowF_boot; tmp];
        tmp = highF(tmpboot);
        tmp = tmp(1:count);
        highF_boot = [highF_boot; tmp];
        
        
        tmpboot = bootN(subsN(bootN)==q);
        tmp = speedN(tmpboot);
        tmp = tmp(1:count);
        speedN_boot = [speedN_boot; tmp];
        tmp = lowN(tmpboot);
        tmp = tmp(1:count);
        lowN_boot = [lowN_boot; tmp];
        tmp = highN(tmpboot);
        tmp = tmp(1:count);
        highN_boot = [highN_boot; tmp];
        clear tmpboot tmp
    end
    
    %Recenter data at origin to isolate changes in slope
    speedF_boot = speedF_boot - mean(speedF_boot);
    lowF_boot = lowF_boot - mean(lowF_boot);
    highF_boot = highF_boot - mean(highF_boot);
    
    speedN_boot = speedN_boot - mean(speedN_boot);
    lowN_boot = lowN_boot - mean(lowN_boot);
    highN_boot = highN_boot - mean(highN_boot);
    
    BF_low(s,:) = regress(lowF_boot,[ones(length(speedF_boot),1) speedF_boot]);
    BF_high(s,:) = regress(highF_boot,[ones(length(speedF_boot),1) speedF_boot]);
       
    BN_low(s,:) = regress(lowN_boot,[ones(length(speedN_boot),1) speedN_boot]);
    BN_high(s,:) = regress(highN_boot,[ones(length(speedN_boot),1) speedN_boot]);
end

LN = mean(BN_low(:,2)); HN = mean(BN_high(:,2));
LF = mean(BF_low(:,2)); HF = mean(BF_high(:,2));
errLN = prctile(BN_low(:,2),[2.5 97.5]) -LN;
errLF = prctile(BF_low(:,2),[2.5 97.5]) -LF;
errHN = prctile(BN_high(:,2),[2.5 97.5]) -HN;
errHF = prctile(BF_high(:,2),[2.5 97.5]) -HF;

figure
bar([1 2 3 4],[LN LF HN HF])
hold on
errorbar2([1 2 3 4],[LN LF HN HF],[errLN' errLF' errHN' errHF'],0.5,'k')
set(gca,'xtick',[1.5 3.5],'xticklabel',{'Slow gamma' 'Fast gamma'},'ylim',[-0.3 0.5],'ytick',[-0.3:0.1:0.5])
ylabel('Slope speed vs.gamma power')

%Save Slow Gamma Figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_withinday_gamma_slope.pdf', m, d, y);
print('-dpdf', savestring)

%% Plot slopes for each animal
bin = [1/4 1 4 16];
for an =1:length(eegA)
    lowN = []; lowF = [];
    highN = []; highF = [];
    speedN = []; speedF = [];

	tmp_lowN = []; tmp_highN = []; tmp_lowF = []; tmp_highF = [];
    speedF = [speedF; eegA(an).output{1}(1).speed];
    speedN = [speedN; eegB(an).output{1}(1).speed];

    for tet = 1:length(eegA(an).output{1})
        tmp_lowF = [tmp_lowF eegA(an).output{1}(tet).lowgamma_power];
        tmp_highF = [tmp_highF eegA(an).output{1}(tet).highgamma_power];
    end
    tmp_lowF = mean(tmp_lowF,2);
    tmp_highF = mean(tmp_highF,2);
    for tet = 1:length(eegB(an).output{1})
        tmp_lowN = [tmp_lowN eegB(an).output{1}(tet).lowgamma_power];
        tmp_highN = [tmp_highN eegB(an).output{1}(tet).highgamma_power];
    end
    tmp_lowN = mean(tmp_lowN,2);
    tmp_highN = mean(tmp_highN,2);
  
    lowF = [lowF; (tmp_lowF - mean(tmp_lowF))./std(tmp_lowF)];
    highF = [highF; (tmp_highF - mean(tmp_highF))./std(tmp_highF)];
          
    lowN = [lowN; (tmp_lowN - mean(tmp_lowF))./std(tmp_lowF)];
    highN = [highN; (tmp_highN - mean(tmp_highF))./std(tmp_highF)];

    invalid = speedN<1/8 | lowN>12 | highN >12;
    speedN(invalid)=[]; highN(invalid)=[]; lowN(invalid)=[];
    speedN = log(speedN);

    invalid = speedF<1/8 | lowF>10 | highF >10;
    speedF(invalid)=[]; highF(invalid)=[]; lowF(invalid)=[];
    speedF = log(speedF);

    %Recenter data at origin to isolate changes in slope
    speedF = speedF - mean(speedF);
    lowF = lowF - mean(lowF);
    highF = highF - mean(highF);
    
    speedN = speedN - mean(speedN);
    lowN = lowN - mean(lowN);
    highN = highN - mean(highN);
    
    [BF_low Flowint] = regress(lowF,[ones(length(speedF),1) speedF])
        [BN_low Nlowint] = regress(lowN,[ones(length(speedN),1) speedN])

        [BF_high Fhighint] = regress(highF,[ones(length(speedF),1) speedF])
       
    [BN_high Nhighint] = regress(highN,[ones(length(speedN),1) speedN])
    
    
%     figure(1)
%     hold on
%     plot([1 2],[BN_low(2) BF_low(2)],'k', [1 2], [BN_low(2) BF_low(2)],'ko','MarkerFace','k')
%     plot([3 4],[BN_high(2) BF_high(2)],'r', [3 4], [BN_high(2) BF_high(2)],'ro','MarkerFace','r')
%     set(gca,'xtick',[1 2 3 4],'xticklabel','Novel|Familiar|Novel|Familiar','xlim',[0.8 4.2],'FontSize',18)
%     set(gca,'yLim',[-0.75 0.75],'ytick',-0.75:0.25:0.75,'FontSize',18)
%     ylabel('Slope','FontSize',18)
%     title({'Within day comparison of slope' , 'speed vs. gamma power'})
%     box off
end
