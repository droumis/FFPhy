%% RUN FILTER FOR CA1 ACROSS DAYS

%animal selection
animals = {'Conley','Corriander','Eight','Five','Miles','Ten'};


% Epoch selection
epochfilter = [];
for i = 1:10
    epochfilter{i} = ['($exposure == ',num2str(i),') & isequal($description,''TrackA'')'];
end

% Tetrode selection
ca1tetfilter =  '(isequal($area, ''CA1'') & $numcells>1 )';

% Time selection
timefilter = {};

%Select iterator
iterator = 'epocheegnonreferenceanal';

% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter, 'eegtetrodes',ca1tetfilter,'iterator', iterator);
f = setfilterfunction(f,'calceegpowerspeed',{'eeg','pos','linpos'});
f = runfilter(f);

power = f;

save('/data13/mcarr/VelocityPaper/power.mat','power')

%% LOAD FILTER
load '/data13/mcarr/VelocityPaper/power.mat'

%% PLOT DATA FOR DAY 1, CONTROL FOR BEHAVIOR

bin = [1/4 1/2 1 2 4 8 16 32];

lowgam = []; highgam = []; speed = [];
day = 1;
count = 100;
for an = 1:length(power)
    if length(power(an).output)>=day && ~isempty(power(an).output{day})
        tmp_lowgam = []; tmp_highgam = [];
        speed = [speed; power(an).output{day}(1).speed];
        for tet = 1:length(power(an).output{day})
            tmp_lowgam = [tmp_lowgam power(an).output{day}(tet).lowgamma_power];
            tmp_highgam = [tmp_highgam power(an).output{day}(tet).highgamma_power];
        end
        tmp_lowgam = mean(tmp_lowgam,2);
        tmp_highgam = mean(tmp_highgam,2);
  
        lowgam = [lowgam; (tmp_lowgam - mean(tmp_lowgam))./std(tmp_lowgam)];
        highgam = [highgam; (tmp_highgam - mean(tmp_highgam))./std(tmp_highgam)];
    end
end

speed(lowgam>10 | highgam>10)=[];
highgam(lowgam>10)=[]; lowgam(lowgam>10)=[];
lowgam(highgam>10)=[]; highgam(highgam>10)=[];
lowgam(speed<1/8) = []; highgam(speed<1/8) = []; speed(speed<1/8) = [];

subs = lookup(speed,bin);

nboot = 1000;
Bboot_lowgam = nan(nboot,2);
Bboot_highgam = nan(nboot,2);

Aboot_lowgam = nan(nboot,length(bin));
Aboot_highgam = nan(nboot,length(bin));

for s = 1:nboot
    speed_boot = []; lowgam_boot = []; highgam_boot = [];
    boot = ceil(length(subs)*rand(2*length(subs),1));
    for q = 1:max(subs)
        tmpboot = find(subs(boot)==q);
        tmpboot = boot(tmpboot);
        tmp = speed(tmpboot);
        tmp = tmp(1:count);
        speed_boot = [speed_boot; tmp];
        tmp = lowgam(tmpboot);
        tmp = tmp(1:count);
        lowgam_boot = [lowgam_boot; tmp];
        tmp = highgam(tmpboot);
        tmp = tmp(1:count);
        highgam_boot = [highgam_boot; tmp];
        clear tmpboot tmp
    end
    Bboot_lowgam(s,:) = regress(lowgam_boot,[ones(length(speed_boot),1) log(speed_boot)]);
    Bboot_highgam(s,:) = regress(highgam_boot,[ones(length(speed_boot),1) log(speed_boot)]);
       
    subs_boot = lookup(speed_boot,bin);
    Aboot_lowgam(s,:) = accumarray(subs_boot,lowgam_boot,[length(bin) 1],@(x) nanmean(x),NaN);
    Aboot_highgam(s,:) = accumarray(subs_boot,highgam_boot,[length(bin) 1],@(x) nanmean(x),NaN);
end


L_lowgam = nan(size(bin)); U_lowgam = nan(size(bin));
L_highgam = nan(size(bin)); U_highgam = nan(size(bin));

for s = 1:length(bin)
    tmp = prctile(Bboot_lowgam(:,2).*log(bin(s)) + Bboot_lowgam(:,1),[2.5 97.5]);
    L_lowgam(s) = tmp(1); U_lowgam(s) = tmp(2);
    tmp = prctile(Bboot_highgam(:,2).*log(bin(s)) + Bboot_highgam(:,1),[2.5 97.5]);
    L_highgam(s) = tmp(1); U_highgam(s) = tmp(2);
end


figure
fill(log([bin fliplr(bin)]), [U_lowgam fliplr(L_lowgam)],'k','FaceAlpha',1,'EdgeColor','none')
hold on
plot(log(bin),log(bin)*mean(Bboot_lowgam(:,2))+mean(Bboot_lowgam(:,1)),'k')
plot(log(bin),mean(Aboot_lowgam),'ko','MarkerFace','k')
errorbar2(log(bin),mean(Aboot_lowgam),std(Aboot_lowgam),0.001,'k','plottype','semilogx')
set(gca,'xtick',log(bin),'xtickLabel',bin)
set(gca,'yLim',[-0.8 0.8],'ytick',-0.8:0.2:.8,'FontSize',20)
xlabel('Speed (cm/sec)','FontSize',22)
ylabel('Low Gamma Power','FontSize',22)
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_speed_lowgamma_day1.pdf', m, d, y);
print('-dpdf', savestring)


figure
fill(log([bin fliplr(bin)]), [U_highgam fliplr(L_highgam)],'k','FaceAlpha',1,'EdgeColor','none')
hold on
plot(log(bin),log(bin)*mean(Bboot_highgam(:,2))+mean(Bboot_highgam(:,1)),'k')    
plot(log(bin),mean(Aboot_highgam),'ko','MarkerFace','k')
errorbar2(log(bin),mean(Aboot_highgam),std(Aboot_highgam),0.001,'k','plottype','semilogx')
set(gca,'xtick',log(bin),'xtickLabel',bin)
set(gca,'yLim',[-1.4 1.4],'ytick',-1.25:0.25:1.25,'FontSize',20)
xlabel('Speed (cm/sec)','FontSize',22)
ylabel('High Gamma Power','FontSize',22)
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_speed_highgamma_day1.pdf', m, d, y);
print('-dpdf', savestring)

%% LOOK AT GAMMA POWER VS. SPEED ACROSS DAYS, CONTROL FOR BEHAVIOR


lowgam = cell(10,1); highgam = cell(10,1); speed = cell(10,1);

for an = 1:length(power)
    for day = 1:length(power(an).output)
        if ~isempty(power(an).output{day})
            tmp_lowgam = []; tmp_highgam = [];
            speed{day} = [speed{day}; power(an).output{day}(1).speed];
            for tet = 1:length(power(an).output{day})
                tmp_lowgam = [tmp_lowgam power(an).output{day}(tet).lowgamma_power];
                tmp_highgam = [tmp_highgam power(an).output{day}(tet).highgamma_power];
            end
            tmp_lowgam = mean(tmp_lowgam,2);
            tmp_highgam = mean(tmp_highgam,2);
  
            lowgam{day} = [lowgam{day}; (tmp_lowgam - mean(tmp_lowgam))./std(tmp_lowgam)];
            highgam{day} = [highgam{day}; (tmp_highgam - mean(tmp_highgam))./std(tmp_highgam)];
        end
    end
end

count = 100;
nboot = 1000;

Bboot_lowgam = nan(nboot,10);
Int_lowgam = nan(nboot,10);
Bboot_highgam = nan(nboot,10);
Int_highgam = nan(nboot,10);

bin = [1/4 1/2 1 2 4 8 16 32];

for day = 1:10
    speed{day}(lowgam{day}>10 | highgam{day}>10)=[];
    highgam{day}(lowgam{day}>10)=[]; lowgam{day}(lowgam{day}>10)=[];
    lowgam{day}(highgam{day}>10)=[]; highgam{day}(highgam{day}>10)=[];
    lowgam{day}(speed{day}<1/8) = []; highgam{day}(speed{day}<1/8) = []; speed{day}(speed{day}<1/8) = [];
    subs = lookup(speed{day},bin);
    
    for s = 1:nboot
        speed_boot = []; lowgam_boot = []; highgam_boot = [];
        boot = ceil(length(subs)*rand(2*length(subs),1));
        for q = 1:max(subs)
            tmpboot = find(subs(boot)==q);
            tmpboot = boot(tmpboot);
            tmp = speed{day}(tmpboot);
            tmp = tmp(1:count);
            speed_boot = [speed_boot; tmp];
            tmp = lowgam{day}(tmpboot);
            tmp = tmp(1:count);
            lowgam_boot = [lowgam_boot; tmp];
            tmp = highgam{day}(tmpboot);
            tmp = tmp(1:count);
            highgam_boot = [highgam_boot; tmp];
            clear tmpboot tmp
        end
        b = regress(lowgam_boot - mean(lowgam_boot), [ones(length(speed_boot),1) log(speed_boot)-mean(log(speed_boot))]);
        Bboot_lowgam(s,day) = b(2);
        Int_lowgam(s,day) = b(1);
        
        b = regress(highgam_boot - mean(highgam_boot), [ones(length(speed_boot),1) log(speed_boot)-mean(log(speed_boot))]);
        Bboot_highgam(s,day) = b(2);
        Int_highgam(s,day) = b(1);
    end

end

figure
hold on
plot(1:10,mean(Bboot_lowgam),'k-o','MarkerFace','K')
errorbar2(1:10,mean(Bboot_lowgam),prctile(Bboot_lowgam-repmat(mean(Bboot_lowgam),nboot,1),[2.5 97.5]),0.001,'k')
set(gca,'xlim',[0 11],'xtick',1:10,'FontSize',18)
set(gca,'ylim',[-0.25 0.1],'ytick',-0.25:0.05:0.1,'FontSize',18)
xlabel('Exposure to Environment','FontSize',24)
ylabel('Slope of Speed vs. Slow Gamma Power','FontSize',24)
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_speed_lowgamma_overdays.pdf', m, d, y);
print('-dpdf', savestring)


figure
hold on
plot(1:10,mean(Bboot_highgam),'k-o','MarkerFace','k')
errorbar2(1:10,mean(Bboot_highgam),prctile(Bboot_highgam-repmat(mean(Bboot_highgam),nboot,1),[2.5 97.5]),0.001,'k')
set(gca,'xlim',[0 11],'xtick',1:10,'FontSize',18)
set(gca,'ylim',[0.1 0.45],'ytick',0:0.05:0.5,'FontSize',18)
xlabel('Exposure to Environment','FontSize',24)
ylabel('Slope of Speed vs. Fast Gamma Power','FontSize',24)
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_speed_highgamma_overdays.pdf', m, d, y);
print('-dpdf', savestring)

%% LOOK AT GAMMA POWER VS. SPEED ACROSS DAYS, INDIVIDUAL ANIMALS

count = 25;
nboot = 500;

Bboot_lowgam = cell(length(power),1);
Bboot_highgam = cell(length(power),1);
bin = [1/4 1 4 16];

for an = 1:length(power)
    Bboot_lowgam{an} = nan(nboot,10);
    Bboot_highgam{an} = nan(nboot,10);
    lowgam = cell(length(power(an).output),1); highgam = cell(size(lowgam)); speed = cell(size(lowgam));
    for day = 1:length(power(an).output)
        if ~isempty(power(an).output{day})
            tmp_lowgam = []; tmp_highgam = [];
            speed{day} = [speed{day}; power(an).output{day}(1).speed];
            for tet = 1:length(power(an).output{day})
                tmp_lowgam = [tmp_lowgam power(an).output{day}(tet).lowgamma_power];
                tmp_highgam = [tmp_highgam power(an).output{day}(tet).highgamma_power];
            end
            tmp_lowgam = mean(tmp_lowgam,2);
            tmp_highgam = mean(tmp_highgam,2);
  
            lowgam{day} = [lowgam{day}; (tmp_lowgam - mean(tmp_lowgam))./std(tmp_lowgam)];
            highgam{day} = [highgam{day}; (tmp_highgam - mean(tmp_highgam))./std(tmp_highgam)];
        end

        speed{day}(lowgam{day}>10 | highgam{day}>10)=[];
        highgam{day}(lowgam{day}>10)=[]; lowgam{day}(lowgam{day}>10)=[];
        lowgam{day}(highgam{day}>10)=[]; highgam{day}(highgam{day}>10)=[];
        lowgam{day}(speed{day}<1/8) = []; highgam{day}(speed{day}<1/8) = []; speed{day}(speed{day}<1/8) = [];
        subs = lookup(speed{day},bin);
    
        for s = 1:nboot
            speed_boot = []; lowgam_boot = []; highgam_boot = [];
            boot = ceil(length(subs)*rand(2*length(subs),1));
            for q = 1:max(subs)
                tmpboot = find(subs(boot)==q);
                tmpboot = boot(tmpboot);
                tmp = speed{day}(tmpboot);
                tmp = tmp(1:count);
                speed_boot = [speed_boot; tmp];
                tmp = lowgam{day}(tmpboot);
                tmp = tmp(1:count);
                lowgam_boot = [lowgam_boot; tmp];
                tmp = highgam{day}(tmpboot);
                tmp = tmp(1:count);
                highgam_boot = [highgam_boot; tmp];
                clear tmpboot tmp
            end
            b = regress(lowgam_boot - mean(lowgam_boot), [ones(length(speed_boot),1) log(speed_boot)-mean(log(speed_boot))]);
            Bboot_lowgam{an}(s,day) = b(2);


            b = regress(highgam_boot - mean(highgam_boot), [ones(length(speed_boot),1) log(speed_boot)-mean(log(speed_boot))]);
            Bboot_highgam{an}(s,day) = b(2);

        end

    end
end
pvalue = nan(length(power),1);
for an = 1:length(power)
    subs = repmat(1:size(Bboot_highgam{an},2),nboot,1);
    subs = reshape(subs,nboot*size(Bboot_highgam{an},2),1);
    bootH = reshape(Bboot_highgam{an},nboot*size(Bboot_highgam{an},2),1);
    bootL = reshape(Bboot_lowgam{an},nboot*size(Bboot_lowgam{an},2),1);
    bH = regress(bootH,[ones(length(subs),1) subs],0.05);
    bL = regress(bootL,[ones(length(subs),1) subs],0.05);
    permH = nan(nboot,1);
    permL = nan(nboot,1);
    for q = 1:nboot
        perm = subs(randperm(length(subs)));
        tmpb = regress(bootH,[ones(length(subs),1) perm],0.05);
        permH(q) = tmpb(2);
        tmpb = regress(bootL,[ones(length(subs),1) perm],0.05);
        permL(q) = tmpb(2);
    end
    pvalue(an,1) = 1-sum(abs(bH(2))>abs(permH))./nboot;
    pvalue(an,2) = 1-sum(abs(bL(2))>abs(permL))./nboot;
end
