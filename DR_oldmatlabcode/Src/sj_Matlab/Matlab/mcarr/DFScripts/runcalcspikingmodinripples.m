%% RUN FILTER FOR ALL CELLS
%Animal selection
animals = {'Bond','Frank','Ten'};

%Define epochs
epochfilter = [];
epochfilter{1} = 'isequal($type, ''run'')';
%epochfilter{1} = 'isequal($type, ''sleep'')';

cellfilter = '($meanrate<7)';

%Define iterator
iterator = 'multicellanal';

timefilter = {{'get2dstate', '$velocity<4'}};
%timefilter = {{'get2dstate', '$velocity<4 & $immobilitytime>60'}};

%Define iterator
iterator = 'multicellanal';

f3 = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);
f3 = setfilterfunction(f3, 'calcspikingmodinripples', {'spikes','ripples','cellinfo'},'min_cells',5,'minthresh',3);
f3 = runfilter(f3);

save('/data13/mcarr/RipplePaper/gammaspikingmodulation.mat','f3')
save('/data13/mcarr/RipplePaper/gammaspikingmodulation_sleep.mat','f3')
%% Plot gamma modulation of CA3 and CA1 during ripples, global gamma phase
bin = -pi:pi/4:pi;
gam1 = []; gam1_first = []; gam3 = []; gam3_first = [];
gam1_preceding = []; gam3_preceding = [];
for an = 1:length(f3)
    for d = 1:length(f3(an).output)
        for e = 1:length(f3(an).output{d})
            if ~isempty(f3(an).output{d}(e).celldata_first)
                gam3 = [gam3; f3(an).output{d}(e).celldata(f3(an).output{d}(e).celldata(:,6)==3,4)];
                gam3_preceding = [gam3_preceding; f3(an).output{d}(e).preceding(f3(an).output{d}(e).preceding(:,6)==3,4)];
                gam3_first =[gam3_first; f3(an).output{d}(e).celldata_first(f3(an).output{d}(e).celldata_first(:,6)==3,4)];
                gam1 = [gam1; f3(an).output{d}(e).celldata(f3(an).output{d}(e).celldata(:,6)==1,4)];
                gam1_preceding = [gam1_preceding; f3(an).output{d}(e).preceding(f3(an).output{d}(e).preceding(:,6)==1,4)];
                gam1_first =[gam1_first; f3(an).output{d}(e).celldata_first(f3(an).output{d}(e).celldata_first(:,6)==1,4)];
            end
            
        end
    end
end
gam1(isnan(gam1)) =[]; gam3(isnan(gam3))=[]; gam1_first(isnan(gam1_first))=[]; gam3_first(isnan(gam3_first))=[];
gam1_preceding(isnan(gam1_preceding)) = []; gam3_preceding(isnan(gam3_preceding)) = [];
gam1(gam1<-pi) = -pi; gam1(gam1>pi) = pi; gam1_first(gam1_first<-pi) = -pi; gam1_first(gam1_first>pi) = pi;
gam3(gam3<-pi) = -pi; gam3(gam3>pi) = pi; gam3_first(gam3_first<-pi) = -pi; gam3_first(gam3_first>pi) = pi;
gam1_preceding(gam1_preceding<-pi) = -pi; gam1_preceding(gam1_preceding>pi) = pi; gam3_preceding(gam3_preceding<-pi) = -pi; gam3_preceding(gam3_preceding>pi) = pi;

h1 = histc(gam1,bin); h3= histc(gam3,bin); h1_first = histc(gam1_first,bin); h3_first = histc(gam3_first,bin);
h1_preceding = histc(gam1_preceding,bin); h3_preceding = histc(gam3_preceding,bin);

h1 = [h1(1:end-1)./length(gam1); h1(1:end-1)./length(gam1); h1(1)./length(gam1)];
h3 = [h3(1:end-1)./length(gam3); h3(1:end-1)./length(gam3); h3(1)./length(gam3)];
h1_preceding = [h1_preceding(1:end-1)./length(gam1_preceding); h1_preceding(1:end-1)./length(gam1_preceding); h1_preceding(1)./length(gam1_preceding)];
h3_preceding = [h3_preceding(1:end-1)./length(gam3_preceding); h3_preceding(1:end-1)./length(gam3_preceding); h3_preceding(1)./length(gam3_preceding)];
h1_first = [h1_first(1:end-1)./length(gam1_first); h1_first(1:end-1)./length(gam1_first); h1_first(1)./length(gam1_first)];
h3_first = [h3_first(1:end-1)./length(gam3_first); h3_first(1:end-1)./length(gam3_first); h3_first(1)./length(gam3_first)];

x = [bin(1:end-1) 2*pi+bin(1:end)];

figure
plot(x,h1,'b',x,h3,'r')
set(gca,'xlim',x([1 end]),'ylim',[0.1 0.15],'xtick',x(1:2:end),'xticklabel',x(1:2:end)*180/pi)
legend([{'CA1'},{'CA3'}])
xlabel('Gamma Phase')
ylabel('Proportion of spikes')

%Rayleigh test on gam3 and gam1 to test for significant modulation:

%RUN:
%   gam3: theta = 0.24, p = 0.001
%   gam1: theta = 1.96, p = 1e-10

%SLEEP:
%   gam3: theta = -0.49, p < 0.01
%   gam1: theta = 2.36,  p <0.01

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_4a.pdf', m, d, y);
print('-dpdf', savestring)

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_7i.pdf', m, d, y);
print('-dpdf', savestring)


figure
plot(x,h1_first,'c',x,h3_first,'m')
set(gca,'xlim',x([1 end]),'ylim',[0.1 0.15],'xtick',x(1:2:end),'xticklabel',x(1:2:end)*180/pi)
legend([{'CA1'},{'CA3'}])
xlabel('Gamma Phase')
ylabel('Proportion of spikes')

%Rayleigh test on gam3_first and gam1_first to test for significant modulation:

%RUN:
%   gam3_first: theta = -0.06, p = 0.001
%   gam1_first: theta = 0.94, p = 0.001

%SLEEP:
%   gam3_first: theta = -1.7, p = 0.5
%   gam1_first: theta = 2,  p = 0.05

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_4d.pdf', m, d, y);
print('-dpdf', savestring)

%Run permutation test to determine whether CA3 comes before CA1
nperm = 1000;
q = [anglemean(gam3)-anglemean(gam1) anglemean(gam3_first)-anglemean(gam1_first)];
Q = nan(length(nperm),2);
Z = [gam1; gam3]; N1 = length(gam1);
Z1 = [gam1_first; gam3_first]; N11 = length(gam1_first);
for s = 1:nperm
    ind = randperm(length(Z));
    x = Z(ind(1:N1)); y = Z(ind(N1+1:end));
    Q(s,1) = anglemean(x)-anglemean(y);
    ind = randperm(length(Z1));
    x = Z1(ind(1:N11)); y = Z1(ind(N11+1:end));
    Q(s,2) = anglemean(x)-anglemean(y);
end
pall = sum(abs(Q(:,1))>abs(q(1)))./size(Q,1); pfirst = sum(abs(Q(:,2))>abs(q(2)))./size(Q,1);
%RUN: pall<1e-5 pfirst<0.01
%SLEEP: p>0.2

%% COMPARE DEPTH OF MODULATION BETWEEN SWRS and PRECEDING
q = [(max(h1)-min(h1))./(max(h1)+min(h1)) (max(h3)-min(h3))./(max(h3)+min(h3))...
    (max(h1_preceding)-min(h1_preceding))./(max(h1_preceding)+min(h1_preceding)) (max(h3_preceding)-min(h3_preceding))./(max(h3_preceding)+min(h3_preceding)) ];

nboot = 1000;
Q = nan(length(nboot),4);
for s = 1:nboot
    boot1 = gam1(ceil(length(gam1)*rand(length(gam1),1)));
    boot3 = gam3(ceil(length(gam3)*rand(length(gam3),1)));
    q1 = histc(boot1,bin);  q3 = histc(boot3,bin);
    q1 = [q1(1:end-1)./length(gam1); q1(1:end-1)./length(gam1); q1(1)./length(gam1)];
    q3 = [q3(1:end-1)./length(gam3); q3(1:end-1)./length(gam3); q3(1)./length(gam3)];
    
    Q(s,[1 2]) = [(max(q1)-min(q1))./(max(q1)+min(q1)) (max(q3)-min(q3))./(max(q3)+min(q3))];
    
    boot1 = gam1_preceding(ceil(length(gam1_preceding)*rand(length(gam1_preceding),1)));
    boot3 = gam3_preceding(ceil(length(gam3_preceding)*rand(length(gam3_preceding),1)));
    q1 = histc(boot1,bin);  q3 = histc(boot3,bin);
    q1 = [q1(1:end-1)./length(gam1_preceding); q1(1:end-1)./length(gam1_preceding); q1(1)./length(gam1_preceding)];
    q3 = [q3(1:end-1)./length(gam3_preceding); q3(1:end-1)./length(gam3_preceding); q3(1)./length(gam3_preceding)];
    
    Q(s,[3 4]) = [(max(q1)-min(q1))./(max(q1)+min(q1)) (max(q3)-min(q3))./(max(q3)+min(q3))];
    
end

%CA1:   sum(Q(:,1)<mean(Q(:,3)))./nboot
%RUN:   p<1e-5
%SLEEP: p<0.02

%CA3:   sum(Q(:,2)<mean(Q(:,4)))./nboot
%RUN:   p>0.5
%IMMOBILE: p>0.4

%Rayleigh test on gam3_preceding and gam1_preceding to test for significant modulation:

%RUN:
%   gam3_preceding: theta = -0.66, p = 1e-5
%   gam1_preceding: theta = 0.66, p = 0.01

%RUN CA1 depth of modulation > sleep depth of modulation p<0.05
%depth of modulation no different between RUN & immobile for CA3 swrs or 
%CA1 or CA3 spikes preceding SWRs

figure
bar([1 2 4 5],mean(Q(:,[4 2 3 1])),1,'b')
hold on
errorbar2([1 2 4 5],mean(Q(:,[4 2 3 1])),std(Q(:,[4 2 3 1])),'k')
set(gca,'xtick',[1 2 4 5],'xticklabel',[{'CA3 preceding'},{'CA3 SWRs'},{'CA1 preceding'},{'CA1 SWR'}])
ylabel('Depth of modulation')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_4b.pdf', m, d, y);
print('-dpdf', savestring)

%% COMPARE DEPTH OF MODULATION BETWEEN ALL SPIKES AND FIRST SPIKES
q = [(max(h1)-min(h1))./(max(h1)+min(h1)) (max(h3)-min(h3))./(max(h3)+min(h3))...
    (max(h1_first)-min(h1_first))./(max(h1_first)+min(h1_first)) (max(h3_first)-min(h3_first))./(max(h3_first)+min(h3_first))];

nboot = 5000;
Q = nan(length(nboot),4);
for s = 1:nboot
    boot1 = gam1(ceil(length(gam1)*rand(length(gam1),1)));
    boot3 = gam3(ceil(length(gam3)*rand(length(gam3),1)));
    q1 = histc(boot1,bin);  q3 = histc(boot3,bin);
    q1 = [q1(1:end-1)./length(gam1); q1(1:end-1)./length(gam1); q1(1)./length(gam1)];
    q3 = [q3(1:end-1)./length(gam3); q3(1:end-1)./length(gam3); q3(1)./length(gam3)];
    
    Q(s,[1 2]) = [(max(q1)-min(q1))./(max(q1)+min(q1)) (max(q3)-min(q3))./(max(q3)+min(q3))];
    
    boot1 = gam1_first(ceil(length(gam1_first)*rand(length(gam1_first),1)));
    boot3 = gam3_first(ceil(length(gam3_first)*rand(length(gam3_first),1)));
    q1 = histc(boot1,bin);  q3 = histc(boot3,bin);
    q1 = [q1(1:end-1)./length(gam1_first); q1(1:end-1)./length(gam1_first); q1(1)./length(gam1_first)];
    q3 = [q3(1:end-1)./length(gam3_first); q3(1:end-1)./length(gam3_first); q3(1)./length(gam3_first)];
    
    Q(s,[3 4]) = [(max(q1)-min(q1))./(max(q1)+min(q1)) (max(q3)-min(q3))./(max(q3)+min(q3))];
    
end


%CA1:   sum(Q(:,3)<mean(Q(:,1)))./nboot
%RUN:   p>0.15

%CA3:   sum(Q(:,4)<mean(Q(:,2)))./nboot
%RUN:   p<0.05

figure
bar([1 2 4 5],mean(Q(:,[4 2 3 1])),1,'b')
hold on
errorbar2([1 2 4 5],mean(Q(:,[4 2 3 1])),std(Q(:,[4 2 3 1])),'k')
set(gca,'xtick',[1 2 4 5],'xticklabel',[{'CA3 first'},{'CA3 all'},{'CA1 first'},{'CA1 all'}])
ylabel('Depth of modulation')

%% Run permutation test to determine whether there is a difference in preferred phase between:
% preceding and all spikes, preceding and first spikes, first spikes and all spikes

nperm = 1000;
q1 = [anglemean(gam1_preceding)-anglemean(gam1) anglemean(gam1_preceding)-anglemean(gam1_first) anglemean(gam1_first)-anglemean(gam1)];
Q1 = nan(length(nperm),3);
Z = [gam1_preceding; gam1]; N1 = length(gam1);
Z2 = [gam1_preceding; gam1_first]; N2 = length(gam1_first);
Z3 = [gam1_first; gam1]; N3 = length(gam1_first);

for s = 1:nperm
    ind = randperm(length(Z));
    x = Z(ind(1:N1)); y = Z(ind(N1+1:end));
    Q1(s,1) = anglemean(x)-anglemean(y);
    ind = randperm(length(Z2));
    x = Z2(ind(1:N2)); y = Z2(ind(N2+1:end));
    Q1(s,2) = anglemean(x)-anglemean(y);
    ind = randperm(length(Z3));
    x = Z3(ind(1:N3)); y = Z3(ind(N3+1:end));
    Q1(s,3) = anglemean(x)-anglemean(y);
end
p1_preceding_all = sum(abs(Q1(:,1))>abs(q1(1)))./nperm
% RUN: p< 1e-5 SLEEP: p > 0.1 
p1_preceding_first = sum(abs(Q1(:,2))>abs(q1(2)))./nperm
% RUN: p>0.5 SLEEP: p > 0.1 
p1_first_all = sum(abs(Q1(:,2))>abs(q1(3)))./nperm
% RUN: p<0.01 SLEEP: p > 0.5 

q3 = [anglemean(gam3_preceding)-anglemean(gam3) anglemean(gam3_preceding)-anglemean(gam3_first) anglemean(gam3_first)-anglemean(gam3)];
Q3 = nan(length(nperm),3);
Z = [gam3_preceding; gam3]; N1 = length(gam3);
Z2 = [gam3_preceding; gam3_first]; N2 = length(gam3_first);
Z3 = [gam3_first; gam3]; N3 = length(gam3_first);

for s = 1:nperm
    ind = randperm(length(Z));
    x = Z(ind(1:N1)); y = Z(ind(N1+1:end));
    Q3(s,1) = anglemean(x)-anglemean(y);
    ind = randperm(length(Z2));
    x = Z2(ind(1:N2)); y = Z2(ind(N2+1:end));
    Q3(s,2) = anglemean(x)-anglemean(y);
    ind = randperm(length(Z3));
    x = Z3(ind(1:N3)); y = Z3(ind(N3+1:end));
    Q3(s,3) = anglemean(x)-anglemean(y);
end
p3_preceding_all = sum(abs(Q3(:,1))>abs(q3(1)))./nperm
% RUN: p<0.01 SLEEP: p<0.05
p3_preceding_first = sum(abs(Q3(:,2))>abs(q3(2)))./nperm
% RUN: p>0.05 SLEEP: p>0.5
p3_first_all = sum(abs(Q3(:,2))>abs(q3(3)))./nperm
% RUN: p>0.3 SLEEP: p<0.05

%% Plot gamma modulation of CA3 and CA1 during ripples, local gamma phase
bin = -pi:pi/4:pi;
gam1 = []; gam1_first = []; gam3 = []; gam3_first = [];
gam1_preceding = []; gam3_preceding = [];
for an = 1:length(f3)
    for d = 1:length(f3(an).output)
        for e = 1:length(f3(an).output{d})
            if ~isempty(f3(an).output{d}(e).celldata_first)
                gam3 = [gam3; f3(an).output{d}(e).celldata(f3(an).output{d}(e).celldata(:,6)==3,3)];
                gam3_preceding = [gam3_preceding; f3(an).output{d}(e).preceding(f3(an).output{d}(e).preceding(:,6)==3,3)];
                gam3_first =[gam3_first; f3(an).output{d}(e).celldata_first(f3(an).output{d}(e).celldata_first(:,6)==3,3)];
                gam1 = [gam1; f3(an).output{d}(e).celldata(f3(an).output{d}(e).celldata(:,6)==1,3)];
                gam1_preceding = [gam1_preceding; f3(an).output{d}(e).preceding(f3(an).output{d}(e).preceding(:,6)==1,3)];
                gam1_first =[gam1_first; f3(an).output{d}(e).celldata_first(f3(an).output{d}(e).celldata_first(:,6)==1,3)];
            end
            
        end
    end
end
gam1(isnan(gam1)) =[]; gam3(isnan(gam3))=[]; gam1_first(isnan(gam1_first))=[]; gam3_first(isnan(gam3_first))=[];
gam1_preceding(isnan(gam1_preceding)) = []; gam3_preceding(isnan(gam3_preceding)) = [];
gam1(gam1<-pi) = -pi; gam1(gam1>pi) = pi; gam1_first(gam1_first<-pi) = -pi; gam1_first(gam1_first>pi) = pi;
gam3(gam3<-pi) = -pi; gam3(gam3>pi) = pi; gam3_first(gam3_first<-pi) = -pi; gam3_first(gam3_first>pi) = pi;
gam1_preceding(gam1_preceding<-pi) = -pi; gam1_preceding(gam1_preceding>pi) = pi; gam3_preceding(gam3_preceding<-pi) = -pi; gam3_preceding(gam3_preceding>pi) = pi;

h1 = histc(gam1,bin); h3= histc(gam3,bin); h1_first = histc(gam1_first,bin); h3_first = histc(gam3_first,bin);
h1_preceding = histc(gam1_preceding,bin); h3_preceding = histc(gam3_preceding,bin);

h1 = [h1(1:end-1)./length(gam1); h1(1:end-1)./length(gam1); h1(1)./length(gam1)];
h3 = [h3(1:end-1)./length(gam3); h3(1:end-1)./length(gam3); h3(1)./length(gam3)];
h1_preceding = [h1_preceding(1:end-1)./length(gam1_preceding); h1_preceding(1:end-1)./length(gam1_preceding); h1_preceding(1)./length(gam1_preceding)];
h3_preceding = [h3_preceding(1:end-1)./length(gam3_preceding); h3_preceding(1:end-1)./length(gam3_preceding); h3_preceding(1)./length(gam3_preceding)];
h1_first = [h1_first(1:end-1)./length(gam1_first); h1_first(1:end-1)./length(gam1_first); h1_first(1)./length(gam1_first)];
h3_first = [h3_first(1:end-1)./length(gam3_first); h3_first(1:end-1)./length(gam3_first); h3_first(1)./length(gam3_first)];
x = [bin(1:end-1) 2*pi+bin(1:end)];

figure
plot(x,h1,'b',x,h3,'r')
set(gca,'xlim',x([1 end]),'ylim',[0.11 0.14],'xtick',x(1:2:end),'xticklabel',x(1:2:end)*180/pi)
legend([{'CA1'},{'CA3'}])
xlabel('Gamma Phase')
ylabel('Proportion of spikes')

figure
plot(x ,h1_first,'c',x,h3_first,'m')
set(gca,'xlim',x([1 end]),'ylim',[0.11 0.14],'xtick',x(1:2:end),'xticklabel',x(1:2:end)*180/pi)
legend([{'CA1'},{'CA3'}])
xlabel('Gamma Phase')
ylabel('Proportion of spikes')

%% HOW MANY CELLS?
epoch = cell(length(f3),1);
for an = 1:length(f3)
    
    for d = 1:length(f3(an).output)
        for e = 1:length(f3(an).output{d})
            if ~isempty(f3(an).output{d}(e).celldata_first)
                epoch{an}(e) = f3(an).epochs{d}(e,1);
            end
            
        end
    end
    epoch{an} = diff([0 epoch{an}])>0;
end

gam1 = []; gam3 = []; 
for an = 1:length(f3)
    for d = 1:length(f3(an).output)
        for e = 1:length(f3(an).output{d})
            if ~isempty(f3(an).output{d}(e).celldata) && epoch{an}(e)
                gam3 = [gam3; unique(f3(an).output{d}(e).celldata(f3(an).output{d}(e).celldata(:,6)==3,5))];
                gam1 = [gam1; unique(f3(an).output{d}(e).celldata(f3(an).output{d}(e).celldata(:,6)==1,5))];
            end
            
        end
    end
end


%% COMPARE PHASE OF ALL SPIKES DURING WAKING AND SLEEPING
agam1 = gam1; agam3 = gam3;
qgam1 = gam1; qgam3 = gam3;

nboot = 1000;
Q = nan(length(nboot),4);
for s = 1:nboot
    boot1 = agam1(ceil(length(agam1)*rand(length(agam1),1)));
    boot3 = agam3(ceil(length(agam3)*rand(length(agam3),1)));
    Q(s,1) = anglemean(boot1);  Q(s,2) = anglemean(boot3);

    boot1 = qgam1(ceil(length(qgam1)*rand(length(qgam1),1)));
    boot3 = qgam3(ceil(length(qgam3)*rand(length(qgam3),1)));
    Q(s,3) = anglemean(boot1);  Q(s,4) = anglemean(boot3);
end
mean(Q)
prctile(Q,[2.5 97.5])
