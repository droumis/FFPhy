%% RUN FILTER FOR ALL CELLS
%Animal selection
animals = {'Bond','Frank','Ten'};

%Define epochs
epochfilter = [];
epochfilter{1} = 'isequal($type, ''run'')';

cellfilter = '($meanrate<7) & ($numspikes > 100)';

%Define iterator
iterator = 'multicellanal';

timefilter = {{'get2dstate', '($velocity<4)'}};

%Define iterator
iterator = 'multicellanal';

%create training data by calulating the linearized rates of all cells
f3 = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);
f3 = setfilterfunction(f3, 'calcspikingmodinripples', {'spikes','ripples','cellinfo'},'min_cells',5);
f3 = runfilter(f3);


%% Plot gamma modulation of CA3 and CA1 during ripples, global gamma phase
bin = -pi:pi/4:pi;
gam1 = []; gam1_first = []; gam3 = []; gam3_first = [];
for an = 1:length(f3)
    for d =1:length(f3(an).output)
        for e = 1:length(f3(an).output{d})
            if ~isempty(f3(an).output{d}(e).celldata_first)
                gam3 = [gam3; f3(an).output{d}(e).celldata(f3(an).output{d}(e).celldata(:,6)==3,4)];
                gam3_first =[gam3_first; f3(an).output{d}(e).celldata_first(f3(an).output{d}(e).celldata_first(:,6)==3,4)];
                gam1 = [gam1; f3(an).output{d}(e).celldata(f3(an).output{d}(e).celldata(:,6)==1,4)];
                gam1_first =[gam1_first; f3(an).output{d}(e).celldata_first(f3(an).output{d}(e).celldata_first(:,6)==1,4)];

            end
            
        end
    end
end
gam1(isnan(gam1)) =[]; gam3(isnan(gam3))=[]; gam1_first(isnan(gam1_first))=[]; gam3_first(isnan(gam3_first))=[];
gam1(gam1<-pi) = -pi; gam1(gam1>pi) = pi; gam1_first(gam1_first<-pi) = -pi; gam1_first(gam1_first>pi) = pi;
gam3(gam3<-pi) = -pi; gam3(gam3>pi) = pi; gam3_first(gam3_first<-pi) = -pi; gam3_first(gam3_first>pi) = pi;
h1 = histc(gam1,bin); h3= histc(gam3,bin); h1_first = histc(gam1_first,bin); h3_first = histc(gam3_first,bin);

h1 = [h1(1:end-1)./length(gam1); h1(1:end-1)./length(gam1); h1(1)./length(gam1)];
h3 = [h3(1:end-1)./length(gam3); h3(1:end-1)./length(gam3); h3(1)./length(gam3)];
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
%   gam3: theta = 0.4371, p = 0.0114
%   gam1: theta = 1.9834, p = 1e-5

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_probability_spikes_gamma_phase.pdf', m, d, y);
print('-dpdf', savestring)


figure
plot(x,h1_first,'c',x,h3_first,'m')
set(gca,'xlim',x([1 end]),'ylim',[0.1 0.15],'xtick',x(1:2:end),'xticklabel',x(1:2:end)*180/pi)
legend([{'CA1'},{'CA3'}])
xlabel('Gamma Phase')
ylabel('Proportion of spikes')

%Rayleigh test on gam3_first and gam1_first to test for significant modulation:
%   gam3_first: theta = -0.2519, p = 0.002
%   gam1_first: theta = 0.9, p = 0.005

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_proportion_first_spike_gamma_phase.pdf', m, d, y);
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
%on average, CA3 comes before CA1
%   all spikes: p = 1e-10
%   first spikes: p = 0.02

%% Plot preferred phase of all neurons
bin = -pi:pi/4:pi;
ca1 = nan(1000,4);
ca3 = nan(1000,4);
count = 1;
for an = 1:length(f3)
    for d = 1:length(f3(an).output)
        for e = 1:length(f3(an).output{d})
            if ~isempty(f3(an).output{d}(e).celldata)
                for c = unique(f3(an).output{d}(e).celldata(:,5))'
                    [m r] = anglemean(f3(an).output{d}(e).celldata(f3(an).output{d}(e).celldata(:,5)==c,3));
                    [m1 r1] = anglemean(f3(an).output{d}(e).celldata_first(f3(an).output{d}(e).celldata_first(:,5)==c,3));
                    if f3(an).output{d}(e).celldata(f3(an).output{d}(e).celldata(:,5)==c,6)==1
                    	ca1(count,:) = [m r m1 r1];
                    elseif f3(an).output{d}(e).celldata(f3(an).output{d}(e).celldata(:,5)==c,6)==3
                        ca3(count,:) =[m r m1 r1];
                    end
                    count = count+1;   
                end 
            end
        end
    end
end


invalid = isnan(ca1(:,1)) | isnan(ca1(:,2));
ca1(invalid,:) = [];
invalid = isnan(ca3(:,1)) | isnan(ca3(:,2));
ca3(invalid,:) = [];

x = [bin(1:end-1) 2*pi+bin(1:end)];
h1 = histc(ca1(:,1),bin); h3= histc(ca3(:,1),bin); h1_first = histc(ca1(:,3),bin); h3_first = histc(ca3(:,3),bin);

h1 = [h1(1:end-1)./size(ca1,1); h1(1:end-1)./size(ca1,1); h1(1)./size(ca1,1)];
h3 = [h3(1:end-1)./size(ca3,1); h3(1:end-1)./size(ca3,1); h3(1)./size(ca3,1)];
h1_first = [h1_first(1:end-1)./size(ca1,1); h1_first(1:end-1)./size(ca1,1); h1_first(1)./size(ca1,1)];
h3_first = [h3_first(1:end-1)./size(ca3,1); h3_first(1:end-1)./size(ca3,1); h3_first(1)./size(ca3,1)];

figure
plot(x,h1,'b',x,h3,'r')
set(gca,'xlim',x([1 end]),'ylim',[0.05 0.2],'xtick',x(1:2:end),'xticklabel',x(1:2:end)*180/pi)
legend([{'CA1'},{'CA3'}])
xlabel('Prefered Gammma Phase')
ylabel('Proportion of neurons')

%Rayleigh test on ca3(:,1) and ca1(:,1) to test for significant modulation:
%   ca3: theta = -0.626, p = 0.01
%   ca1: theta = 1.316, p = 0.03

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_neurons_perferred_gamma_phase.pdf', m, d, y);
print('-dpdf', savestring)

figure
plot(x ,h1_first,'c',x,h3_first,'m')
set(gca,'xlim',x([1 end]),'ylim',[0.05 0.2],'xtick',x(1:2:end),'xticklabel',x(1:2:end)*180/pi)
legend([{'CA1'},{'CA3'}])
xlabel('Prefered Gammma Phase')
ylabel('Proportion of neurons')

%Rayleigh test on ca3(:,3) and ca1(:,3) to test for significant modulation:
%   ca3: theta = -1.1077, p = 0.005
%   ca1: theta = -0.3462, p = 0.1

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_neurons_preferred_phase_firstspike.pdf', m, d, y);
print('-dpdf', savestring)


%Run permutation test to determine whether CA3 comes before CA1
nperm = 10000;
q = [anglemean(ca3(:,1))-anglemean(ca1(:,1)) anglemean(ca3(:,3))-anglemean(ca1(:,3))];
Q = nan(length(nperm),2);
Z = [ca1(:,1); ca3(:,1)]; N1 = size(ca1,1);
Z1 = [ca1(:,3); ca3(:,3)]; N11 = size(ca1,1);
for s = 1:nperm
    ind = randperm(length(Z));
    x = Z(ind(1:N1)); y = Z(ind(N1+1:end));
    Q(s,1) = anglemean(x)-anglemean(y);
    ind = randperm(length(Z1));
    x = Z1(ind(1:N11)); y = Z1(ind(N11+1:end));
    Q(s,2) = anglemean(x)-anglemean(y);
end
pall = sum(abs(Q(:,1))>abs(q(1)))./size(Q,1); pfirst = sum(abs(Q(:,2))>abs(q(2)))./size(Q,1);

% while CA3 has preferred phase before CA1, first spike is not significantly different
%   all spikes: p = 0.03
%   first spikes: p = 0.2

%% Plot gamma modulation of CA3 and CA1 during ripples, local gamma phase
bin = -pi:pi/4:pi;
gam1 = []; gam1_first = []; gam3 = []; gam3_first = [];
for an = 1:length(f3)
    for d =1:length(f3(an).output)
        for e = 1:length(f3(an).output{d})
            if ~isempty(f3(an).output{d}(e).celldata_first)
                gam3 = [gam3; f3(an).output{d}(e).celldata(f3(an).output{d}(e).celldata(:,6)==3,3)];
                all3 =[gam3_first; f3(an).output{d}(e).celldata_first(f3(an).output{d}(e).celldata_first(:,6)==3,3)];
                gam1 = [gam1; f3(an).output{d}(e).celldata(f3(an).output{d}(e).celldata(:,6)==1,3)];
                gam1_first =[gam1_first; f3(an).output{d}(e).celldata_first(f3(an).output{d}(e).celldata_first(:,6)==1,3)];

            end
            
        end
    end
end
gam1(isnan(gam1)) =[]; gam3(isnan(gam3))=[]; gam1_first(isnan(gam1_first))=[]; gam3_first(isnan(gam3_first))=[];
h1 = histc(gam1,bin); h3= histc(gam3,bin); h1_first = histc(gam1_first,bin); h3_first = histc(gam3_first,bin);

h1 = [h1(1:end-1)./length(gam1); h1(1:end-1)./length(gam1); h1(1)./length(gam1)];
h3 = [h3(1:end-1)./length(gam3); h3(1:end-1)./length(gam3); h3(1)./length(gam3)];
h1_first = [h1_first(1:end-1)./length(gam1_first); h1_first(1:end-1)./length(gam1_first); h1_first(1)./length(gam1_first)];
h3_first = [h3_first(1:end-1)./length(gam3_first); h3_first(1:end-1)./length(gam3_first); h3_first(1)./length(gam3_first)];

x = [bin(1:end-1) 2*pi+bin(1:end)];

figure
plot(x,h1,'b',x,h3,'r')
set(gca,'xlim',x([1 end]),'ylim',[0.065 0.105],'xtick',x(1:2:end),'xticklabel',x(1:2:end)*180/pi)
legend([{'CA1'},{'CA3'}])
xlabel('Gamma Phase')
ylabel('Proportion of spikes')

figure
plot(x ,h1_first,'c',x,h3_first,'m')
set(gca,'xlim',x([1 end]),'ylim',[0.065 0.105],'xtick',x(1:2:end),'xticklabel',x(1:2:end)*180/pi)
legend([{'CA1'},{'CA3'}])
xlabel('Gamma Phase')
ylabel('Proportion of spikes')
