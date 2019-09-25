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

f = createfilter('animal',animals,'epochs',epochfilter,'cells',cellfilter,'excludetime', timefilter,'iterator',iterator);
f = setfilterfunction(f, 'calcripplespikingmodulation', {'spikes','ripples','cellinfo'},'min_cells',0,'minthresh',3);
f = runfilter(f);

%% PLOT MODULATION BY LOCAL RIPPLE OSCILLATION

bin = -pi:pi/8:pi;
rip1 = []; rip3 = [];
thresh1 = []; thresh3 = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if ~isempty(f(an).output{d}(e).celldata)
                rip1 = [rip1; f(an).output{d}(e).celldata(f(an).output{d}(e).celldata(:,7)==1,5)];
                thresh1 = [thresh1; f(an).output{d}(e).celldata(f(an).output{d}(e).celldata(:,7)==1,3)];
                rip3 = [rip3; f(an).output{d}(e).celldata(f(an).output{d}(e).celldata(:,7)==3,5)];
                thresh3 = [thresh3; f(an).output{d}(e).celldata(f(an).output{d}(e).celldata(:,7)==3,3)];
            end
        end
    end
end

invalid = isnan(rip1) | rip1<-pi | rip1>pi; rip1(invalid) = []; thresh1(invalid) = [];
invalid = isnan(rip3) | rip3<-pi | rip3>pi; rip3(invalid) = []; thresh3(invalid) = [];

threshbin =[0 3.5; 0 5; 0 6; 0 7; 0 100];
x = [bin(1:end-1) 2*pi+bin(1:end)];

R1 = zeros(length(x),size(threshbin,1)); R3 = zeros(size(R1));
count = 1;
for i = threshbin'
    r1 = histc(rip1(thresh1>i(1) & thresh1<i(2)),bin);
    r3= histc(rip3(thresh3>i(1) & thresh3<i(2)),bin);
    r1 = [r1(1:end-1)./sum(thresh1>i(1) & thresh1<i(2)); r1(1:end-1)./sum(thresh1>i(1) & thresh1<i(2)); r1(1)./sum(thresh1>i(1) & thresh1<i(2))];
    r3 = [r3(1:end-1)./sum(thresh3>i(1) & thresh3<i(2)); r3(1:end-1)./sum(thresh3>i(1) & thresh3<i(2)); r3(1)./sum(thresh3>i(1) & thresh3<i(2))];
    R1(:,count) = r1; R3(:,count) = r3;
    count = count+1;
end
    
figure
plot(x,R1(:,1),'r', x,R1(:,2),'m',x,R1(:,3),'y',x,R1(:,4),'g',x,R1(:,5),'b')
set(gca,'xlim',x([1 end]),'ylim',[0 0.2],'xtick',x(1:4:end),'xticklabel',(x(1:4:end))*180/pi)
legend([{'CA1 < 4std'},{'CA1<5std'},{'CA1<6std'},{'CA1<7std'},{'CA1<20std'}])
xlabel('Ripple Phase')
ylabel('Proportion of spikes')


% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_CA1_ripplemodulatedspiking.pdf', m, d, y);
print('-dpdf', savestring)


figure
plot(x,R3(:,1),'r', x,R3(:,2),'m',x,R3(:,3),'y',x,R3(:,4),'g',x,R3(:,5),'b')
set(gca,'xlim',x([1 end]),'ylim',[0 0.2],'xtick',x(1:4:end),'xticklabel',(x(1:4:end))*180/pi)
legend([{'CA3 < 4std'},{'CA3<5std'},{'CA3<6std'},{'CA3<7std'},{'CA3<20std'}])
xlabel('Ripple Phase')
ylabel('Proportion of spikes')


% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_CA3_ripplemodulatedspiking.pdf', m, d, y);
print('-dpdf', savestring)



% Use rayleigh test to determine whether spiking relative to local
% oscillation is significant and the preferred phase
%   rip1: theta = 130 degrees, p<1e-5
%   rip3: theta = 130 degrees, p<1e-5

