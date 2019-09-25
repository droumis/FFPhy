%% RIPPLE PAPER FIGURES: POWER

%Animals
animals = {'Bond','Frank','Ten'};

%Define epochs
epochfilter = [];
epochfilter{1} = 'isequal($type, ''run'')';

%Define iterator
iterator = 'epochbehaveanal';

%Define time filter
timefilter = {{'get2dstate', '($velocity < 4)'}};

f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'iterator', iterator);
f = setfilterfunction(f, 'calcripsstats', {'rips'});
f = runfilter(f);

time = f(1).output{1}(1).time;
frequency = f(1).output{1}(1).frequency;
%% CA1 and CA3 spectrum
ca1 = zeros(length(time),length(frequency)); ca3 = zeros(length(time),length(frequency)); count1 = 0; count3 = 0;
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if ~isempty(f(an).output{d}(e).ca1_spectrum)
                ca1 = ca1 + f(an).output{d}(e).ca1_spectrum;
                count1 = count1+f(an).output{d}(e).nrips;
            end
            if ~isempty(f(an).output{d}(e).ca3_spectrum)
                ca3 = ca3 + f(an).output{d}(e).ca3_spectrum;
                count3 = count3+f(an).output{d}(e).nrips;
            end
        end
    end
end
ca1 = ca1./count1; ca3 = ca3./count3;

%Plot CA1 rip triggered spectrum
figure
surf(time,frequency,ca1')
view(0,90); shading interp; colormap hot
set(gca,'clim',[0 1.5],'xlim',[time(26) time(86)],'ylim',f(1).output{1}(1).frequency([1 end]),'xtick',-0.2:0.1:0.4)
colorbar('ytick',0:0.25:1.5)

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1_riptriggered_spectrum_allruns.png', m, d, y);
print('-dpng', savestring)

%Plot CA3 rip triggered spectrum
figure
surf(time,frequency,ca3')
view(0,90); shading interp;  colormap hot
set(gca,'clim',[0 1.5],'xlim',[time(26) time(86)],'ylim',f(1).output{1}(1).frequency([1 end]),'xtick',-0.2:0.1:0.4)
colorbar('ytick',0:0.25:1.5)

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca3_riptriggered_spectrum_allruns.png', m, d, y);
print('-dpng', savestring)

%% RIPPLE PAPER FIGURES: COHERENCE

%Animals
animals = {'Bond','Frank','Ten'};

%Define epochs
epochfilter = [];
epochfilter{1} = 'isequal($type, ''run'')';

%Define iterator
iterator = 'epochbehaveanal';

%Define time filter
timefilter = {{'get2dstate', '($velocity < 4)'}};

f = createfilter('animal',animals,'epochs',epochfilter,'excludetime',timefilter,'iterator', iterator);
f = setfilterfunction(f, 'calcripcstats', {'ripc'});
f = runfilter(f);

time = f(1).output{1}(1).time;
frequency = f(1).output{1}(1).frequency;
bin = -pi+pi/40:pi/20:pi-pi/40;

%% CA1 and CA3 phase coherence
%NOTE: FIX THIS SECTION!!! MAKE calcripstats compute angular phase.
%Show images of examples
figure
imagesc(time,bin,f(3).output{1}(2).ca1_ca1_phase')
set(gca,'clim',[-.015 0.045],'xlim',[time(1) time(end)],'ylim',bin([1 end]))
colorbar('ytick',-.015:0.015:.045); colormap('hot')
xlabel('Time since ripple detection (s)')
ylabel('Gamma phase')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_example_ca1_ca1_phase.png', m, d, y);
print('-dpng', savestring)

figure
imagesc(time,bin,f(3).output{1}(2).ca1_ca3_phase')
set(gca,'clim',[-.015 0.045],'xlim',[time(1) time(end)],'ylim',bin([1 end]))
colorbar('ytick',-.015:0.015:.045); colormap('hot')
xlabel('Time since ripple detection (s)')
ylabel('Gamma phase')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_example_ca1_ca3_phase.png', m, d, y);
print('-dpng', savestring)

figure
imagesc(time,bin,f(3).output{1}(2).ca3_ca3_phase')
set(gca,'clim',[-.015 0.045],'xlim',[time(1) time(end)],'ylim',bin([1 end]))
colorbar('ytick',-.015:0.015:.045); colormap('hot')
xlabel('Time since ripple detection (s)')
ylabel('Gamma phase')

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_example_ca3_ca3_phase.png', m, d, y);
print('-dpng', savestring)

%ALTER PHASE CALCULATION!
p11 = []; p13 = []; p33 = [];
for e = 1:length(ripc)
	tmp11 = zeros(size(ripc{e}.ca1_ca1_phase,1),1);
    tmp13 = zeros(size(tmp11));
    tmp33 = zeros(size(tmp11));
    for i = 1:size(ripc{e}.ca1_ca1_phase,1)
        tmp11(i) = anglestd(ripc{e}.ca1_ca1_phase(i,:));
        tmp13(i) = anglestd(ripc{e}.ca1_ca3_phase(i,:));
        tmp33(i) = anglestd(ripc{e}.ca3_ca3_phase(i,:));
    end
    p11 = [p11 tmp11];
    p13 = [p13 tmp13];
    p33 = [p33 tmp33];
end

x = -0.4:0.1:0.4;
mean_11 = mean(ca11(6:10:end,:),2); se_11 = std(ca11(6:10:end,:),[],2)./sqrt(size(ca11,2)-1);
mean_13 = mean(ca13(6:10:end,:),2); se_13 = std(ca13(6:10:end,:),[],2)./sqrt(size(ca13,2)-1);
mean_33 = mean(ca33(6:10:end,:),2); se_33 = std(ca33(6:10:end,:),[],2)./sqrt(size(ca33,2)-1);

figure
hold on
bar(1:9,mean_11,'b')
bar(11:19,mean_13,'m')
bar(21:29,mean_33,'r')
legend([{'CA1-CA1'},{'CA1-CA3'},{'CA3-CA3'}])
errorbar2([1:9 11:19 21:29],[mean_11 mean_13 mean_33],[se_11 se_13 se_33],'k')
set(gca,'xtick',[1:9 11:19 21:29],'xticklabel',[x x x])

set(gca,'xtick',time,'xticklabel',time(1):0.1:time(end),'ylim',[-0.1 0.6],'xlim',[-0.5 0.65])
ylabel('Correlation between gamma power and SWR amplitude')
xlabel('Time since ripple detection (s)')
box off
%Note: p1 is significant p<1e-5


%% CA1 and CA3 phase coherence
ca11 = zeros(length(time),length(frequency));
ca13 = zeros(length(time),length(frequency));
ca33 = zeros(length(time),length(frequency));
count11 = 0; count13 = 0; count33 = 0;
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            if ~isempty(f(an).output{d}(e).ca1_ca1_coherence)
                if ~any(any(isnan(f(an).output{d}(e).ca1_ca1_coherence)))
                    ca11 = ca11 + f(an).output{d}(e).ca1_ca1_coherence;
                    count11 = count11+f(an).output{d}(e).nrips;
                end
            end
            if ~isempty(f(an).output{d}(e).ca1_ca3_coherence)
                ca13 = ca13 + f(an).output{d}(e).ca1_ca3_coherence;
                count13 = count13+f(an).output{d}(e).nrips;
            end
            if ~isempty(f(an).output{d}(e).ca3_ca3_coherence)
                ca33 = ca33 + f(an).output{d}(e).ca3_ca3_coherence;
                count33 = count33+f(an).output{d}(e).nrips;
            end
        end
    end
end
ca11 = ca11./count11; ca13 = ca13./count13; ca33 = ca33./count33;

%Plot CA1-CA1 coherence
figure
surf(time,frequency,ca11')
view(0,90); shading interp; colormap hot
set(gca,'clim',[0.5 0.75],'xlim',[time(26) time(86)],'ylim',frequency([1 end]))
colorbar('ytick',0.5:0.05:0.75)


% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1_ca1_coherence_allruns.png', m, d, y);
print('-dpng', savestring)

%Plot CA1-CA3 coherence
figure
surf(time,frequency,ca13')
view(0,90); shading interp; colormap hot
set(gca,'clim',[0.45 0.7],'xlim',[time(26) time(86)],'ylim',frequency([1 end]))
colorbar('ytick',0.45:0.05:0.7)

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca1_ca3_coherence_allruns.png', m, d, y);
print('-dpng', savestring)

%Plot CA3 - CA3 coherence
figure
surf(time,f(1).output{1}(1).frequency,ca33')
view(0,90); shading interp;  colormap hot
set(gca,'clim',[0.5 1],'xlim',[time(26) time(86)],'ylim',f(1).output{1}(1).frequency([1 end]))
colorbar('ytick',0:0.25:1.75)

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_ca3_ca3_coherence_allruns.png', m, d, y);
print('-dpng', savestring)

