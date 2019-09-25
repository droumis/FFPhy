%% PLOT COHERENCE

%Animals
animals = {'Bond','Frank','Ten'};

%Define epochs
epochfilter = [];
epochfilter{1} = 'isequal($type, ''run'')';
epochfilter{2} = 'isequal($type,''sleep'')';

%Define iterator
iterator = 'epocheegnonreferenceanal';

%Define time filter
timefilter = {{'get2dstate', '($velocity < 4)'}};

%Run filter for CA1 across hemispheres
CA1tetfilter = 'isequal($area,''CA1'') & $numcells> 0 & isequal($hemisphere,''right'')';
fCA1 = createfilter('animal',animals,'epochs',epochfilter,'eegtetrodes',CA1tetfilter,'excludetime',timefilter,'iterator', iterator);
fCA1 = setfilterfunction(fCA1, 'loadriptriggeredcoherence', {'coherence','tetinfo'},'pair','isequal($hemisphere,''left'') & isequal($area,''CA1'') & $numcells>0');
fCA1 = runfilter(fCA1);

%% Plot average coherence
C11 = [];
for an = 1:length(fCA1)
    for d = 1:length(fCA1(an).output)
        for e = 1:length(fCA1(an).output{d})
            if isempty(C11)
                C11 = fCA1(an).output{d}(e).coherence;
            else
                C11 = cat(3,C11,fCA1(an).output{d}(e).coherence);
            end
        end
    end
end

figure
surf(time,freq,squeeze(mean(C11,3))');
view(0,90);
shading interp
set(gca,'xlim',[time(1) time(end)],'ylim',[freq(1) freq(end)])
set(gca,'clim',[0.4 0.8])
xlabel('Time relative to CA1 ripple')
ylabel('Frequency')

%Look at phase
time = fCA1(1).output{1}(1).time;
bin = -pi+pi/40:pi/20:pi-pi/40;
gam = [];
for an = 1:length(fCA1)
    animal_count = 0;
    for d = 1:length(fCA1(an).eegdata{1})
        tmp_gam = [];
        for e = 1:length(fCA1(an).eegdata{1}{d})
            animal_count = animal_count+1;
            if isempty(tmp_gam)
                tmp_gam = fCA1(an).output{1}(animal_count).gamma_phase;
                tet_count = 1;
            else
                tmp_gam = tmp_gam+fCA1(an).output{1}(animal_count).gamma_phase;
                tet_count = tet_count+1;
            end
        end
        if isempty(gam)
            gam = tmp_gam./tet_count;
        elseif ~isempty(tmp_gam)
            gam = cat(3,gam,tmp_gam./tet_count);
        end
    end
end

figure
surf(bin,time,mean(gam,3));
view(0,90);
shading interp
set(gca,'clim',[0 0.2])
set(gca,'xlim',[bin(2) bin(end-1)],'ylim',[time(1) time(end)])
xlabel('Time relative to CA1 ripple')
ylabel('Frequency')


%% Run filter for CA1-CA3 across hemispheres

CA13tetfilter = 'isequal($area,''CA1'') & $numcells> 1';
fCA13 = createfilter('animal',animals,'epochs',epochfilter,'eegtetrodes',CA13tetfilter,'excludetime',timefilter,'iterator', iterator);
fCA13 = setfilterfunction(fCA13, 'loadriptriggeredcoherence', {'coherence','tetinfo','ripples','cellinfo'},'pair','isequal($area,''CA3'') & $numcells>1','hemisphere',1);
fCA13 = runfilter(fCA13);

save('/data21/mcarr/RipplePaper/CA1CA3coherenceacrosshemispheres.mat','fCA13')

%% Plot gamma coherence for across hemispheres
C13 = []; baseline = [];
for an = 1:length(fCA13)
    ind = [0 0 0];
    for d = 1%:length(fCA13(an).output)
        for e = 1:length(fCA13(an).output{d})
            if ~isempty(fCA13(an).output{d}(e).coherence)                
                if all(ind) == 0
                    freq = fCA13(an).output{d}(e).frequency;
                    ind = [an fCA13(an).output{d}(e).index([1 2])];
                    tmp_baseline = fCA13(an).output{d}(e).coherence_baseline;
                    tmp = fCA13(an).output{d}(e).coherence;
                    tetcount = 1;
                elseif all(ind == [an fCA13(an).output{d}(e).index([1 2])])
                    if size(tmp,2) == size(fCA13(an).output{d}(e).coherence,3)
                        tmp = tmp + fCA13(an).output{d}(e).coherence;
                        tmp_baseline = tmp_baseline + fCA13(an).output{d}(e).coherence_baseline;
                    elseif size(tmp,2) < size(fCA13(an).output{d}(e).coherence,3)
                        index = lookup(fCA13(an).output{d}(e-1).ripples,fCA13(an).output{d}(e).ripples);
                        tmp = tmp + fCA13(an).output{d}(e).coherence(:,index);
                        tmp_baseline = tmp_baseline + fCA13(an).output{d}(e).coherence_baseline;
                   else
                        index = lookup(fCA13(an).output{d}(e).ripples,fCA13(an).output{d}(e-1).ripples);
                        tmp = tmp(:,index) + fCA13(an).output{d}(e).coherence;
                        tmp_baseline = tmp_baseline + fCA13(an).output{d}(e).coherence_baseline;
                   end
                    tetcount = tetcount+1;
                else
                    baseline = [baseline tmp_baseline./tetcount];
                    C13 = cat(2,C13,tmp./tetcount);
                    ind = [an fCA13(an).output{d}(e).index([1 2])];
                    tmp_baseline = fCA13(an).output{d}(e).coherence_baseline;
                    tmp = fCA13(an).output{d}(e).coherence;
                    tetcount = 1;
                end
            end
        end
    end
end

baseline = mean(baseline);
%Run baseline: 0.5453
%Sleep baseline: 0.5743

mean_13 = mean(C13,2);
se_13 = std(C13,[],2)./sqrt(size(C13,2)-1);

x = -0.4:0.1:0.4;

figure
hold on
bar(x,mean_13,'b')
legend('CA1-CA3','location','NorthWest')
errorbar2(x,mean_13,se_13,'k')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.01 0.075],...
    'xlim',[-0.5 0.5],'ytick',0:0.025:0.075)
ylabel('delta CA1-CA3 coherence')
xlabel('Time since ripple detection (s)')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_3f.pdf', m, d, y);
print('-dpdf', savestring)

%Test significance
group = ones(size(C13)); 
for i = 1:size(C13,1)
    group(i,:) =i;
end
[p table stats] = kruskalwallis(reshape(C13,size(C13,1)*size(C13,2),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);
%RUN: Time 0-400ms > baseline at p<1e-5
%SLEEP: Time -100-400ms > baseline at p<1e-5

%% Plot gamma phase locking for across hemispheres
C13 = []; baseline = [];
for an = 1:length(fCA13)
    ind = [0 0 0];
    for d = 1%:length(fCA13(an).output)
        for e = 1:length(fCA13(an).output{d})
            if ~isempty(fCA13(an).output{d}(e).coherence)                
                if all(ind) == 0
                    freq = fCA13(an).output{d}(e).frequency;
                    ind = [an fCA13(an).output{d}(e).index([1 2])];
                    tmp_baseline =  fCA13(an).output{d}(e).gamma_phase_baseline;
                    tmp = fCA13(an).output{d}(e).gamma_phase;
                    tetcount = 1;
                elseif all(ind == [an fCA13(an).output{d}(e).index([1 2])])
                	tmp = cat(2,tmp,fCA13(an).output{d}(e).gamma_phase);
                    tmp_baseline = cat(2,tmp_baseline,fCA13(an).output{d}(e).gamma_phase_baseline);
                    tetcount = tetcount+1;
                else
                    C13 = cat(2,C13,mean(tmp,2));
                    baseline = [baseline; mean(tmp_baseline)];
                    ind = [an fCA13(an).output{d}(e).index([1 2])];
                    tmp = fCA13(an).output{d}(e).gamma_phase;
                    tmp_baseline =  fCA13(an).output{d}(e).gamma_phase_baseline;
                    tetcount = 1;
                end
            end
        end
    end
end
baseline = mean(baseline);
%RUN baseline: 0.698
%SLEEP baseline: 0.7016

mean_13 = mean(C13,2);
se_13 = std(C13,[],2)./sqrt(size(C13,2)-1);

x = -0.4:0.1:0.4;

figure
hold on
bar(x,mean_13,'b')
legend('CA1-CA3','location','NorthWest')
errorbar2(x,mean_13,se_13,'k')
set(gca,'xtick',x,'xticklabel',[-0.4:0.1:0.4 100],'ylim',[-0.025 0.125],...
    'xlim',[-0.5 0.5],'ytick',0:0.025:0.125)
ylabel('delta CA1-CA3 phase locking')
xlabel('Time since ripple detection (s)')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/RipplePaper/%d_%d_%d_3c.pdf', m, d, y);
print('-dpdf', savestring)

%Test significance
group = ones(size(C13)); 
for i = 1:size(C13,1)
    group(i,:) =i;
end
[p table stats] = kruskalwallis(reshape(C13,size(C13,1)*size(C13,2),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);
%RUN: Time 100ms > baseline at p<1e-5
% 0ms > baseline at p<0.01    
%SLEEP: Time 0-400ms > baseline at p<1e-5

%% Run filter for CA3-CA3 across hemispheres
CA33tetfilter = 'isequal($area,''CA3'') & $numcells> 1 & isequal($hemisphere,''right'')';
fCA33 = createfilter('animal',animals,'epochs',epochfilter,'eegtetrodes',CA33tetfilter,'excludetime',timefilter,'iterator', iterator);
fCA33 = setfilterfunction(fCA33, 'loadriptriggeredcoherence', {'coherence','tetinfo'},'pair','isequal($hemisphere,''left'') & isequal($area,''CA3'') & $numcells>1');
fCA33 = runfilter(fCA33);

% Plot average coherence
C33 = [];
for an = 1:length(fCA33)
    for d = 1:length(fCA33(an).output)
        for e = 1:length(fCA33(an).output{d})
            if isempty(C33)
                time = fCA33(an).output{d}(e).time;
                freq = fCA33(an).output{d}(e).frequency;
                C33 = fCA33(an).output{d}(e).coherence;
            else
                C33 = cat(3,C33,fCA33(an).output{d}(e).coherence);
            end
        end
    end
end
C33 = squeeze(mean(C33,3));
figure
surf(time,freq,C33');
view(0,90);
shading interp
set(gca,'clim',[0.4 0.8])
set(gca,'xlim',[time(1) time(end)],'ylim',[freq(1) freq(end)])
xlabel('Time relative to CA1 ripple')
ylabel('Frequency')

%Look at phase
time = fCA33(1).output{1}(1).time;
bin = -pi+pi/40:pi/20:pi-pi/40;
gam = [];
for an = 1:length(fCA33)
    animal_count = 0;
    for d = 1:length(fCA33(an).eegdata{1})
        tmp_gam = [];
        for e = 1:length(fCA33(an).eegdata{1}{d})
            animal_count = animal_count+1;
            if isempty(tmp_gam)
                tmp_gam = fCA33(an).output{1}(animal_count).gamma_phase;
                tet_count = 1;
            else
                tmp_gam = tmp_gam+fCA33(an).output{1}(animal_count).gamma_phase;
                tet_count = tet_count+1;
            end
        end
        if isempty(gam)
            gam = tmp_gam./tet_count;
        else
            gam = cat(3,gam,tmp_gam./tet_count);
        end
    end
end

figure
surf(bin,time,mean(gam,3));
view(0,90);
shading interp
set(gca,'clim',[0 0.3])
set(gca,'xlim',[bin(2) bin(end-1)],'ylim',[time(1) time(end)])
xlabel('Time relative to CA1 ripple')
ylabel('Frequency')
