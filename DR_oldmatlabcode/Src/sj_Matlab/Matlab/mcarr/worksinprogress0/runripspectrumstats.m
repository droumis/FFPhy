%Animals
animals = {'Bond','Frank','Ten'};

%Define epochs
epochfilter = [];
epochfilter{1} = 'isequal($type, ''run'')';

CA1tetfilter = 'isequal($area,''CA1'') & $numcells>1';
CA1tetfilter = 'isequal($area,''CA1'') & $numcells> 0 & isequal($hemisphere,''right'')';
CA3tetfilter = {'isequal($area, ''CA3'') & $numcells>1', 'isequal($area, ''CA1'') & $numcells>1'};

%Define iterator
iterator = 'eegnonreferenceanal';

%Define time filter
timefilter = {{'get2dstate', '($velocity < 4)'}};

%Run filter for CA1
f1 = createfilter('animal',animals,'epochs',epochfilter,'eegtetrodes',CA1tetfilter,'excludetime',timefilter,'iterator', iterator);
f1 = setfilterfunction(f1, 'calcripspectrumstats', {'spectrum'});
f1 = runfilter(f1);

% %Run filter for CA1 across hemispheres
% fCA1 = createfilter('animal',animals,'epochs',epochfilter,'eegtetrodes',CA1tetfilter,'excludetime',timefilter,'iterator', iterator);
% fCA1 = setfilterfunction(fCA1, 'loadriptriggeredcoherence', {'coherence','tetinfo'},'pair','isequal($hemisphere,''left'') & isequal($area,''CA1'') & $numcells>0');
% fCA1 = runfilter(fCA1);

%Run filter for CA3
f3 = createfilter('animal',animals,'epochs',epochfilter,'eegtetrodepairs',CA3tetfilter,'excludetime',timefilter,'iterator', iterator);
f3 = setfilterfunction(f3, 'calcripspectrumstats', {'spectrum'});
f3 = runfilter(f3);


%% Plot average spectrum for CA1
S = []; rip_count = 0;
for an = 1:length(f1)
    animal_count = 0;
    for d = 1:length(f1(an).eegdata{1})
        tmp = []; tmp_peak = [];
        for e = 1:length(f1(an).eegdata{1}{d})
            animal_count = animal_count+1;
            if isempty(tmp)
                tmp = f1(an).output{1}(animal_count).total;
                rip_count = rip_count+f1(an).output{1}(animal_count).nrips;
                tet_count = 1;
            else
                tet_count = tet_count+1;
                tmp = tmp+f1(an).output{1}(animal_count).total;
            end
        end
        if isempty(S)
        	S = tmp./tet_count;
            time = f1(an).output{d}(e).time;
            freq = f1(an).output{d}(e).freq;
        elseif ~isempty(tmp)
        	S = S +tmp./tet_count;
        end
    end
end
S = S./rip_count;
figure
surf(time,freq,S');
view(0,90);
shading interp
colormap bone
set(gca,'clim',[0 1.25],'xlim',[time(1) time(end)],'ylim',[freq(1) freq(end)])

%% Look at gamma power as compared to ripples for CA1

power = []; peak = [];
for an = 1:length(f1)
    animal_count = 0;
    for d = 1:length(f1(an).eegdata{1})
        tmp = []; tmp_peak = [];
        for e = 1:length(f1(an).eegdata{1}{d})
            animal_count = animal_count+1;
            if isempty(tmp)
                tmp = f1(an).output{1}(animal_count).power;
                tmp_peak = f1(an).output{1}(animal_count).peak;
                tet_count = 1;
            else
                tet_count = tet_count+1;
                tmp = tmp+f1(an).output{1}(animal_count).power;
                tmp_peak = tmp_peak+f1(an).output{1}(animal_count).peak;
            end
        end
        peak = [peak; tmp_peak./tet_count];
        power = [power; tmp./tet_count];
    end
end

%Plot gamma power before ripple onset, at the peak ripple power, and after.
figure
hold on
bar([1 2 3],mean(power),'b')
errorbar2([1 2 3],mean(power),std(power)./sqrt(size(power,1)-1),'k')
set(gca,'ylim',[0 1],'xlim',[0.5 3.5],'xtick',[1 2 3],'xticklabel',{'Before','Peak of Ripple','After'})
ylabel('Normalized gamma power')

%Statistics
group = ones(size(power)); group(:,2) = 2; group(:,3) = 3;
[p table stats] = kruskalwallis(reshape(power,size(power,1)*size(power,2),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);

%Plot correlation between gamma power and peak ripple power
r_power = nan(3,4);
for i = 1:3
    [r p rlo rup] = corrcoef(peak,power(:,i));
    r_power(i,:) = [rlo(1,2)-r(1,2) r(1,2) rup(1,2)-r(1,2) p(1,2)];
end
figure
hold on
bar([1 2 3],r_power(:,2))
errorbar2([1 2 3],r_power(:,2),r_power(:,[1 3])')
set(gca,'ylim',[0 0.5],'xlim',[0.5 3.5],'xtick',[1 2 3],'xticklabel',{'Before','Peak of Ripple','After'})
ylabel('Correlation between gamma power and peak ripple power')

%% Look at correlation between gamma and ripple power for CA1
c = []; lags = [];
for an = 1:length(f1)
    animal_count = 0;
    for d = 1:length(f1(an).eegdata{1})
        tmp = [];
        for i = 1:length(f1(an).eegdata{1}{d})
            animal_count = animal_count+1;
            if isempty(tmp)
                tmp = f1(an).output{1}(animal_count).correlation;
                tet_count = 1;
            else
                tet_count = tet_count+1;
                tmp = tmp + f1(an).output{1}(animal_count).correlation;
            end
        end
        if isempty(lags)
        	lags = f1(an).output{1}(count).lags;
        end
        c = [c; tmp./tet_count];
    
    end
end

c_se = std(c)./sqrt(size(c,1)-1);
c_L = mean(c) - c_se;
c_U = mean(c) + c_se;

figure
plot(lags,mean(c),'k')
hold on
fill([lags lags(end:-1:1)],[c_L c_U(end:-1:1)],'m','EdgeColor','none')

%% Plot average spectrum for CA3
S = []; rip_count = 0;
for an = 1:length(f3)
    animal_count = 0;
    for d = 1:length(f3(an).eegdata{1})
        tmp = []; tmp_peak = []; tet = [];
        for e = 1:size(f3(an).eegdata{1}{d},1)
            animal_count = animal_count+1;
            if isempty(tet)
                tet = f3(an).eegdata{1}{d}(e,1);
                tmp = f3(an).output{1}(animal_count).total;
                rip_count = rip_count+f3(an).output{1}(animal_count).nrips;
                tet_count = 1;
            elseif f3(an).eegdata{1}{d}(e,1) ~= tet
                tet_count = tet_count+1;
                tmp = tmp+f3(an).output{1}(animal_count).total;
                tet = f3(an).eegdata{1}{d}(e);
            end
        end
        if isempty(S)
        	S = tmp./tet_count;
            time = f3(an).output{d}(e).time;
            freq = f3(an).output{d}(e).freq;
        elseif ~isempty(tmp)
        	S = S +tmp./tet_count;
        end
    end
end
S = S./rip_count;
figure
surf(time,freq,S');
view(0,90);
shading interp
colormap bone
set(gca,'clim',[0 1.25],'xlim',[time(1) time(end)],'ylim',[freq(1) freq(end)])

%% Look at CA3 gamma power as compared to CA1 ripples
power = []; peak = [];
for an = 1:length(f3)
    animal_count = 0;
    for d = 1:length(f3(an).eegdata{1})
        tmp = []; tmp_peak = []; tet = [];
        for e = 1:size(f3(an).eegdata{1}{d},1)
            animal_count = animal_count+1;
            if isempty(tet)
                tet = f3(an).eegdata{1}{d}(e,1);
                tmp = f3(an).output{1}(animal_count).power;
                tmp_peak = f3(an).output{1}(animal_count).peak;
                tet_count = 1;
            elseif f3(an).eegdata{1}{d}(e,1) ~= tet
                tet_count = tet_count+1;
                tmp = tmp+f3(an).output{1}(animal_count).power;
                tmp_peak = tmp_peak + f3(an).output{1}(animal_count).peak;
                tet = f3(an).eegdata{1}{d}(e);
            end
        end
        peak = [peak; tmp_peak./tet_count];
        power = [power; tmp./tet_count];
    end
end


%Plot gamma power before ripple onset, at the peak ripple power, and after.
figure
hold on
bar([1 2 3],mean(power),'b')
errorbar2([1 2 3],mean(power),std(power)./sqrt(size(power,1)-1),'k')
set(gca,'ylim',[0 1],'xlim',[0.5 3.5],'xtick',[1 2 3],'xticklabel',{'Before','Peak of Ripple','After'})
ylabel('Normalized gamma power')

%Statistics
group = ones(size(power)); group(:,2) = 2; group(:,3) = 3;
[p table stats] = kruskalwallis(reshape(power,size(power,1)*size(power,2),1),reshape(group,size(group,1)*size(group,2),1));
c = multcompare(stats);

%Plot correlation between CA3 gamma power and peak CA1 ripple power
r_power = nan(3,4);
for i = 1:3
    [r p rlo rup] = corrcoef(peak,power(:,i));
    r_power(i,:) = [rlo(1,2)-r(1,2) r(1,2) rup(1,2)-r(1,2) p(1,2)];
end
figure
hold on
bar([1 2 3],r_power(:,2))
errorbar2([1 2 3],r_power(:,2),r_power(:,[1 3])')
set(gca,'ylim',[0 0.5],'xlim',[0.5 3.5],'xtick',[1 2 3],'xticklabel',{'Before','Peak of Ripple','After'})
ylabel('Correlation between gamma power and peak ripple power')

%% Look at correlation between gamma and ripple power for CA1
c = []; lags = [];
for an = 1:length(f3)
    animal_count = 0;
    for d = 1:length(f3(an).eegdata{1})
        if ~isempty(f3(an).eegdata{1}{d})
            ca1_count = length(unique(f3(an).eegdata{1}{d}(:,2)));
            ca3_count = length(unique(f3(an).eegdata{1}{d}(:,1)));
            tmp = cell(ca3_count);
            tet_index = lookup(f3(an).eegdata{1}{d}(:,1),unique(f3(an).eegdata{1}{d}(:,1)));
            for e = 1:length(f3(an).eegdata{1}{d})
                animal_count = animal_count+1;
                if isempty(tmp{tet_index(e)})
                    tmp{tet_index(e)} = f3(an).output{1}(animal_count).correlation;            
                else
                    tmp{tet_index(e)} = [f3(an).output{1}(animal_count).correlation + tmp{tet_index(e)}];
                end
            end
            if isempty(lags)
                lags = f3(an).output{1}(animal_count).lags;
            end
            tmp_corr = tmp{1}./ca1_count;
            for e = 2:ca3_count
                tmp_corr = tmp_corr + tmp{e}./ca1_count;
            end
            c = [c; tmp_corr./ca3_count];
        end
    end
end

c_se = std(c)./sqrt(size(c,1)-1);
c_L = mean(c) - c_se;
c_U = mean(c) + c_se;

figure
plot(lags,mean(c),'k')
hold on
fill([lags lags(end:-1:1)],[c_L c_U(end:-1:1)],'m','EdgeColor','none')
set(gca,'xlim',[lags(1) lags(end)],'ylim',[-0.005 0.12])