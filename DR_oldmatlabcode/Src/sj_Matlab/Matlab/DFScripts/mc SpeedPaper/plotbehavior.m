load '/data13/mcarr/VelocityPaper/power.mat'

speed = []; speed_hist = cell(10,1);
linear = [];
bin = [1/4 1/2 1 2 4 8 16 32];
index = [1 5 10];
for an = 1:length(power)
    for day = 1:length(power(an).output)
        if ~isempty(power(an).output{day})
            tmp = power(an).output{day}(1).speed;
            tmpl = power(an).output{day}(1).linearspeed;
            tmpl = tmpl(tmp>1/8);
            tmp = tmp(tmp>1/8);
            speed_hist{day} = [speed_hist{day}; hist(tmp,bin)./length(tmp)];
            speed = [speed; tmp];
            linear = [linear; tmpl];
        end
    end
end

mean_speed = nan(10,length(bin));
se_speed = nan(10,length(bin));

for day = 1:10
    mean_speed(day,:) = mean(speed_hist{day});
    se_speed(day,:) = std(speed_hist{day})./sqrt(size(speed_hist{day},1)-1);
end

figure
hold on
plot(bin,mean_speed(index(1),:),'r',bin,mean_speed(index(2),:),'b',bin,mean_speed(index(3),:),'g')
legend(num2str(index(1)),num2str(index(2)),num2str(index(3)))
fill([bin bin(end:-1:1)],[mean_speed(index(1),:)+se_speed(index(1),:) mean_speed(index(1),end:-1:1)-se_speed(index(1),end:-1:1)],'r','EdgeColor','none')
fill([bin bin(end:-1:1)],[mean_speed(index(2),:)+se_speed(index(2),:) mean_speed(index(2),end:-1:1)-se_speed(index(2),end:-1:1)],'b','EdgeColor','none')
fill([bin bin(end:-1:1)],[mean_speed(index(3),:)+se_speed(index(3),:) mean_speed(index(3),end:-1:1)-se_speed(index(3),end:-1:1)],'g','EdgeColor','none')
set(gca,'xscale','log','xtick',bin,'xticklabel',bin,'xlim',[bin(1) bin(end)])
xlabel('Speed cm/sec')
ylabel('Proportion of time spent')
box off

speedbin = lookup(log(speed),log(bin));
speedcorr = nan(length(bin),3);
for b = 1:length(bin)
     [R P RL RU]= corrcoef(speed(speedbin==b),linear(speedbin==b));
     speedcorr(b,:) = [R(1,2) RL(1,2) RU(1,2)];
end

figure
plot(bin,speedcorr(:,1),'k')
hold on
fill([bin bin(end:-1:1)],[speedcorr(:,2)' speedcorr(end:-1:1,3)'],'r','EdgeColor','none')
set(gca,'xscale','log','xtick',bin,'xticklabel',bin,'xlim',[bin(1) bin(end)])
xlabel('Speed cm/sec')
ylabel('Correlation with linear speed')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_linear_2d_speedcorrelation.pdf', m, d, y);
print('-dpdf', savestring)

%% RUN FILTER TO LOOK AT SPEED AUTO CORRELATION

%Animal Selection
animals = {'Conley','Corriander','Eight','Five','Miles','Ten'};

% Epoch selection
epochfilter = [];

for i = 1:4
    epochfilter{i} = ['($exposure == ',num2str(i),') & isequal($description,''TrackA'')'];
end

iterator = 'singleepochanal';

f = createfilter('animal',animals,'epochs',epochfilter,'iterator', iterator);

time_bin = 2.5;
f = setfilterfunction(f, 'calcspeedautocorrelation', {'pos'},'bin',time_bin);

f = runfilter(f);

x = f(1).output{1}(1).lags; y = [];
for an = 1:length(f)
    for d = 1:length(f(an).output)
        for e = 1:length(f(an).output{d})
            ind = lookup(x,f(an).output{d}(e).lags);
            y = [y f(an).output{d}(e).xcorr(ind)];
        end
    end
end

x = 1/30*x;
y_mean = mean(y,2);
y_U = y_mean + std(y')'./sqrt(size(y,2)-1);
y_L = y_mean - std(y')'./sqrt(size(y,2)-1);

figure
hold on
plot(x,y_mean,'k')
fill([x x(end:-1:1)],[y_U' y_L(end:-1:1)'],'r','EdgeColor','none')
set(gca,'xtick',-2.5:0.5:2.5,'ylim',[0 1],'xlim',[-2.5 2.5])
xlabel('Time Delay (sec)')
ylabel('Correlation')
box off

% Save figure
[y, m, d] = datevec(date);
savestring = sprintf('/home/mcarr/Figures/SpeedPaper/%d_%d_%d_speedautocorrelation.pdf', m, d, y);
print('-dpdf', savestring)

