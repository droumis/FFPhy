%% Define Filter Options

% Animal selection
animals = {'Dudley','Miles','Bond','Frank','Ten','Five','Eight'};
%animals = {'Five','Eight'};

% Epoch selection
epochfilter = [];
for i = 1:23
    epochfilter{i} = ['($exposure == ',num2str(i),') & isequal($description, ''TrackA'')'];
end

% Tetrode pair selection
tetfilter = {'(isequal($area, ''CA1'') & ($representative==1))','(isequal($area, ''CA1'') )'};
%tetfilter = {'(isequal($area, ''CA1'') & ($representative == 1))','($layer == 3) & ($representative==1)' };

%Iterator Selection
iterator = 'eeganal';

% Create and Run Filter
f = createfilter('animal',animals,'epochs',epochfilter, 'eegtetrodepairs',tetfilter,'iterator', iterator);

f = setfilterfunction(f, 'calceegxcorr', {'eeg'});

f = runfilter(f);

%% Plot Example

y = zeros(length(f),1);
for an = 1:length(f)
    y(an) = length(f(an).output);
end 
x = find(y==max(y));
temp = cell(length(f(x(1)).output),1);
for an = 1:length(f)
    for d=1:length(f(an).output)
        for e=1:length(f(an).output{d})
            if isempty(temp{d})
                temp{d} = f(an).output{d}(e).xcorr;   
            else
                temp{d} = [temp{d}; f(an).output{d}(e).xcorr];
            end
        end
    end
end

test= [];
for d=1:length(temp)
    test = [test; temp{d}];
end
figure
plot(f(1).output{1}(1).lags/1500,mean(test),'k','LineWidth',2)

%% Plot each epochs
figure
color = ['r'; 'm'; 'y'; 'g'; 'c'; 'b'; 'k';'r';'m';'y';'g';'c';'b';'k'];
for d=1:10
    plot(f(1).output{1}.lags/1500,mean(temp{d}),color(d),'LineWidth',2)
    hold on
end
set(gca,'yLim',[-0.45 0.45])
set(gca,'xLim',[f(1).output{1}.lags(1)/1500 f(1).output{1}.lags(end)/1500])
set(gca,'FontSize',18)
xlabel('Time Lag (seconds)','fontsize',18)
ylabel('Cross Correlation','fontsize',18)
title('Cross Correlation between CA1 and CA3','fontsize',18)
legend('Exposure 1','Exposure 2', 'Exposure 3', 'Exposure 4', 'Exposure 5','Exposure 6', 'Exposure 7','Exposure 8')
