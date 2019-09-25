trialtimes = cumsum(ntimesteps);
validtrial = find(testID);
nvalidtrials = length(validtrial);
delayarea = zeros(nvalidtrials,1);
fixarea = zeros(nvalidtrials,1);

% calculate the five trial moving average of percent correct for this stimulus
trialcorrect = correct(find(testID));

%apply a five trial box car filter
%trialfilt = gaussian(1,9);
trialfilt = ones(1,5) / 5;

%smoothcorrect = smoothvect(trialcorrect, trialfilt);
%smoothcorrect = smoothcorrect(1:24);
for i = 1:(nvalidtrials-4)
    smoothcorrect(i) = sum(trialcorrect(i:(i+4))/5);
end

% compute the area in the delay period for all trials

xvals = zeros(nvalidtrials, length(CP.x));
x = zeros(nvalidtrials, 200);
for (i = 1:nvalidtrials)
    [xvals(i,:) tvals] = instantspline(CP.x, CP.t, thetainit, thetahat, trialtimes(validtrial(i)+1));
    tmp = cardinal(CP.x(1:40), xvals(i,1:40), 200);
    x(i,:) = tmp(:,2)';
    for cseg = 2:7
        % note that the segments start at 1 
        fixarea(i) = fixarea(i) - xvals(i,cseg-1)/24 + xvals(i,cseg) * 13/24 + ...
                                    + xvals(i,cseg+1) * 13/24 - xvals(i,cseg)/24;
    end
    fixarea(i) = fixarea(i) / 6;
    for (cseg = 18:31)
        % note that the segments start at 1 
        delayarea(i) = delayarea(i) - xvals(i,cseg-1)/24 + xvals(i,cseg) * 13/24 + ...
                                    + xvals(i,cseg+1) * 13/24 - xvals(i,cseg)/24;
    end
    delayarea(i) = delayarea(i) / 14;
end


delayarea = delayarea(1:(nvalidtrials-4));


% get the number of spikes from the original spike train 

times = 0:50:2000;
histspikes = zeros(nvalidtrials, length(times)-1);
for (i = 1:nvalidtrials)
    spikes = double(cobj.spikes(cobj.spike_start(validtrial(i)):cobj.spike_end(validtrial(i))));
    stimes = spikes - fixationtimes(validtrial(i));
    [histspikes(i,:) tmp] = hist(stimes, times);
    fixrate(i) = length(find(stimes < 300)) / .3;
    delayrate(i) = length(find((stimes >= 800) & (stimes < 1240))) / .7 - fixrate(i);
end
histspikes = histspikes * 1000 / (times(2) - times(1));
for i = 1:(nvalidtrials-4)
    smoothrate(i) = sum(delayrate(i:(i+4))/5);
end

% get the histogram estimate for the selected trials
figure
orient tall
% make a list of times for the plots
adapttimes = 0:10:1990;
subplot(4,2,1)
cmap = colormap(jet);
hold on
i = 1;
h = plot(adapttimes, x(i,:));
set(gca, 'YLim', [0 120], 'LineWidth', 1.5);
set(h, 'Color', cmap(200,:));
i = 18;
h = plot(adapttimes, x(i,:));
set(h, 'Color', cmap(350,:));
i = 26;
%i = 24;
h = plot(adapttimes, x(i,:));
set(h, 'Color', cmap(700,:));
i = 40;
h = plot(adapttimes, x(i,:));
set(h, 'Color', cmap(950,:));
set(gca, 'XTick', [0 300 800 1500]);

ylabel('Firing rate');
[a b c d] = legend('0% correct', '20% correct', '100% (early)', '100% (late)', 2);
h = line([300 300], [0 120]);
set(h, 'Color', [0 0 0]);
h = line([800 800], [0 120]);
set(h, 'Color', [0 0 0]);
h = line([1500 1500], [0 120]);
set(h, 'Color', [0 0 0]);

subplot(4,2,3)
cmap = colormap(jet);
hold on
i = 1;
h = plot(times(1:end-1)+24, histspikes(i,:));
set(gca, 'YLim', [0 120], 'LineWidth', 1.5);
set(h, 'Color', cmap(200,:));
i = 18;
h = plot(times(1:end-1)+24, histspikes(i,:));
set(h, 'Color', cmap(350,:));
i = 26;
%i = 24;
h = plot(times(1:end-1)+24, histspikes(i,:));
set(h, 'Color', cmap(700,:));
i = 40;
h = plot(times(1:end-1)+24, histspikes(i,:));
set(h, 'Color', cmap(950,:));
set(gca, 'XTick', [0 300 800 1500]);
set(a, 'FontSize', 10);

%legend('0% correct', '20% correct', '80% correct', '120% correct', 1);
ylabel('Firing rate');
xlabel('Time (ms)');
h = line([300 300], [0 120]);
set(h, 'Color', [0 0 0]);
h = line([800 800], [0 120]);
set(h, 'Color', [0 0 0]);
h = line([1500 1500], [0 120]);
set(h, 'Color', [0 0 0]);

