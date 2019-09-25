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
%smoothcorrect = smoothcorrect(1:25);
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
        fixarea(i) = fixarea(i) - xvals(cseg-1)/24 + xvals(cseg) * 13/24 + ...
                                    + xvals(cseg+1) * 13/24 - xvals(cseg)/24;
    end
    fixarea(i) = fixarea(i) / 6;
    for (cseg = 18:31)
        % note that the segments start at 1 
        delayarea(i) = delayarea(i) - xvals(cseg-1)/24 + xvals(cseg) * 13/24 + ...
                                    + xvals(cseg+1) * 13/24 - xvals(cseg)/24;
    end
    delayarea(i) = delayarea(i) / 14;
end


delayarea = delayarea(1:(nvalidtrials-4));


histtimes = 0:50:2000;
histspikes = zeros(nvalidtrials, length(histtimes)-1);
for (i = 1:nvalidtrials)
    spikes = double(cobj.spikes(cobj.spike_start(validtrial(i)):cobj.spike_end(validtrial(i))));
    stimes = spikes - fixationtimes(validtrial(i));
    [histspikes(i,:) tmp] = hist(stimes, histtimes);
    fixrate(i) = length(find(stimes < 300)) / .3;
    delayrate(i) = length(find((stimes >= 800) & (stimes < 1250))) / .7 - fixrate(i);
end
histspikes = histspikes * 1000 / (histtimes(2) - histtimes(1));
for i = 1:(nvalidtrials-4)
    smoothrate(i) = sum(delayrate(i:(i+4))/5);
end

% make a lsit of adapttimes for the plots
adapttimes = 0:10:1990;
subplot('Position', [.15 .15 .8 .8]);
cmap = colormap(jet);
hold on
% trials 3, 5, 8, 20
i = 1;
h = plot(adapttimes, x(i,:));
set(gca, 'YLim', [0 100], 'LineWidth', 2);
set(h, 'Color', cmap(200,:));
i = 5;
h = plot(adapttimes, x(i,:));
set(h, 'Color', cmap(350,:));
i = 15;
h = plot(adapttimes, x(i,:));
set(h, 'Color', cmap(700,:));
i = 23;
h = plot(adapttimes, x(i,:));
set(h, 'Color', cmap(950,:));
set(gca, 'XTick', [0 300 800 1500]);
[a b c d] = legend('20% correct', '40% correct', '100% (early)', '100% (late)', 1);
set(a, 'FontSize', 18);
set(gca, 'FontSize', 24);

%ylabel('Firing rate');
xlabel('Time (ms)');
h = line([300 300], [0 100]);
set(h, 'Color', [0 0 0]);
h = line([800 800], [0 100]);
set(h, 'Color', [0 0 0]);
h = line([1500 1500], [0 100]);
set(h, 'Color', [0 0 0]);


% get the histogram estimate for the selected trials
%figure
%orient tall
%subplot(4,2,4);
%cmap = colormap(jet);
%hold on
%i = 1;
%h = plot(histtimes(1:end-1)+25, histspikes(i,:));
%set(gca, 'YLim', [0 120], 'LineWidth', 1.5);
%set(h, 'Color', cmap(200,:));
%i = 5;
%h = plot(histtimes(1:end-1)+25, histspikes(i,:));
%set(h, 'Color', cmap(350,:));
%i = 15;
%h = plot(histtimes(1:end-1)+25, histspikes(i,:));
%set(h, 'Color', cmap(700,:));
%i = 23;
%h = plot(histtimes(1:end-1)+25, histspikes(i,:));
%set(h, 'Color', cmap(950,:));
%set(gca, 'XTick', [0 300 800 1500]);
%%set(gca, 'FontSize', 10);

%legend('0% correct', '20% correct', '80% correct', '100% correct', 1);
%ylabel('Firing rate');
%xlabel('Time (ms)');
%h = line([300 300], [0 120]);
%set(h, 'Color', [0 0 0]);
%h = line([800 800], [0 120]);
%set(h, 'Color', [0 0 0]);
%h = line([1500 1500], [0 120]);
%set(h, 'Color', [0 0 0]);

