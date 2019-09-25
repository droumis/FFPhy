trialtimes = cumsum(ntimesteps);
validtrial = find(testID);
nvalidtrials = length(validtrial);
delayarea = zeros(nvalidtrials,1);

% calculate the five trial moving average of percent correct for this stimulus
trialcorrect = correct(find(testID));

%apply a five trial box car filter
%trialfilt = gaussian(1,5)
trialfilt = ones(1,5) / 5;

%smoothcorrect = smoothvect(trialcorrect, trialfilt);
for i = 1:(nvalidtrials-4)
    smoothcorrect(i) = sum(trialcorrect(i:(i+4))/5);
end


for (i = 1:nvalidtrials)
    [xvals tvals] = instantspline(CP.x, CP.t, thetainit, thetahat, trialtimes(validtrial(i)+1));
    for (cseg = 17:31)
        % note that the segments start at 1 
        delayarea(i) = delayarea(i) - xvals(cseg-1)/24 + xvals(cseg) * 13/24 + ...
                                    + xvals(cseg+1) * 13/24 - xvals(cseg)/24;
    end
    %plot(CP.x(1:20), xvals(1:20));
    %pause
    delayarea(i) = delayarea(i) / 8;
end

% smooth the delay area in the same way
for i = 1:(nvalidtrials-4)
    smoothdelay(i) = sum(delayarea(i:(i+4))/5);
end


% get the number of spikes from the original spike train 
for (i = 1:nvalidtrials)
    spikes = double(cobj.spikes(cobj.spike_start(validtrial(i)):cobj.spike_end(validtrial(i))));
    stimes = spikes - fixationtimes(validtrial(i));
    fixrate(i) = length(find((stimes > 0) & (stimes < 300))) / .3;
    fixrate(i) = length(find((stimes > 0) & (stimes < 300))) / .3;
    %delayrate(i) = length(find((stimes >= 800) & (stimes < 1500))) / .7 - fixrate(i);
    stimrate(i) = length(find((stimes >= 300) & (stimes < 800))) / .5; 
    delayrate(i) = length(find((stimes >= 800) & (stimes < 1500))) / .7; 
end
for i = 1:(nvalidtrials-4)
    smoothfixrate(i) = sum(fixrate(i:(i+4))/5);
    smoothdelayrate(i) = sum(delayrate(i:(i+4))/5);
    smoothstimrate(i) = sum(stimrate(i:(i+4))/5);
end
trialnum = 1:(nvalidtrials-4);

