function rates = simulateRateDecline2(changedata,startdata, steps)
%rates = simulateRateDecline(changedata,startdata, steps)
%
%Picks a value at random from STARTDATA and uses CHANGEDATA as a lookup
%table for how the value should be shanged.  This is done by picking the
%closest value in the first column of CHANGEDATA.  After each iteration,
%the new value is saved in the output and is again compared to CHANGEDATA
%for the next iteration. 
%CHANGEDATA - N by 2, first column is peak rate, second column is change in peak rate
%STARTDATA - a cell with an R by 1 matrix of startvalues for peak rate
%STEPS - the number of steps to iterate through

bounds = [-.01 70]; %Runaway bounds - any values that go past the 2nd bound are excluded
rateind = randsample(1:size(startdata{1},1),1); %the start cell is picked at random

rates(1,1) = noNanMean(startdata{1}(rateind,:)); %output the mean rate of the cell 

tmprates = startdata{1}(rateind,:);
for i = 1:steps %iterate through the model
    for j = 1:length(tmprates)
        tmpdata = changedata;
        currrate = tmprates(j);
        if (~isnan(currrate) & (currrate > .1))

            [trash, changeind] = min(abs(tmpdata(:,1)-currrate));
            currrate = currrate+tmpdata(changeind,2);
            if (currrate < bounds(1)) %minimum crossings are reset to the minimum bound
                currrate = bounds(1);
            end
            if (currrate > bounds(2)) 
                 currrate = nan; %exclude runaway values
            
            end
            tmprates(j) = currrate;
        end
    end
end
rates(1,2) = noNanMean(tmprates);
