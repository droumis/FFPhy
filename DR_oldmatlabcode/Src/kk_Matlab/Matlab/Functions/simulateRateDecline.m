function rates = simulateRateDecline(changedata,startdata, steps)
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

rates(1,1) = startdata{1}(rateind,2); %also output the mean rate of the cell 
rates(1,2) = startdata{1}(rateind,1); %initial peak rate

for i = 3:steps+1 %iterate through the model
    tmpdata = changedata;
    [trash, changeind] = min(abs(tmpdata(:,1)-rates(i-1)));
    tmprate = rates(i-1)+tmpdata(changeind,2);
    if (tmprate < bounds(1)) %minimum crossings are reset to the minimum bound
        tmprate = bounds(1);
    end
    if ~( (tmprate > bounds(2)) )
        rates(1,i) = tmprate;
    else
        rates(1,i) = nan; %exclude runaway values
    end
    rateind = changeind;
end
