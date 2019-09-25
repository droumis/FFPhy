function info = calclininfo(trajData)
%          info - calclininfo(trajData)
%
%          Returns the place information carried by the cell in bits using the
%          formula info = sum(P * (Ri/Rm) * log(Ri/Rm)) from Skaggs (1993)
%          The place information is computed based on the trajectory firing
%          rate maps
%          trajData = a n by 5 matrix from a single trajectory (output of
%          calclinfieds)


goodValues = find(~isnan(trajData(:,5)));
trajData = trajData(goodValues,:);

occupancy = trajData(:,2); %occupancy for each bin
Poccupancy = occupancy/sum(occupancy);
firingRate = trajData(:,5); %binned occ-normalized firing rate
meanrate = sum(firingRate.*occupancy)/sum(occupancy);

goodvalues2 = find(firingRate > 0); %we ignore bins with 0 firing rate because they add nothing to the final sum

firingRate = firingRate(goodvalues2);
Poccupancy = Poccupancy(goodvalues2);


rateRatio = firingRate/meanrate;

info =  sum((Poccupancy.*rateRatio) .* log2(rateRatio));
