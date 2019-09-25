function out = ratevschange(pairdata)
%out = ratevschange(pairdata)
%Calculates the change in rate in each bin, as a function of the original
%rate in the bin


rates = [];
ratechange = [];
for i = 1:length(pairdata)
    
    tmprates = pairdata(i).run1rates;
    tmpratechange = pairdata(i).run2rates-pairdata(i).run1rates;
    validchange = find(~isnan(tmpratechange));
    
    rates = [rates; tmprates(validchange)];
    ratechange = [ratechange; tmpratechange(validchange)];
        
end

out.rates = rates;
out.ratechange = ratechange;