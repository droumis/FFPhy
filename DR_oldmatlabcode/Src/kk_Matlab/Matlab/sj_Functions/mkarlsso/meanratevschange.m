function out = meanratevschange(pairdata)

rates = [];
ratechange = [];
for i = 1:length(pairdata)
    
    rate1 = noNanMean(pairdata(i).run1rates);
    rate2 = noNanMean(pairdata(i).run2rates);
    
    %peakrate2 = pairdata(i).run2rates(peakind);
    tmpratechange = rate2-rate1;
    if ~isempty(tmpratechange)
        rates = [rates; rate1];
        ratechange = [ratechange; tmpratechange];
    end
        
end

out.rates = rates;
out.ratechange = ratechange;