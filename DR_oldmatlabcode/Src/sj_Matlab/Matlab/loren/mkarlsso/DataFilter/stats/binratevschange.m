function out = binratevschange(pairdata)


ratechange = [];
rates = [];
indchange = [];

for i = 1:length(pairdata)
    
    
    %validind = find((pairdata(i).run1rates > 1) & (pairdata(i).run1rates < 12));
    validind = find((pairdata(i).run1rates > 12));
    %validind = ((pairdata(i).run1rates > 1));
    
    %validind2 = find((pairdata(i).run2rates > 1) & (pairdata(i).run2rates < 12) );
    validind2 = find((pairdata(i).run2rates > 12));
    tmpratechange = pairdata(i).run2rates(validind) - pairdata(i).run1rates(validind);
    tmprate = pairdata(i).run1rates(validind);
    
    
    if ~isempty(tmpratechange)
        indchange = [indchange; (length(validind2)-length(validind))/length(validind2)];
        ratechange = [ratechange; tmpratechange];
        rates = [rates; tmprate];
        
    end
        
end


out.ratechange = ratechange;
out.rates = rates;
out.indchange = indchange;
