function out = peakratevschange(pairdata)

rates = [];
ratechange = [];
meanrate = [];
fields = [];
fieldoverlap = [];
for i = 1:length(pairdata)
    
    [peakrate1, peakind] = max(pairdata(i).run1rates);
    peakrate2 = max(pairdata(i).run2rates);
    
    
    tmpfieldoverlap = overlap(pairdata(i).run1rates,pairdata(i).run2rates);
%     if isnan(tmpfieldoverlap)
%         tmpfieldoverlap = 0;
%     end
    
    tmpmeanrate = noNanMean(pairdata(i).run1rates);
    tmpfield = pairdata(i).run1rates(:)';
    %peakrate2 = pairdata(i).run2rates(peakind);
    tmpratechange = peakrate2-peakrate1;
    if ~isnan(tmpratechange)
        rates = [rates; peakrate1];
        ratechange = [ratechange; tmpratechange];
        meanrate = [meanrate; tmpmeanrate];
        fields = stack(fields,tmpfield);
        fieldoverlap = [fieldoverlap; tmpfieldoverlap];
    end
        
end

out.rates = rates;
out.ratechange = ratechange;
out.meanrate = meanrate;
out.fields = fields;
out.overlap = fieldoverlap;