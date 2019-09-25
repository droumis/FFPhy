function out = overlap_threeepochs(pairdata)

ratesB1 = [];
ratesB2 = [];
ratesA = [];
ratechange = [];
meanrate = [];
fields = [];
fieldoverlapBB = [];
fieldoverlapBA = [];
for i = 1:length(pairdata)
    
    [peakrate1, peakind] = max(pairdata(i).run1rates);
    peakrate2 = max(pairdata(i).run2rates);
    peakrate3 = max(pairdata(i).run3rates);
    
    
    tmpfieldoverlap1 = overlap(pairdata(i).run1rates,pairdata(i).run2rates);
    tmpfieldoverlap2 = overlap(pairdata(i).run2rates,pairdata(i).run3rates);
%     if isnan(tmpfieldoverlap)
%         tmpfieldoverlap = 0;
%     end
    
    tmpmeanrate = noNanMean(pairdata(i).run1rates);
    tmpfield = pairdata(i).run1rates(:)';
    %peakrate2 = pairdata(i).run2rates(peakind);
    tmpratechange = peakrate2-peakrate1;
    if ~isnan(tmpratechange)
        ratesB1 = [ratesB1; peakrate1];
        ratesB2 = [ratesB2; peakrate2];
        ratesA = [ratesA; peakrate3];
        ratechange = [ratechange; tmpratechange];
        meanrate = [meanrate; tmpmeanrate];
        fields = stack(fields,tmpfield);
        fieldoverlapBB = [fieldoverlapBB; tmpfieldoverlap1];
        fieldoverlapBA = [fieldoverlapBA; tmpfieldoverlap2];
    end
        
end

out.ratesA = ratesA;
out.ratesB1 = ratesB1;
out.ratesB2 = ratesB2;
out.ratechange = ratechange;
out.meanrate = meanrate;
out.fields = fields;
out.overlapBB = fieldoverlapBB;
out.overlapBA = fieldoverlapBA;