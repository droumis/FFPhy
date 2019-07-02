function [index, allNTDataCat, dataByDay] = gatherRipSnips(andays, Foutput)

        index = []; iEpTetData = []; dataByDay = [];
        for iday = 1:length(andays)
            day = andays(iday);
            eps = find(~cellfun(@isempty,{Foutput{day}.index})); %get nonempty eps
            index = [index; Foutput{day}(eps(1)).index];
            for iep = eps
                iEpTetData = cat(2,iEpTetData,Foutput{day}(iep).data{1,1}); %cat the event snips across eps, days
                dataByDay = cat(1,dataByDay, repmat(day, length(Foutput{day}(iep).data{1,1}),1));
            end
        end
        allNTDataCat = cat(3, iEpTetData{:});
end

