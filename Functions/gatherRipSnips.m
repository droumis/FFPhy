function [index, allNTDataCat, dataByDay] = gatherRipSnips(andays, Foutput)
        index = []; iEpTetData = []; dataByDay = [];
        for iday = 1:length(andays)
            day = andays(iday);
            eps = find(~cellfun(@isempty,{Foutput{day}.index})); %get nonempty eps
            index = [index; Foutput{day}(eps(1)).index];
            for iep = eps;
                iEpTetData = cat(2,iEpTetData,Foutput{day}(iep).data{1,1}); %cat the event snips across eps, days
                dataByDay = cat(1,dataByDay, repmat(day, length(Foutput{day}(iep).data{1,1}),1));
            end
        end
        allNTDataCat = cat(3, iEpTetData{:});
        
        
        %         for iday = 1:length(ixpc.andays{ian})
%             day = ixpc.andays{ian}(iday);
%             eps = find(~cellfun(@isempty,{F(ian).output{day}.index})); %get nonempty eps
%             ixpc.index{ian} = [ixpc.index{ian}; F(ian).output{day}(eps(1)).index];
%             for iep = eps;
%                 iEpTetData = cat(2,iEpTetData,F(ian).output{day}(iep).data{1,1}); %cat the event snips across eps, days
%                 ixpc.dataByDay = cat(1,dataByDay, repmat(day, length(F(ian).output{day}(iep).data{1,1}),1));
%             end
%             %             % trim events too close to epoch boundaries
%             %             idayInds = ixpc.eventstate{ian}.state(:,8) == day;
%             %             numidayevents = length(idayInds);
%             %             if numidayevents ~= size(iEpTetData(iEpTetDataByDay == day),2)
%             %                 error('this needs to be fixed')
%             %                 iEpTetData = iEpTetData(logical(ixpc.removevec{ian}));
%             %                 if size(ixpc.eventstate{ian}.state(:,1),1) ~= size(iEpTetData,2)
%             %                     error('mismatch between number of state info and lfp sanples')
%             %                 end
%             %               end
%         end
%         %         if max(abs(diff(ixpc.eventstate{ian}.state(:,1:2),[],2))) > 0.033;
%         %             disp(sprintf('%d max event-time offset between pos and lfp times is more than 33ms (1 cam frame).. removing offset events', max(abs(diff(eventstate.state(:,1:2),[],2)))))
%         %             removeevents = find(abs(diff(eventstate{ian}.state(:,1:2),[],2)) > 0.033);
%         %             %         removevec = ones(length(ixpc.eventstate{ian}.state(:,1)),1);
%         %             removevec = ones(length(removeevents),1);
%         %             removevec(removeevents) = 0;
%         %             eventstate.state = eventstate.state(logical(removevec),:);
%         %             eventstate.eventStartTimes = eventstate.eventStartTimes(logical(removevec),:);
%         allNTDataCat = cat(3, iEpTetData{:});
%         clear iEpTetData F % save much needed space 
end

