% function [index, allNTDataCat, dataByDay] = gatherRipSnips_spikes(Foutput)
%         index = []; iEpTetData = []; dataByDay = [];
%         
%         % spikes are saved as a struct array for each individual SU/MU
%         % MU is broken up into nElectrode clusters for each ntrode.. so if
%         
%         % for each epoch, add a cell that reshapes an ntrode's MU's 2nd dim
%         % the third dim.. then dim1 (vert) stack the ntrodes
%         
%         % i can use my combinedataby day .. but i should actually make a
%         % combine mu data by epoch instead now
%         
%         for iday = 1:length(andays)
%             day = andays(iday);
%             eps = find(~cellfun(@isempty,{Foutput{day}.index})); %get nonempty eps
%             index = [index; Foutput{day}(eps(1)).index];
%             for iep = eps;
%                 iEpTetData = cat(2,iEpTetData,Foutput{day}(iep).data{1,1}); %cat the event snips across eps, days
%                 dataByDay = cat(1,dataByDay, repmat(day, length(Foutput{day}(iep).data{1,1}),1));
%             end
%         end
%         allNTDataCat = cat(3, iEpTetData{:});
% end
% 
