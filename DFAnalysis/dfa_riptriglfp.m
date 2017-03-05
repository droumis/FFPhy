function out = dfa_riptriglfp(index, excludeperiods, eeg, ripple, events, varargin)

% gathers lfp around ripple times
% out is a structure with the following fields


win = [0.5 0.5];
appendindex = 1;
eventtype = 'rippleskons';

% process varargin and overwrite default values
if (~isempty(varargin))
    assign(varargin{:});
end



if strcmp(eventtype,'ripples')
    eventtimes = events{index(1)}{index(2)}{index(3)}.starttime;
else % ripplecons, ripplekons, etc
    eventtimes = [];
    eventtimes(:,1) = events{index(1,1)}{index(1,2)}{1}.starttime(:,1);
    eventtimes(:,2) = events{index(1,1)}{index(1,2)}{1}.endtime(:,1);
end
    eventtimes = eventtimes(~isExcluded(eventtimes(:,1),excludeperiods),:);

% Define EEG
e = eeg{index(1,1)}{index(1,2)}{index(1,3)}.data';
num_samp = length(e);
starttime = double(eeg{index(1,1)}{index(1,2)}{index(1,3)}.starttime);
endtime_eeg = double(eeg{index(1,1)}{index(1,2)}{index(1,3)}.endtime);
clockrate = eeg{index(1,1)}{index(1,2)}{index(1,3)}.clockrate;
samprate = eeg{index(1,1)}{index(1,2)}{index(1,3)}.samprate;
endtime = starttime + (num_samp-1) * (1 / samprate);
disp(sprintf('endtime_eeg: %.04f calculated_endtime: %.04f', endtime_eeg/60, endtime/60));

% Define triggering events as the start of each ripple
triggers = eventtimes(:,1)-starttime;

%Remove triggering events that are too close to the beginning or end
while triggers(1)<win(1)
    triggers(1) = [];
end
while triggers(end)> endtime-win(2)
    triggers(end) = [];
end

%update... nan filled any gaps in the lfp data so i can now actually use this method of getting timestamp
    LFPtimes=(starttime:1/samprate:(starttime+(1/samprate)*num_samp))'; % prob should be using the adjusted timestamps instead, no?
    eventStartIndices = lookup(eventtimes(:,1),LFPtimes);
    eventEndIndices = lookup(eventtimes(:,2),LFPtimes);
    windowStartIndices = lookup(eventtimes(:,1)-win(1),LFPtimes);
    windowEndIndices = lookup(eventtimes(:,2)+win(2),LFPtimes);
%     eventStartLFPtimes = LFPtimes(eventStartIndices);
%     windowStartEndTimes = [eventStartLFPtimes(:)-win(1) eventStartLFPtimes(:)+win(2)]; 

    %% STEP 3: Gather the ripple window data from all the regions 
    for currrip=1:length(windowStartIndices) %for each ripple from source region
        clear YripLFPdata
        %         for iArea = 1:length(regions);
        for iNTrode = 1:length(index(:,3)); %for each ntrode
%             nTrodeID = index(iNTrode,3);
            % Gather lfp data for each ntrode within the current rip window
            out.data{currrip}(iNTrode,:) = eeg{index(iNTrode,1)}{index(iNTrode,2)}{index(iNTrode,3)}.data(windowStartIndices(currrip):windowEndIndices(currrip));
        end
%         YripLFPdataMAT{iLFPtype}{currrip} = cell2mat(YripLFPdata); %stack all the traces from all the regions next to each other
%         Xwindowtimes{iLFPtype}{currrip} = WindowStartEndTimes(currrip,1):1/samprate:WindowStartEndTimes(currrip,2);
%         lfptraceLUTregion = cell2mat(lfptraceregion{iLFPtype}'); % [1 x length(nTrodes)] vec mapping each lfp trace to a region in 'regions' by index for visualization
    end
    out.LFPtimes = LFPtimes;
    out.eventStartIndices = eventStartIndices;
    out.eventEndIndices = eventEndIndices;
    out.windowStartIndices = windowStartIndices;
    out.windowEndIndices = windowEndIndices;
    out.index = index;

end