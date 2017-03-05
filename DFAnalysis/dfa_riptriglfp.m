function out = dfa_riptriglfp(index, excludeperiods, eeg, ripple, events, varargin)

% gathers lfp around ripple times

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
wininds = win*samprate;
endtime = (num_samp/samprate) + starttime;
disp(sprintf('endtime_eeg: %.04f calculated_endtime: %.04f', endtime_eeg/60, endtime/60));

%Remove triggering events that are too close to the beginning or end
while eventtimes(1)<(starttime+win(1))
    eventtimes(1) = [];
end
while eventtimes(end)>(endtime-win(2))
    eventtimes(end) = [];
end

    LFPtimes=(starttime:1/samprate:endtime)'; % prob should be using the adjusted timestamps instead, no?
    eventStartIndices = lookup(eventtimes(:,1),LFPtimes);
    eventEndIndices = lookup(eventtimes(:,2),LFPtimes);
%     windowStartIndices = lookup(eventtimes(:,1)-win(1),LFPtimes);
%     windowEndIndices = lookup(eventtimes(:,2)+win(2),LFPtimes);
%     eventStartLFPtimes = LFPtimes(eventStartIndices);
%     windowStartEndTimes = [eventStartLFPtimes(:)-win(1) eventStartLFPtimes(:)+win(2)]; 

    %% STEP 3: Gather the ripple window data from all the regions 
    for currrip=1:length(eventStartIndices) %for each ripple from source region
        clear YripLFPdata
        %         for iArea = 1:length(regions);
        for iNTrode = 1:length(index(:,3)); %for each ntrode
%             nTrodeID = index(iNTrode,3);
            % Gather lfp data for each ntrode within the current rip window
            out.data{currrip}(iNTrode,:) = eeg{index(iNTrode,1)}{index(iNTrode,2)}{index(iNTrode,3)}.data(eventStartIndices(currrip)-wininds(1):eventStartIndices(currrip)+wininds(2));
        end
%         YripLFPdataMAT{iLFPtype}{currrip} = cell2mat(YripLFPdata); %stack all the traces from all the regions next to each other
%         Xwindowtimes{iLFPtype}{currrip} = WindowStartEndTimes(currrip,1):1/samprate:WindowStartEndTimes(currrip,2);
%         lfptraceLUTregion = cell2mat(lfptraceregion{iLFPtype}'); % [1 x length(nTrodes)] vec mapping each lfp trace to a region in 'regions' by index for visualization
    end
    out.LFPtimes = LFPtimes;
    out.eventStartIndices = eventStartIndices;
    out.eventEndIndices = eventEndIndices;
    out.win = win;
%     out.windowStartIndices = windowStartIndices;
%     out.windowEndIndices = windowEndIndices;
    out.index = index;

end