function totaltimes = calcSegmentTimeOn(segments, excludeperiods)
%totaltimes = calcSegmentTimeOn(segments, excludeperiods)
%Calculates the total time within SEGMENTS vector that is not
%excluded by the list of exclusion periods [startexclude endexclude];


for i = 1:length(segments)-1
    
    segment = [segments(i) segments(i+1)];
    totaltime = segment(2)-segment(1);
    if ~isempty(excludeperiods)
        containedperiods = find( (excludeperiods(:,1) >= segment(1)) & (excludeperiods(:,2) <= segment(2)) );
        totaltime = totaltime - sum(excludeperiods(containedperiods,2) - excludeperiods(containedperiods,1));

        openstartperiods = find( (excludeperiods(:,1) < segment(1)) & (excludeperiods(:,2) <= segment(2)) & (excludeperiods(:,2) >= segment(1)) );
        totaltime = totaltime - sum(excludeperiods(openstartperiods,2) - segment(1));

        openendperiods = find( (excludeperiods(:,1) >= segment(1)) & (excludeperiods(:,1) <= segment(2)) & (excludeperiods(:,2) > segment(2)) );
        totaltime = totaltime - sum(segment(2) - excludeperiods(openendperiods,1));

        %If any exclusiopn periods encompass the segment, the total time is 0
        completeperiods = find( (excludeperiods(:,1) <= segment(1)) & (excludeperiods(:,2) >= segment(2)) );
        if ~isempty(completeperiods)
            totaltime = 0;
        end
    end
    totaltimes(i) = totaltime;
end