function out = sj_isTimePtExcluded(timepoints, excludePeriods)
% Shantanu - Jun2011

% isExcluded in DataFilter framework takes in a timerange vector in ms
% resolution and a list of excludedPeriods with start and end times (nX2
% size), and returns a out vector of length times with 1s and 0s based on
% exclude Periods

% sj_isTimePtExcluded takes in a list of timepoints to check and exclude periods
% It checks whether each timepoint in the input list belongs to any of the
% excludePeriods and returns a list of 1s and zeros. 1 means that
% timepoint belongs to exclude time and should be discarded, just like
% isExcluded

% Used to evaluate DIO stimtime for sj_calcriprate, which is called by
% DFAsj_getriprate_noDIO.m  I need to know which DIOs to keep to subtract
% from ripcount to accurate riprate. This is different from a DIOtimefilter
% where I exclude time period of entire DIO stimulation artifact, which is
% not the right/accurate thing to do for stimrate. More accurate is to just
% subtract the nDIOs from nRipples

% out = isTimePtExcluded(timepoints, excludePeriods)
% out = isExcluded(times, excludePeriods)
% Returns a list of 1s and 0s signifying inclusion or exclusion, as defined
% by the exclude start and end times in excludePeriods.  ExcludePeriods is
% an n by 2 matrix where the columns are the start and end times for each
% exclusion period.
% 1 signifies times that fall in between start and end times of
% excludePeriods

out = zeros(size(timepoints));

if ~isempty(excludePeriods)
    
    for pt = 1:length(timepoints),
        currtimept = timepoints(pt);   
        for n=1:size(excludePeriods,1) 
            if currtimept>=excludePeriods(n,1) & currtimept<=excludePeriods(n,2) 
                out(pt)=1; % If hit, exit the excludePeriods loop
                break
            end
        end
    end
    
end  
        
