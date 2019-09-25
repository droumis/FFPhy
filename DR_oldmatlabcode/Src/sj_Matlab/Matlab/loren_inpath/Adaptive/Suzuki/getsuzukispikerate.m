function [rate] = getspikerate(cp, heights, starttime, endtime)
% area = GETSPIKERATE(cp, heights, starttime, endtime)
%        Computes the firing rate in an interval, given the following arguments:
%        cp         locations of control points
%        heights    corresponding heights of control points
%        starttime  the start time
%        end        the ending time
%         
%        Note that the area is unnormalized.
i1 = max(find(cp == starttime));
i2 = min(find(cp == endtime));

if (isempty(i1) | isempty(i2))
    error('starttime and endtime must fall on control points');
end

spacing = cp(i1+1) - cp(i1);
area = 0;
for i = i1:(i2-1),
    area = area + (-heights(i-1)/24 + heights(i) * 13/24 + ...
                + heights(i+1) * 13/24 - heights(i+2)/24) * spacing;
end 

% to convert to firing rate, we need to divide by the total number of ms in the interval
rate = area / (cp(i2) - cp(i1));
