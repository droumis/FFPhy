function out = calcripratehist(rip, excludeperiods, numtetrodes,binsize)
%out = calcriprate(rip, excludetimes, numtetrodes)
%
%   rip.times: time of each sample in nripples, 1msec bins
%   rip.nripples: number of electrodes with a ripple recorded at that time
%   point
%
%   out is [ rate proportiontime] 
%       proprotiontime is proportion of included time during which ripples
%   were recorded
%       rate is number ripples/sec during included time

%apply excludetimes to nripples 
includetimes = ~isExcluded(rip.times, excludeperiods); %list of ones and zeros sampled every millisecond, ones = included, zeros = excluded
if size(rip.nripples) ~= size(includetimes)
    includetimes = includetimes';
end
includerips = rip.nripples .* includetimes;

%calculate proportion of time spent in ripples
%riptimecount = length( find(includerips >=numtetrodes) );
%rippercenttime = riptimecount / sum(includetimes); %proportion total recorded time that included ripples

%calculate number of ripples per time 
a = zeros(1, length(includerips));
a(find(includerips >= numtetrodes)) = 1;
ripcount = rip.times(diff(a)==1);
% calculate ripples per bin
%[ripcount xout] = hist(ripcount, [rip.times(1):binsize:rip.times(length(rip.times))]);
[ripcount xout] = hist(ripcount, [1329:binsize:rip.times(length(rip.times))]);
% calculate included time per bin
%[goodtimes xout] = hist(rip.times(includetimes), [rip.times(1):binsize:rip.times(length(rip.times))]);
[goodtimes xout] = hist(rip.times(includetimes), [1329:binsize:rip.times(length(rip.times))]);


%calculate rate, divide by 1000 change from msec to seconds
riprate = ripcount./(goodtimes/1000);

out.goodtimes = goodtimes/1000;
out.ripcount = ripcount;
out.riprate = riprate;
out.xout = xout;
end