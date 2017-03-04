function [out] = DR_timestampskipcheck(timestamps)

%zero if there are no gaps in the timestampds out to the 0.01ms
anyskipped = sum(diff(round(round(diff(round(timestamps,7)), 7),6)));

out = anyskipped;
