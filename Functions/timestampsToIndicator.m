function indicator = timestampsToIndicator(ts, time)
% make binary time series of from time stamps
time = time(:);
ts = ts(:);
tsidx = lookup(ts, time);
indicator = zeros(length(time),1);
indicator(tsidx)=1;
end