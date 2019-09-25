function [avg,stdev] = stat(x)
           %STAT Interesting statistics.
           n = length(x);
           avg = mean(x);
           %stdev = sqrt(sum((x-mean(x,n)).^2)/n);
           stdev=std(x);
