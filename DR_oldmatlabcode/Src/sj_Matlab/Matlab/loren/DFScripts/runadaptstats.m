


load /data/mkarlsso/Adaptive/ca1trackA.mat
ca1f = setfilterfunction(ca1f, 'calclinadaptstats', {'spikes', 'linpos', 'filteroutput'});
ca1trackAstats = runfilter(ca1f);

load /data/mkarlsso/Adaptive/ca1trackB.mat
ca1f = setfilterfunction(ca1f, 'calclinadaptstats', {'spikes', 'linpos', 'filteroutput'});
ca1trackBstats = runfilter(ca1f);


load /data/mkarlsso/Adaptive/ca3trackA.mat
ca3f = setfilterfunction(ca3f, 'calclinadaptstats', {'spikes', 'linpos', 'filteroutput'});
ca3trackAstats = runfilter(ca3f);

load /data/mkarlsso/Adaptive/ca3trackB.mat
ca3f = setfilterfunction(ca3f, 'calclinadaptstats', {'spikes', 'linpos', 'filteroutput'});
ca3trackBstats = runfilter(ca3f);


