function out = getrawdata_DR(index, excludetimes, linpos, spikes, includestates, minV)
%trying to get all the data for the filtered cells so that that I can
%combine them within a day


%things I need
out.spikes = spikes;
out.state = excludetimes;
out.linpos = linpos;
out.index = index;
out.includestates = includestates;
out.minV = minV;

