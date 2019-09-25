function out = calcepochreplayspikestats(index, placefieldindex, placefieldfilter, decodefilter, sigthresh)
%  out = calcepochreplayspikestats(index, placefieldindex, placefieldindex, %  placefieldplacefieldfilter, decodefilter, sigthresh)
%
% This function goes through each of the significant events in the decodefilter
% as indicated by the p value in replaystat field in eventdata.  sigthresh is
% the minimum value to be considered significant.  This function then 
% determines, for each spike, the peak rate of that cell's place field and the 
% spatial rate of that cell at the location where the replay event occurred.
% 
% This function must be run after calcepochreplaystats which should set the
% ripplestat field in the decode filter.  Also, the placefieldindex and
% placefieldfilter are not the same as the standard placefield filters, as they
% have firing rate information for all cells; not just cells with place fields.


animalnum = index(1);
epochnum = index(2);
trajmapping = [1 1 2 2];
binsize = .015; %default temporal bin
nevents = length(decodefilter(animalnum).output{1}(epochnum).eventdata);
out = {};
outindex = 1;

for eventindex = 1:length(decodefilter(animalnum).output{1}(epochnum).eventdata)
    df = decodefilter(animalnum).output{1}(epochnum);
    pf = placefieldfilter(animalnum).output{placefieldindex}(epochnum);
    disp(eventindex)
    placefielddata = [];
    spikedata = [];
    decodedata = [];
    indexlist = [];

    %pick out all the matching cells from the placefield data and the
    %decoding data
    %placefieldfilter contains linear rates, and is n by x, where n is the
    %number of cells and x is the number of spatial bins
    %spikedata contains spikecounts, and is n by t, where t is the
    %number of temporal bins in the data to be decoded.
    %matches = rowfind(placefieldfilter(animalnum).output{placefieldindex}(epochnum).index(:,[1 3 4]),decodefilter(animalnum).output{1}(epochnum).index(:,[1 3 4])); %find the matching cell indices

    % check to see if this was a significant event
    rstat = df.eventdata(eventindex).replaystat;
    if ((length(rstat) == 3) & (df.eventtraj(eventindex) > 0) & ...
	(df.eventdata(eventindex).replaystat(3) < sigthresh))
	% find, for each cell in the decodefilter, the index in the place field
	% filter
	matches = rowfind(df.index(:,[1 3 4]), pf.index(:,[1 3 4])); 
	% check to see if any of the matches are 0, if so there is an error
	if (sum(matches == 0))
	    error('all cells should be represented in the place field filter');
	end

	% we need to be able to look up the firing rate at the animal's
	% location, so we first figure out the boundaries for the trajectories
	trajbound = [find(pf.dist == 0) length(pf.dist)];

	nspikes = length(df.eventdata(eventindex).cellindex);
	out{outindex} = zeros(nspikes,2);

	% go through each of the cells in order within the event get the peak 
 	% rate and rate at the animals location for that place field
	for c = 1:nspikes
	    pfcell = matches(df.eventdata(eventindex).cellindex(c));
	    peakrate = max(pf.rates(pfcell,:));

	    % get the rate on the appropriate trajectory
	    traj = df.eventtraj(eventindex);
	    ind = lookup(df.eventdist(eventindex), pf.dist(trajbound(traj):(trajbound(traj+1)-1)));
	    localrate = pf.rates(pfcell, ind+trajbound(traj)-1);
	    out{outindex}(c,1) = peakrate;
	    out{outindex}(c,2) = localrate;
	    out{outindex}(c,3) = df.eventdata(eventindex).spiketimes(c);
	end
	outindex = outindex + 1;
    end
end
