function statsum= collectStatsSeq(varname, maxindex)

minspikes=10;

global fmaux

fname= ['stats-' fmaux.celllist '-' fmaux.timelist '.mat'];
load(fname);

statsum= [];
statsum.name= varname;

n=0;
[d,e,t,c]= startCellList;
while ~isempty(d)
    for traj=1:4
	% check whether statistics was valid variable
	stmp= stats{d}{e}{t}{c}{traj};
	if isempty(stmp); continue; end
	eval(['svar= stmp.' varname ';']);

	if length(svar) < maxindex; continue; end
	if find(isnan(svar(1:maxindex))); continue; end

	if sum(stats{d}{e}{t}{c}{traj}.nspikes(1:maxindex)>=minspikes)~= maxindex; continue; end

% statistic was valid add to list	
	if size(svar,2)==1; svar= svar'; end
	n= n+1;
	statsum.cellnum(n,:)= [d e t c];
	statsum.traj(n)= traj;
	statsum.val(n,:)= svar(1:maxindex);
    end
    [d,e,t,c]= getNextCell;
end

