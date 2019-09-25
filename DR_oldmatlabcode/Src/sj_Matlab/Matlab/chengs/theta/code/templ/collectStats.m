function statsum= collectStats(varname, index)

minspikes=0;

global fmaux

fname= ['stats-' fmaux.celllist '-' fmaux.timelist '.mat'];
load(fname);

statsum= [];
statsum.name= varname;
statsum.index= index;

n=0;
[d,e,t,c]= startCellList;
while ~isempty(d)
    for traj=1:4
	% check whether statistics was valid variable
	stmp= stats{d}{e}{t}{c}{traj};
	if isempty(stmp); continue; end
	eval(['svar= stmp.' varname ';']);

	if length(svar) < index; continue; end
	if isnan(svar(index)); continue; end
	if stats{d}{e}{t}{c}{traj}.nspikes(index) < minspikes;
	    continue; 
	end
	
% statistic was valid add to list	
	n= n+1;
	statsum.cellnum(n,:)= [d e t c];
	statsum.traj(n)= traj;
	statsum.val(n)= svar(index);
    end
    [d,e,t,c]= getNextCell;
end

