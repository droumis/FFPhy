function statsum= extractStats(varname, indices, minspikes)
%function statsum= extractStats(varname, indices, minspikes)
% 
%  Extract summary of statistic for one animal.
%  Each statistic is usually computed for each cell at several time points and
%  saved per cell, i.e. stat{d}{e}{t}{c}.xyz= [1 2 3].
%  This function extract the statistics across all cells at particular time
%  indices. 
%  
% Output:
%  statsum{i}   stat for each element of indices
%
% Input:
%  prefix    animals for which to collect statistics
%  varname   name of statistics
%  index:            index of time selection (e.g. pass, occupancy time)


global fmaux select tetloc

fname= ['stats-' fmaux.selectid '.mat'];
load(fname);
load(fmaux.select);
loadVar(fmaux.datadir, 'tetloc', 0, fmaux.prefix, 1);
nindices= length(indices);
n= zeros(nindices,1);
statsum= cell(nindices, 1);

resetCellList;
while 1
    [d,e,t,c]= getNextCell;
    if isempty(d) | isempty(e) | isempty(t) | isempty(c); break; end

    % check whether statistics is valid variable
    stmp= stats{d}{e}{t}{c};
    if isempty(stmp); continue; end
%    if ~isfield(stmp, varname)
%        error(['stats doesn''t have field named ' varname]);
%    end
    svar= stmp.(varname);
    if isempty(svar); continue; end
    for i= 1:length(indices)
        for k= 1:length(svar)
            if length(svar{k}) < indices(i); continue; end
            if ~isfinite(svar{k}(indices(i))); continue; end
            if nargin >= 3 & stats{d}{e}{t}{c}.nspikes{k}(indices(i)) < minspikes; continue; end
            % select trajectories
%            traj= select.a{fmaux.currentCell}{k}.traj;
%            if(traj~=3); continue; end

            % statistic was valid add to list	
            n(i)= n(i)+1;
            statsum{i}.cellnum(n(i),:)= [d e t c];
            statsum{i}.aindex(n(i))= k;
            statsum{i}.val(n(i))= svar{k}(indices(i));
        end
    end
end
