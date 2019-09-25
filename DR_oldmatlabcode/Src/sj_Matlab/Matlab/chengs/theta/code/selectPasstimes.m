function selectTimes
% figure out times, at which one pass ends

global fmaux adapt select behavdata
epoch= 2; 

loadFile(fmaux.select);
loadVar(fmaux.adaptdir, 'adapt', 0, fmaux.prefix, true);

traj= 0;
for n=1:size(adapt(epoch).cellnum,1)
    d= adapt(epoch).cellnum(n,1);
    e= adapt(epoch).cellnum(n,2);

    % skip invalid cells
    if d== -1; continue; end
    
    if (e ~= epoch); error('selection gone wrong'); end;
    
    % adapt lists 8 entries per cell: 4 trajectories X 2 arms
    % we are only interested in the ends of the 4 trajs
    if adapt(epoch).traj(n,2) ~= adapt(epoch).armnum(n)
        continue;
    end

    selid= getSelectId(adapt(epoch).cellnum(n,:));
    if isempty(selid); continue; end

    select.a{selid}{traj+1}.traj= traj;
    select.x{selid}{traj+1}.time=  adapt(epoch).passtime{n};

%    loadVar(fmaux.data2dir, 'behavdata', d);
%    select.x{selid}{traj+1}.timeindex=  ...
%    findIndex(adapt(epoch).passtime{n},behavdata{d}{e}.time);    
    
    traj= mod(traj+1,4);
end

save([fmaux.select '-passtimes'], 'select');
