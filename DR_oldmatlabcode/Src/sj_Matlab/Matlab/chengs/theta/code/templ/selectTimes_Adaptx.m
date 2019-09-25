function selectTimes_Adaptx 

global fmaux
epoch= 2; 

load /bach/Ter/data/Adaptx3.0t0.05d/teradapt.mat

timeselect= [];
traj= 0;
global behavdata
for n=1:size(teradapt(epoch).cellnum,1)
    d= teradapt(epoch).cellnum(n,1);
    e= teradapt(epoch).cellnum(n,2);
    t= teradapt(epoch).cellnum(n,3);
    c= teradapt(epoch).cellnum(n,4);
    if d== -1; continue; end
    
    if (e ~= epoch); error('selection gone wrong'); end;
    
    if teradapt(epoch).traj(n,2) ~= teradapt(epoch).armnum(n)
	continue;
    end
    
    loadVar(fmaux.data2dir, 'behavdata', d);
    timeselect{d}{e}{traj+1}.time=  teradapt(epoch).passtime{n};
    timeselect{d}{e}{traj+1}.index= ...
	findTimeIndex(timeselect{d}{e}{traj+1}.time,behavdata{d}{e}.time);    
    
    traj= mod(traj+1,4);
end

save /bach/Ter/data2/timeselect-passAdaptx3-0t0-05d timeselect