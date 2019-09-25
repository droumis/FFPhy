function tmp
% write files with behavioral states for Murat Okatan
% 0,1,2,3:  valid trajectories and running
% -1:       not running
% -2:       ripples, but running
% -3:       ripples and not running



global behavdata 
setRoot

[tet, epochs]= collectTet('all');
ne= length(epochs.rat);

for ie=1:ne
    rat= epochs.rat{ie}; d= epochs.num(ie,1); e= epochs.num(ie,2);

    data2dir= fullfile(root,rat,'data2');
    loadVar(data2dir, 'behavdata', d);
    bd= behavdata{d}{e};

    data{d}{e}.time= bd.time;

    % behavioral state
    bd.traj(bd.ripple==1 & bd.traj>=0)= -2;
    bd.traj(bd.ripple==1 & bd.traj==-1)= -3;
    data{d}{e}.traj= bd.traj;

    % theta phase angle
    data{d}{e}.theta= bd.phase;

    if ie==ne | ~strcmp(rat, epochs.rat{ie+1}) | d~= epochs.num(ie+1,1)
        save(sprintf('%sdata%.2d', rat, d), 'data');
        data= {};
    end

end
