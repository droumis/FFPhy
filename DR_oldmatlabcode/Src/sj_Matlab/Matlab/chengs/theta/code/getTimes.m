function tm= getTimes(timeid, rat, d, e, traj, maxt)
%function tm= getTimes(timeid, rat, d, e, traj, maxt)
%
% Return times at which certain behavioral landmarks are reached using saved
% table if available.
% timeid: passes occ minocc acc minacc
%
%function tm= getTimes(timeid, maxt)
% Calculate and return the times for all epochs and trajectories.
%
% 

dt= 0.002;
debugging= 1;

if nargin==2; maxt= rat; end
if nargin==1 | nargin<6; 
    nomaxt= true; 
else
    nomaxt= false; 
end

setRoot

fname= sprintf('%s/data/times-%s.mat', root, timeid);

if exist(fname, 'file') & nargin>2
    load(fname);
else
    time_based= true; times=[];

    switch timeid
    case {'minocc'}
        timeinterval= 60;
    case {'occ'}
        timeinterval= 1;
    case {'pass'}
        time_based= false;
    otherwise
        error(['unknown timeid ' timeid]);
    end

    global spikedata behavdata info lindistpos  task
    [tet, epochs]= collectTet('all');
    ne= length(epochs.rat);

    for ie=1:ne
        rat= epochs.rat{ie}; d= epochs.num(ie,1); e= epochs.num(ie,2);
    
        datadir= fullfile(root,rat,'data');
        data2dir= fullfile(root,rat,'data2');
        loadVar(datadir, 'task', 0, rat, 1);
        loadVar(data2dir, 'info', 0);

        if time_based
            minpos= info{d}{e}.centerlinpos/task{d}{e}.pixelsize;
            maxpos= info{d}{e}.maxlinpos/task{d}{e}.pixelsize;

            loadVar(datadir, 'lindistpos', d, rat, 1);

            for j=0:3
                if strfind(timeid, 'occ')
                    ontraj= find(minpos<=lindistpos{d}{e}.data(:,2) & lindistpos{d}{e}.estinfo==j);
                    ntraj= length(ontraj);
                    if nomaxt
                        maxt= floor(ntraj*dt/timeinterval);
                    end
                    len= round(timeinterval/dt);
                    ind=[1:maxt]'*len;
                    valid= ind<=ntraj;
                    times.(rat){d}{e}{j+1}= nan*ones(maxt,1);
                    times.(rat){d}{e}{j+1}(valid)= lindistpos{d}{e}.data(ontraj(ind(valid)),1);
                end
            end
        else % passes
            % ordering of arms in adaptest, if task is 1-3-7
            % id=1: *1 3, traj 0
            % id=2: 1 *3, traj 0
            % id=3: 3 *1, traj 1
            % id=4: *3 1, traj 1
            % id=5: *1 7, traj 2
            % id=6: 1 *7, traj 2
            % id=7: 7 *1, traj 3
            % id=8: *7 1, traj 3
            %
            % arms of 8-arm maze are numbered clockwise, the home arm being 1

            global adaptest
            loadVar([datadir '/Adaptx3.0t0.05d'], 'adaptest', d, rat, 1);
            if e==2; ie= 1; else ie= 2; end

            % find active tetrode  (any tetrode)
            for t=1:length(adaptest{d}{e})
                if ~isempty(adaptest{d}{e}{t}); break; end
            end
            % find active cell (any cell)
            for c=1:length(adaptest{d}{e}{t})
                if ~isempty(adaptest{d}{e}{t}{c}); break; end
            end
            sstat= adaptest{d}{e}{t}{c}.sstat;
            st= sstat(1).time; % stat times

            % consistency check
            if debugging
                for i=2:8; ti{i}= sstat(i).time; end
                for i=2:8; 
                    if sum(ti{i}==st) ~= length(st); 
                        error('sstat(i).times not identical for all arms');
                    end
                end
            end

            % get end times of passes and occupancies of arms
            for itraj=0:3  % end of traj in OUTER arms
                times.(rat){d}{e}{itraj+1}= st(sstat(2*(itraj+1)).trajendind)';
            end

        end % if time_based
    end % epoch
    save(fname, 'times');
end

if nargin > 2
    tm= times.(rat){d}{e}{traj+1};
    if size(tm,2)>1 & size(tm,1)==1; tm= tm'; end
    if ~nomaxt
        len= length(tm);
        if len<=maxt; 
            tm(len+1:maxt,:)= nan*ones(maxt-len,1);
        end
        tm= tm(1:maxt); 
    end
end
