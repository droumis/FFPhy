function convertData

%  Timesteps are defined as valid, if the animal is running at instantaneouly
%  velocity of 4cm/s or more.


minvel= 4; % [cm/s]
global fmaux 
global task lindistspikepos lindistpos pos vel

loadVar(fmaux.datadir, 'task', 0, fmaux.prefix, 1);

oldd= -123;
olde= -123;

[d,e,t,c]= startCellList;
%load(fmaux.EEGselect);
while ~isempty(d)
    if oldd~= d | olde~= e 
        % get behavioral data
        loadVar(fmaux.datadir, 'lindistpos', d, fmaux.prefix, 1);
        loadVar(fmaux.datadir, 'lindistspikepos', d, fmaux.prefix, 1);
        loadVar(fmaux.datadir, 'pos', d, fmaux.prefix, 1);
        loadVar(fmaux.datadir, 'vel', d, fmaux.prefix, 1);
        
        times= lindistpos{d}{e}.data(:,1);
        behavdata{d}{e}.time= times;
        behavdata{d}{e}.linpos=lindistpos{d}{e}.data(:,2) ...
            * task{d}{e}.pixelsize;
        
        ts= pos{d}{e}.data(:,1);
        xs= pos{d}{e}.data(:,2)* task{d}{e}.pixelsize;
        ys= pos{d}{e}.data(:,3)* task{d}{e}.pixelsize;
        dirs= pos{d}{e}.data(:,4);

        behavdata{d}{e}.xpos= interp1(ts, xs, times, 'splines');
        behavdata{d}{e}.ypos= interp1(ts, ys, times, 'splines');
        behavdata{d}{e}.dir=  interp1(ts, dirs, times, 'splines');

%        behavdata{d}{e}.traj=lindistpos{d}{e}.estinfothetavel;
        behavdata{d}{e}.traj=lindistpos{d}{e}.estinfo;
        v= interp1(vel{d}{e}.data(:,1), vel{d}{e}.data(:,2), times);
        behavdata{d}{e}.traj(find(v<minvel))= -1;

    % get theta phase
%        EEGfile= sprintf('%s/%stheta%.2d-%d-%.2d',fmaux.EEGdir,fmaux.prefix,...
%                 d,e,EEGselect{d}{e});
%        load(EEGfile);
%        eval(['eeg = ' fmaux.prefix 'theta{d}{e}{EEGselect{d}{e}};']);
%        behavdata{d}{e}.phase= convertEEG2Phase(eeg, times);
    end
    
    % spike data
    sptime= lindistspikepos{d}{e}{t}{c}.data(:,1);
    valid= find(sptime > times(1) & sptime < times(end));
    nvalid= length(valid);
    spikedata{d}{e}{t}{c}.time= sptime(valid);
    % determine time index of spike
    ind= findIndex(spikedata{d}{e}{t}{c}.time, times);
    %@@ consistency check
    dt= mean(diff(times))+1e-15;
    if(sum((sptime(valid)-times(ind))<=dt) ~= nvalid);
        error('spike indices assigned incorrectly');
    end
    spikedata{d}{e}{t}{c}.index= ind;
    spikedata{d}{e}{t}{c}.linpos= behavdata{d}{e}.linpos(ind);
    spikedata{d}{e}{t}{c}.traj= behavdata{d}{e}.traj(ind);
%    spikedata{d}{e}{t}{c}.phase = behavdata{d}{e}.phase(ind);
%    spikedata{d}{e}{t}{c}.dir= lindistspikepos{d}{e}{t}{c}.data(valid,3);
    
    oldd= d;
    olde= e;

    [d,e,t,c]= getNextCell;

    if isempty(d) | (oldd~= d & oldd >= 1)
        dstring= sprintf('%.2d', oldd);
        save(['behavdata' dstring], 'behavdata');
        save(['spikedata' dstring], 'spikedata');
        clear behavdata spikedata
    end
end

