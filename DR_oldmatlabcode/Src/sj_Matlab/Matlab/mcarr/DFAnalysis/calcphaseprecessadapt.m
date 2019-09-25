function out = calcphaseprecessadapt(index, excludeperiods, spikes, linpos, varargin)


model = [];
filter = [];
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'model'
                model = varargin{option+1};
            case 'filter'
                filter = varargin{option+1};    
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

% to fill in the data matrix we need to create a 2 ms time array and
% interpolate the time, lindist and traj elements of the linpos statematrix

timestep = 0.002; % 2ms timestep
time = linpos{index(1)}{index(2)}.statematrix.time;
data.time = time(1):timestep:time(end);

lp = linpos{index(1)}{index(2)}.statematrix.lindist;
data.linpos = interp1(time, lp, data.time, 'linear');

traj = linpos{index(1)}{index(2)}.statematrix.traj;
data.traj = interp1(time, traj, data.time, 'nearest');


% note that the trajectory definitions from the statematrix use -1 for invalid
% elements and 1-n for the valid elements.  The adaptive filtering code expects
% the valid elements to start at 0, so we subtract 1 from all values > 0
tmp = find(data.traj > 0);
data.traj(tmp) = data.traj(tmp) - 1;
data.traj(find(isExcluded(data.time, excludeperiods))) = -1;
data.spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);

if isempty(filter)
    filter.name= 'AscentFilter';
    filter.maxGradient=50;
    % the learning rate (eps) has a single set of four values for the spatial
    % component and another set of four values for the temporal component
    filter.eps= [0.15 * ones(1,16), 0.05 * ones(1,4)];
    %filter.eps= [2 * ones(1,4), 0.05 * ones(1,4)];
    filter.alternatePass= 1;
    filter.niter= 20;
end


if isempty(model)
    small=1e-4;
    minpos=min(data.linpos)-small;
    maxpos=max(data.linpos)+small;
    ncpx= 20; % spacing n=20
    ncpy= 12; % separation 30deg

    model.name= 'PosPhase_Isi';
    model.max_isi= 0.08; %in sec


    model.xp.name='CSplines2d';
    model.xp.periodic= 'y';
    model.xp.cpx{1}= [minpos-small, linspace(minpos,maxpos,ncpx), maxpos+small];
    model.xp.cpx{2}= [minpos-small, linspace(minpos,maxpos,ncpx), maxpos+small];
    model.xp.cpx{3}= [minpos-small, linspace(minpos,maxpos,ncpx), maxpos+small];
    model.xp.cpx{4}= [minpos-small, linspace(minpos,maxpos,ncpx), maxpos+small];
    model.xp.cpy{1}= 2*pi*[0:ncpy]/ncpy;
    model.xp.cpy{2}= 2*pi*[0:ncpy]/ncpy;
    model.xp.cpy{3}= 2*pi*[0:ncpy]/ncpy;
    model.xp.cpy{4}= 2*pi*[0:ncpy]/ncpy;
    model.xp.conv= 0.3;
    model.xp.convper= .1;
    % we need one set of control point for each of the four trajectories 
    model.x.outputCompress= 1;
    model.x.outputInterval= 1000; % timesteps!

    model.isi.name='CSplines';
    %model.isi.name='CSplinesNorm';
    %model.isi.cpx=[-small 0 1:2:25 30:5:45  50:10:70 80+small 81]/1000;
    model.isi.cpx=[-small 0 1:4:25 30:10:70 80+small 81]/1000;
    model.isi.conv= .1;
    model.isi.convper= .1;
    model.isi.outputCompress= 1;
    model.isi.outputInterval= 1000; % timesteps!
end

out= adaptFilter(data,model,filter);
out.index = index;



