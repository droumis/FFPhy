function out = calclinadapt(sind, tind, excludeperiods, spikes, linpos, theta,varargin)


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
cpxspacing = 5;	% approximately 5 cm spacing

data = createphaseprecessadaptdata(sind, tind, excludeperiods, linpos, spikes, theta, timestep);


if isempty(filter)
    filter.name= 'AscentFilter';
    filter.alternatePass= 1;
    % learning 
    filter.eps= [2 * ones(1,16), 0.1 * ones(1,4)];
    filter.niter= 20;
    filter.maxGradient=100;
end


small=1e-4;
if isempty(model)
    model.spatial.conv= 0.005;
    model.spatial.convper= .1;
    model.spatial.outputCompress= 1;
    model.spatial.outputInterval= 1000; % timesteps!

    model.name= 'PosPhase_Isi';

    model.spatial.name='CSplines2d';
    model.spatial.periodic= 'y';
    model.spatial.operator.name= 'Rectify';
    % figure out how many control points there should be for each trajectory
    ncpy = 12;
    for k = 1:(max(data.traj)+1)
	t = find(data.traj == (k-1));
	minpos=min(data.linpos(t))-small;
	maxpos=max(data.linpos(t))+small;
	ncpx = round((maxpos - minpos) / 5);
	model.spatial.cpx{k}= [minpos-small, linspace(minpos,maxpos,ncpx), maxpos+small]; 
	model.spatial.cpy{k}= 2*pi*[0:ncpy]/ncpy - pi;
    end

    model.max_isi= 0.08; %in sec
    model.isi.name='CSplines'; 
    model.isi.operator.name= 'Rectify';
    %model.isi.name='CSplinesNorm'; 
    %model.isi.cpx=[-small 0 1:4:25 30:10:70 80+small 81]/1000;
    model.isi.cpx=[-small 1 3 5 7 11 15:10:65 80+small 81]/1000;
    model.isi.conv= .005;
    model.isi.convper= .1;
    model.isi.outputCompress= 1;
    model.isi.outputInterval= 1000; % timesteps!
end

out= adaptFilter(data,model,filter);
out.spikeindex = sind;
out.thetaindex = tind;

