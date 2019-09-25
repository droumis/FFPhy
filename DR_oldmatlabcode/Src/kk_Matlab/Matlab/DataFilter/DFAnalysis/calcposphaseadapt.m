function out = calclinadapt(index, excludeperiods, spikes, linpos, varargin)


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
data = createadaptdata(index, excludeperiods, linpos, spikes, timestep);

seginfo = linpos{index(1)}{index(2)}.segmentInfo;


if isempty(filter)
    filter.name= 'AscentFilter';
    filter.maxGradient=50;
    % the learning rate (eps) has a single set of four values for the spatial
    % component and another set of four values for the temporal component
    filter.eps= [2 * ones(1,4), 0.05 * ones(1,4)];
    filter.alternatePass= 1;
    filter.niter= 20;
end


if isempty(model)
    small=1e-4;
    minpos=min(data.linpos)-small;
    maxpos=max(data.linpos)+small;
    ncpx= floor((maxpos - minpos) / 5); % spacing at 5 cm

    model.name= 'Pos_Isi';
    model.max_isi= 0.08; %in sec
    model.x.name='CSplines';
    % we need one set of control point for each of the four trajectories 
    model.x.cpx{1} = [minpos-small, linspace(minpos,maxpos,ncpx), maxpos+small];
    model.x.cpx{2} = [minpos-small, linspace(minpos,maxpos,ncpx), maxpos+small];
    model.x.cpx{3} = [minpos-small, linspace(minpos,maxpos,ncpx), maxpos+small];
    model.x.cpx{4} = [minpos-small, linspace(minpos,maxpos,ncpx), maxpos+small];
    model.x.conv= 0.1;
    model.x.convper= .1;
    model.x.outputCompress= 1;
    model.x.outputInterval= 1000; % timesteps!

    % set up the assocation between control points.  To do so we get the
    % segment Length for the first segment of the track corresponding to the
    % home arm. We then set the control points associated with that home arm
    % to be the same
    homelen = seginfo.segmentLength(1);
    model.x.map= [1 minpos homelen 3; 3 minpos homelen 1; ...
    		  0 minpos homelen 2;  2 minpos homelen 0]';


    model.isi.name='CSplines';
    %model.isi.name='CSplinesNorm';
    %model.isi.cpx=[-small 0 1:2:25 30:5:45  50:10:70 80+small 81]/1000;
    model.isi.cpx=[-small 0 1:4:25 30:10:70 80+small 81]/1000;
    model.isi.conv= .1;
    model.isi.convper= .1;
    model.isi.outputCompress= 1;
    model.isi.outputInterval= 1000; % timesteps!
end

out = adaptFilter(data,model,filter);
out.filter = filter;
out.index = index



