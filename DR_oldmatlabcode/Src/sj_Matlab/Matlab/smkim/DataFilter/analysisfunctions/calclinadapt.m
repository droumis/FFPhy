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


data.time = linpos{index(1)}{index(2)}.statematrix.time;
data.linpos = linpos{index(1)}{index(2)}.statematrix.lindist;
data.traj = linpos{index(1)}{index(2)}.statematrix.traj    
data.traj(find(isExcluded(data.time, excludeperiods))) = -1;
data.spiketime = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);


if isempty(filter)
    filter.name= 'AscentFilter';
    filter.maxGradient=50;
    filter.eps= [.15* ones(1,16), 0.05 * ones(1,4)];
    filter.alternatePass= 1;
    filter.niter= 20;
end


if isempty(model)
    small=1e-4;
    minpos=min(data.linpos)-small;
    maxpos=max(data.linpos)+small;
    ncpx= 20; % spacing n=20

    model.name= 'Pos_Isi';
    model.max_isi= 0.08; %in sec
    model.x.name='CSplines';
    model.x.cpx= [minpos-small, linspace(minpos,maxpos,ncpx), maxpos+small];
    model.x.conv= 0.1;
    model.x.convper= .1;
    model.x.outputMode= 'interval'; % other option is 'timesteps' or not setting
    model.x.outputInterval= 2;  % [sec]

    model.isi.name='CSplines';
    %model.isi.name='CSplinesNorm';
    %model.isi.cpx=[-small 0 1:2:25 30:5:45  50:10:70 80+small 81]/1000;
    model.isi.cpx=[-small 0 1:4:25 30:10:70 80+small 81]/1000;
    model.isi.conv= .1;
    model.isi.convper= .1;
    model.isi.outputMode= 'interval'; % other option is 'timesteps' or not setting
    model.isi.outputInterval= 2;  % [sec]
end

out= adaptFilter(data,model,filter);



