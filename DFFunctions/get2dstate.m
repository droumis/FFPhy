function out = get2dstate(animaldir,animalprefix,epochs,varargin)
% out = get2dstate(animaldir,animalprefix,epochs)
% Produces a cell structure with the fields:
% time, headdir, velocity
% EPOCHS - N by 2 matrix, columns are [day epoch]
%

immobility_velocity = 1;  % in cm/sec
immobility_buffer = 0;  % in sec

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'immobility_velocity'
            immobility_velocity = varargin{option+1};
        case 'immobility_buffer'
            immobility_buffer = varargin{option+1};            
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end


loaddays = unique(epochs(:,1));
pos = loaddatastruct(animaldir, animalprefix, 'pos', loaddays);
    
    
for i = 1:size(epochs,1)
    
    timevec = pos{epochs(i,1)}{epochs(i,2)}.data(:,1);
    
    if size(pos{epochs(i,1)}{epochs(i,2)}.data,2) > 5
        if strcmp(pos{epochs(i,1)}{epochs(i,2)}.fields,'time x y dir vel x-sm y-sm dir-sm vel-sm')
            tmpvelocity = pos{epochs(i,1)}{epochs(i,2)}.data(:,9);    % kk 4.19.13
        elseif strcmp(pos{epochs(i,1)}{epochs(i,2)}.fields,'time x y dir vel smooth-x smooth-y smooth-v smooth-dir')
            tmpvelocity = pos{epochs(i,1)}{epochs(i,2)}.data(:,8);  % kk 5.29.13
        elseif strcmp(pos{epochs(i,1)}{epochs(i,2)}.fields,'time x y dir vel smooth-x smooth-y smooth-dir smooth-v')
            tmpvelocity = pos{epochs(i,1)}{epochs(i,2)}.data(:,9);  % kk 5.29.13         
        elseif strcmp(pos{epochs(i,1)}{epochs(i,2)}.fields,'time x-interp y-interp dir-interp vel-gauss x-loess y-loess dir-loess vel-loess')
            tmpvelocity = pos{epochs(i,1)}{epochs(i,2)}.data(:,9);  % dr 4.10.19         
        else
            disp('unknown pos format')
        end
    else
        tmpvelocity = pos{epochs(i,1)}{epochs(i,2)}.data(:,5);
    end
    timestep = (pos{epochs(i,1)}{epochs(i,2)}.data(2,1) - pos{epochs(i,1)}{epochs(i,2)}.data(1,1));
    out{epochs(i,1)}{epochs(i,2)}.time = timevec;
    out{epochs(i,1)}{epochs(i,2)}.headdir = pos{epochs(i,1)}{epochs(i,2)}.data(:,4);
    out{epochs(i,1)}{epochs(i,2)}.velocity = tmpvelocity;
    
    % immobility
    immobilityvec = tmpvelocity < immobility_velocity;
        % apply buffer (to fully eliminate theta during decel and accel periods)
    if immobility_buffer > 0
        immobperiods = vec2list(immobilityvec,timevec);
        immobperiods_buff = [immobperiods(:,1) + immobility_buffer      immobperiods(:,2) - immobility_buffer];
            durations = immobperiods_buff(:,2)-immobperiods_buff(:,1);
            surviveinds = durations > 0;
        immobilityvec = list2vec(immobperiods_buff(surviveinds,:),timevec);
    end
    out{epochs(i,1)}{epochs(i,2)}.immobilitytime = cumsum_reset(immobilityvec)*timestep;
    
    
end
%--------------------------------------------------------------------------

function out = cumsum_reset(inputvect)

resetpoints = find(diff(inputvect) < 0)+1;
out = cumsum(inputvect);
for i = 1:length(resetpoints)
    out(resetpoints(i):end) = out(resetpoints(i):end)-out(resetpoints(i)-1);
end