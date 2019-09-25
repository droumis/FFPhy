function out = get2dstate(animaldir,animalprefix,epochs,varargin)
% out = get2dstate(animaldir,animalprefix,epochs,options)
% Produces a cell structure with the calculated fields about 2D location.
% output fields: immobilitytime
% Options: 
%   immobilecutoff (default 2).  The maximum velocity allowed to be considered immobile 
% EPOCHS - N by 2 matrix, columns are [day epoch]
%

immobilecutoff = 2;
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'immobilecutoff'
            immobilecutoff = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

loaddays = unique(epochs(:,1));
pos = loaddatastruct(animaldir, animalprefix, 'pos', loaddays);
for i = 1:size(epochs,1)
    tmpvelocity = pos{epochs(i,1)}{epochs(i,2)}.vel;
    timestep = (pos{epochs(i,1)}{epochs(i,2)}.time(2,1) - pos{epochs(i,1)}{epochs(i,2)}.time(1,1));
    out{epochs(i,1)}{epochs(i,2)}.time = pos{epochs(i,1)}{epochs(i,2)}.time;
    out{epochs(i,1)}{epochs(i,2)}.immobilitytime = cumsum_reset(tmpvelocity<immobilecutoff)*timestep;   
end
%--------------------------------------------------------------------------

function out = cumsum_reset(inputvect)

resetpoints = find(diff(inputvect) < 0)+1;
out = cumsum(inputvect);
for i = 1:length(resetpoints)
    out(resetpoints(i):end) = out(resetpoints(i):end)-out(resetpoints(i)-1);
end