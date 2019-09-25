function out = JY_getlinvelocity(animaldir,animalprefix,epochs, varargin)
% out = getalllinstate(animaldir,animalprefix,epochs)
% Produces a cell structure with the fields:
% time, velocity (linearized)
%   EPOCHS - N by 2 matrix, columns are [day epoch]
%
%   OPTION: 'smooth', default no smoothing
%                   to compute linear speed, can smooth linear position data.
%                   It is smoothed with a gaussian of length VSW and std VSW/4.
%                   default for lineardayprocess is 2 seconds
%

smooth = [];
if ~isempty(varargin)
    smooth = varargin{2};
end
loaddays = unique(epochs(:,1));
linpos = loaddatastruct(animaldir, animalprefix, 'linpos', loaddays);

% get velocity from data.Pos.correcteddata column 5

pos=loaddatastruct(animaldir, animalprefix, 'data', loaddays);

for i = 1:size(epochs,1)
    out{epochs(i,1)}{epochs(i,2)}= JY_calclinvelocity(pos, epochs(i,:));
    
end

end