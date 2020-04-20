

function out = paramconfig(user, varargin)
% import params per user
% DKR 2020

user = 'Demetris';
if ~isempty(varargin)
   assign(varargin); 
end
    
switch user
    case 'Demetris'
    out.andef = animaldef(user);
    out.savefigs = 1;
    out.pausefigs = 1;
    out.showfigs = 1;
    out.savefigas = {'png', 'pdf'};
end

end
