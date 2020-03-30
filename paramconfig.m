

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
end

end
