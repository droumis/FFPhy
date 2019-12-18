

function out = paramconfig(user, varargin)
% import globals per user

user = 'Demetris';
if ~isempty(varargin)
   assign(varargin); 
end
    
switch user
    case 'Demetris'
    out.andef = animaldef(user);
    
end

end
