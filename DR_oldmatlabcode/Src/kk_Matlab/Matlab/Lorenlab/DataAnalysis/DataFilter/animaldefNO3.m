function [animalinfo, varargout] = animaldef(animalname)
% function animalinfo = animaldef(animalname) [animal info is a cell vector]
% function [animname,fullanimdir,animpref] = animaldef(animalname)
% function [animname,animdir,animpref,finalanimdir] = animaldef(animalname)

load rootdir;

switch lower(animalname)
    case {'dinah','din'}
        animalinfo = {'Dinah', strcat(rootdir,'/Din/'), 'din'};
    case {'phineas','fin'}
        animalinfo = {'Phineas', strcat(rootdir,'/Fin/'), 'fin'};
    case {'ithiel','ith'}
        animalinfo = {'Ithiel', strcat(rootdir,'/Ith/'), 'ith'};
    case {'jabez','jab'}
        animalinfo = {'Jabez', strcat(rootdir,'/Jab/'), 'jab'};
    case {'kohath','ko'}
        animalinfo = {'Kohath', strcat(rootdir,'/Ko/'), 'ko'};
    case {'lemuel','lem'}
        animalinfo = {'Lemuel', strcat(rootdir,'/Lem/'), 'lem'};
    case {'mephibosheth','meph'}
        animalinfo = {'Mephibosheth', strcat(rootdir,'/Meph/'), 'meph'};
    case {'nebuchadnezzar','neb'}
        animalinfo = {'Nebuchadnezzar', strcat(rootdir,'/Neb/'), 'neb'};
    case {'obed','ob'}
        animalinfo = {'Obed', strcat(rootdir,'/Ob/'), 'ob'};
    case {'shantanu','sha'}
        animalinfo = {'Sha', strcat(rootdir,'/Sha/'), 'sha'};
    case {'ten','ten'}
        animalinfo = {'Ten', strcat(rootdir,'/Ten/'), 'ten'};
    otherwise
        error(['Animal ',animalname, ' not defined.']);
end

if nargout >= 3
  varargout{1} = animalinfo{2};
  varargout{2} = animalinfo{3};
  animalinfo = animalinfo{1};
end

if nargout == 4
  andir = varargout{1};
  [aa,bb] = strtok(andir,filesep);
  while ~isempty(bb) & ~strcmp(bb,filesep)
    [aa,bb] = strtok(bb,filesep);
  end
  varargout{3} = aa;
end

