function animalinfo = animaldef(animalname)

switch lower(animalname)
    case 'five'
        animalinfo = {'Five', '/data21/mcarr/Fiv/', 'Fiv'};
    case 'six'
        animalinfo = {'Six', '/data21/mcarr/Six/', 'Six'};
    case 'seven'
        animalinfo = {'Seven', '/data21/mcarr/Sev/', 'Sev'};
    case 'eight'
        animalinfo = {'Eight','/data21/mcarr/Eig/', 'Eig'};
    case 'fear'
        animalinfo = {'Fear', '/data21/mcarr/Fea/', 'Fea'};
    case 'ten'
        animalinfo = {'Ten', '/data21/mcarr/Ten/', 'ten'};
    case 'frank'
        animalinfo = {'Frank', '/data21/mcarr/Fra/', 'fra'};
    case 'bond'
        animalinfo = {'Bond', '/data21/mcarr/Bon/', 'bon'};
    case 'conley'
        animalinfo = {'Conley', '/data21/mcarr/Con/','con'};
    case 'dudley'
        animalinfo = {'Dudley', '/data21/mcarr/Dud/', 'dud'};
    case 'miles'
        animalinfo = {'Miles', '/data21/mcarr/Mil/', 'mil'};
    case 'corriander'
        animalinfo = {'Corriander','/data21/mcarr/Cor/','Cor'};
    case 'corrianderno'
        animalinfo = {'Corriander','/data21/monster/Cor/','Cor'};
    case 'cyclops'
        animalinfo = {'Cyclops','/data21/monster/Cyc/','Cyc'};
    case 'dunphy'
        animalinfo = {'Dunphy','/data21/monster/Dun/','Dun'};
    case 'fafnir'
        animalinfo = {'Fafnir','/data21/monster/Faf/','Faf'};
    case 'godzilla'
        animalinfo = {'Godzilla','/data21/monster/God/','God'};
    case 'grendel'
        animalinfo = {'Grendel','/data21/monster/Gre/','Gre'};
    case 'cml'
        animalinfo = {'Cml','/data21/monster/Cml/','Cml'};
    case 'nico'
        animalinfo = {'Nico','/data21/monster/Nic/','Nic'};        
    otherwise
        error(['Animal ',animalname, ' not defined.']);
end
