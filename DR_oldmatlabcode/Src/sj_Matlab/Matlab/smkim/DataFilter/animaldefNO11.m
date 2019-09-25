function animalinfo = animaldef(animalname)

switch lower(animalname)
    case 'frank'
        animalinfo = {'Frank', '/data/mkarlsso/Fra/', 'fra'};
    case 'miles'
        animalinfo = {'Miles', '/data/mkarlsso/Mil/', 'mil'};
    case 'nine'
        animalinfo = {'Nine', '/data/mkarlsso/Nin/', 'nin'};
    case 'ten'
        animalinfo = {'Ten', '/data/mkarlsso/Ten/', 'ten'};
    case 'dudley'
        animalinfo = {'Dudley', '/data/mkarlsso/Dud/', 'dud'};
    case 'alex'
        animalinfo = {'Alex', '/data/mkarlsso/Ale/', 'ale'};
    case 'conley'
        animalinfo = {'Conley', '/data/mkarlsso/Con/', 'con'};
    case 'bond'
        animalinfo = {'Bond', '/data/mkarlsso/Bon/', 'bon'};
    case 'five'
        animalinfo = {'Five', '/data13/mcarr/Fiv/', 'Fiv'};
    otherwise
        error(['Animal ',animalname, ' not defined.']);
end
