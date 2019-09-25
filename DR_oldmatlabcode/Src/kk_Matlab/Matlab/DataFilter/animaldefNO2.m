function animalinfo = animaldef(animalname)

switch lower(animalname)
    
    % Annabelle's animals
    case 'arnold'
        animalinfo = {'Arnold', '/datatmp/kkay/Arn/', 'arn'};
    case 'barack'
        animalinfo = {'Barack', '/data12/kkay/Bar/', 'bar'};
    case 'calvin'
        animalinfo = {'Calvin', '/data12/kkay/Cal/', 'cal'};
    case 'dwight'
        animalinfo = {'Dwight', '/data12/kkay/Dwi/', 'dwi'};
        
    % Kenny's animals
    case 'chapati'
        animalinfo = {'Chapati','/data12/kkay/Cha/','cha'};
    case 'dave'
        animalinfo = {'Dave','/data12/kkay/Dav/','dav'};
    case 'egypt'
        animalinfo = {'Egypt','/data12/mari/Egy/','egy'};
        

    case 'frank'
        animalinfo = {'Frank', '/data12/kkay/Fra/', 'fra'};
    case 'bond'
        animalinfo = {'Bond', '/data12/kkay/Bon/', 'bon'};
        
    % Maggie's animals
    case 'corriander'
        animalinfo = {'Corriander','/data12/kkay/Cor/','Cor'};
        
        
        
        
    
    % Mattias' animals

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

    case 'five'
        animalinfo = {'Five', '/data13/mcarr/Fiv/', 'Fiv'};

    % Ana's animals
    case 'm01'
        animalinfo = {'m01', '/data/ana/M01/', 'm01'};
    case 'm02'
        animalinfo = {'m02', '/data/ana/M02/', 'm02'};
    otherwise
        error(['Animal ',animalname, ' not defined.']);
end
