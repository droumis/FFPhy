function animalinfo = animaldef(animalname)

switch lower(animalname)
    
    % Shantanu's animals
    
    %%i changed this to reading from gideon's computer bc shantanu's may not be connected forever and it's good to centralize.. note however i'm often saving onto my data 19
    
%     case 'hpa'
%         animalinfo = {'HPa', '/mnt/data25/sjadhav/HPExpt/HPa_direct/', 'HPa'};
%     case 'hpb'
%         animalinfo = {'HPb', '/mnt/data25/sjadhav/HPExpt/HPb_direct/', 'HPb'};
%     case 'hpc'
%         animalinfo = {'HPc', '/mnt/data25/sjadhav/HPExpt/HPc_direct/', 'HPc'};
%     case 'hpd'
%         animalinfo = {'HPd', '/mnt/data25/sjadhav/HPExpt/HPd_direct/', 'HPd'};
        
    case 'hpa'
        animalinfo = {'HPa', '/mnt/data15/gideon/HPa_direct/', 'HPa'};
    case 'hpb'
        animalinfo = {'HPb', '/mnt/data15/gideon/HPb_direct/', 'HPb'};
    case 'hpc'
        animalinfo = {'HPc', '/mnt/data15/gideon/HPc_direct/', 'HPc'};
    case 'hpd'
        animalinfo = {'HPd', '/mnt/data15/gideon/HPd_direct/', 'HPd'};
        
        %Gideon's
%   case 'nadal'
%       animalinfo = {'Nadal', '/mnt/data15/gideon/Ndl/', 'Ndl'};
    case 'ndl'
        animalinfo = {'Ndl', '/mnt/data15/gideon/Ndl/', 'Ndl'};
%       animalinfo = {'Nadal', '/mnt/data25/sjadhav/HPExpt/Ndl_direct/', 'Ndl'};
%           animalinfo = {'Nadal', '/data19/droumis/Ndl_DR/', 'Ndl'};
%     case 'rosenthal'
%         animalinfo = {'Rosenthal', '/mnt/data15/gideon/Rtl/', 'Rtl'};
    case 'rtl'
        animalinfo = {'Rtl', '/mnt/data15/gideon/Rtl/', 'Rtl'};
%     case 'borg'
%         animalinfo = {'Borg', '/mnt/data15/gideon/Brg/', 'Brg'};
    case 'brg'
        animalinfo = {'Brg', '/mnt/data15/gideon/Brg/', 'Brg'};
        
    case 'frank'
        animalinfo = {'Frank', '/home/mkarlsso/datanoeeg/Fra/', 'fra'};
    case 'bond'
        animalinfo = {'Bond', '/home/mkarlsso/datanoeeg/Bon/', 'bon'};
    case 'conley'
        animalinfo = {'Conley', '/home/mkarlsso/datanoeeg/Con/', 'con'};

%   case 'ndl'
%       animalinfo = {'Ndl', '/mnt/data25/sjadhav/HPExpt/Ndl_direct/', 'Ndl'};

        
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
        

%     case 'frank'
%         animalinfo = {'Frank', '/data12/kkay/Fra/', 'fra'};
%     case 'bond'
%         animalinfo = {'Bond', '/data12/kkay/Bon/', 'bon'};
%         
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
%     case 'conley'
%         animalinfo = {'Conley', '/data/mkarlsso/Con/', 'con'};

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
