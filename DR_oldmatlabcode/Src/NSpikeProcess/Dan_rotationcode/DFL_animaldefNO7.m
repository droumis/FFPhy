function animalinfo = animaldef(animalname)

switch lower(animalname)
    % Jai's animals {animalname, processed data directory, processed data directory}
    case 'jyb'
        animalinfo = {'JYB', '/data14/jai/JYB_/', 'JYB'};
    case 'h2'
        animalinfo = {'H2', '/home/daliu/H2_/', 'H2'};
    case 'i1'
        animalinfo = {'I1', '/data14/jai/I1_/', 'I1'};

    otherwise
        error(['Animal ',animalname, ' not defined.']);
end
