function animalinfo = animaldef(animalname)

switch animalname
    
    % Ripple Disruption Expt Animals
    
    case 'sjc'
        animalinfo = {'sjc', '/data25/sjadhav/RippleInterruption/sjc_direct/', 'sjc'};
    case 'RE1'
        animalinfo = {'RE1', '/data25/sjadhav/RippleInterruption/RE1_direct/', 'RE1'};
    case 'RNa'
        animalinfo = {'RNa', '/data25/sjadhav/RippleInterruption/RNa_direct/', 'RNa'};
    case 'RNb'
        animalinfo = {'RNb', '/data25/sjadhav/RippleInterruption/RNb_direct/', 'RNb'};
    case 'RNc'
        animalinfo = {'RNc', '/data25/sjadhav/RippleInterruption/RNc_direct/', 'RNc'};
    case 'RNd'
        animalinfo = {'RNd', '/data25/sjadhav/RippleInterruption/RNd_direct/', 'RNd'};
    case 'RCa'
        animalinfo = {'RCa', '/data25/sjadhav/RippleInterruption/RCa_direct/', 'RCa'};
    case 'RCb'
        animalinfo = {'RCb', '/data25/sjadhav/RippleInterruption/RCb_direct/', 'RCb'};
    case 'RCc'
        animalinfo = {'RCc', '/data25/sjadhav/RippleInterruption/RCc_direct/', 'RCc'};
    case 'RCd'
        animalinfo = {'RCd', '/data25/sjadhav/RippleInterruption/RCd_direct/', 'RCd'};
    case 'REc'
        animalinfo = {'REc', '/data25/sjadhav/RippleInterruption/REc_direct/', 'REc'};
    case 'REd'
        animalinfo = {'REd', '/data25/sjadhav/RippleInterruption/REd_direct/', 'REd'};
    case 'REe'
        animalinfo = {'REe', '/data25/sjadhav/RippleInterruption/REe_direct/', 'REe'};
    case 'REf'
        animalinfo = {'REf', '/data25/sjadhav/RippleInterruption/REf_direct/', 'REf'};
    case 'REg'
        animalinfo = {'REg', '/data25/sjadhav/RippleInterruption/REg_direct/', 'REg'};
    case 'REh'
        animalinfo = {'REh', '/data25/sjadhav/RippleInterruption/REh_direct/', 'REh'};
    case 'HPa'
        animalinfo = {'HPa', '/data25/sjadhav/HPExpt/HPa_direct/', 'HPa'};
    case 'HPb'
        animalinfo = {'HPb', '/data25/sjadhav/HPExpt/HPb_direct/', 'HPb'};
    case 'HPc'
        animalinfo = {'HPc', '/data25/sjadhav/HPExpt/HPc_direct/', 'HPc'};
        
    otherwise
        
        error(['Animal ',animalname, ' not defined.']);
end
