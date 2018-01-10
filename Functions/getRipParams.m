



function rippms = getRipParams(ripSet)
rippms = setParams(ripSet);
end

%add custom ripple param sets as cases below

function rippms = setParams(ripSet)

switch ripSet
    case 'DR1'
        rippms.consensus_numtets = 1;   % minimum # of tets for consensus event detection
        rippms.minstdthresh = 3;        % STD. how big your ripples are
        rippms.exclusion_dur = 0;  % seconds within which consecutive events are eliminated / ignored
        rippms.minvelocity = 0;
        rippms.maxvelocity = 4;
    case 'DR2'
        rippms.consensus_numtets = 1;   % minimum # of tets for consensus event detection
        rippms.minstdthresh = 5;        % STD. how big your ripples are
        rippms.exclusion_dur = 0.5;  % seconds within which consecutive events are eliminated / ignored
        rippms.minvelocity = 0;
        rippms.maxvelocity = 4;
    case 'DR3'
        rippms.consensus_numtets = 1;   % minimum # of tets for consensus event detection
        rippms.minstdthresh = 3;        % STD. how big your ripples are
        rippms.exclusion_dur = 1;  % seconds within which consecutive events are eliminated / ignored
        rippms.minvelocity = 0;
        rippms.maxvelocity = 4;
    case 'DR4'
        rippms.consensus_numtets = 1;   % minimum # of tets for consensus event detection
        rippms.minstdthresh = 3;        % STD. how big your ripples are
        rippms.exclusion_dur = .5;  % seconds within which consecutive events are eliminated / ignored
        rippms.minvelocity = 0;
        rippms.maxvelocity = 4;
    otherwise
        error('pick an existing rip param set or make a new one')
end
end