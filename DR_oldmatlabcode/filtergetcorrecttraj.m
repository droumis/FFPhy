function [out] = filtergetcorrecttraj(animaldir,animalprefix, epochs, includeStates, minV, mintrajl, welldist, maxpasstime, varargin)
% [out] = filtergetcorrecttraj(animaldir,animalprefix, epochs, includeStates, minV)
% selects correct traj; incorrect traj are coded -1 in statevector
%
% assumes animal is performing an alternation task, with a middle well and
% 2, alternating outside wells, though it does not have to be a W-shaped
% track.  assumes task{d}{e}.wellseq exists and gives correct well order
% such that the middle number corresponds to the middle well, the first and third ones correspond to outside wells,
% for instance [2 1 3] would mean the middle well is well 1, the outside
% wells are 2 and 3.  these well numbers are determined in createtaskstruct
%
%   INCLUDESTATES see getlinstate
%   MINV: minimum linear velocity for use in getlinstate
%   MINTRAJL: minimum number of unbroken indices of state for a
%           trajectory to be included, if =0 all traj included, AS
%           recommends min ~30 if looking at well with V==0
%   WELLDIST: distance animal must be from well to calculate pass time, try
%           10 cm
%   MAXPASSTIME: maximum pass duration for included data.  passes that take
%           longer are excluded
%   OPTIONS:
%   correctf, default 1,correctf = 1 to select correct traj, set incorrect == -2,
%           correctf = 0 to select incorrect and set correct ==-2
%           correctf = 2 to set correct = 2 and incorrect = -2
%   trajfollowing, default 0, if 1 selects as following:
%           current traj correct & prev traj correct: C-C = positive traj number
%           current traj correct & prev traj incorrect: IC-C = -2
%           current traj incorrect & prev traj correct: C-IC = -3
%           current traj correct & prev traj incorrect: IC-IC = -4
%   nexttraj, default 0, if 1 selects as following:
%           current traj correct & next traj correct: C-C = positive traj number
%           current traj correct & next traj incorrect: C-IC = -2
%           current traj incorrect & next traj correct: IC-C = -3
%           current traj incorrect & next traj incorrect: IC-IC = -4
%
%    OUT.STATE is a structure with out{d}{e} gives a value for which trajectory the animal was on
% for each time in statematrix.  Odd trajectories are when the animal
% is moving in a positive linear direction, and even trajectories are for
% negative directions, unknown is -1, and incorrect is -2.  this is similar
% to state, the output for getbehavestate, except
% this function assigns incorrect traj -1, whereas getbehavestate does not.

correctf = 1;
trajfollowing = 0;
nexttraj = 0;
correctorder = [];
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'correctf'
            correctf = varargin{option+1};
        case 'nexttraj'
            nexttraj = varargin{option+1};
        case 'trajfollowing'
            trajfollowing = varargin{option+1};
       case 'correctorder'
            correctorder = varargin{option+1};
       otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

%load variables
loaddays = unique(epochs(:,1));
linpos = loaddatastruct(animaldir, animalprefix, 'linpos', loaddays);
task = loaddatastruct(animaldir, animalprefix, 'task', loaddays);

for i = 1:size(epochs,1)
    %get behave state
    [tempstate, lindist] = getbehavestate(linpos, epochs(i,1), epochs(i,2), includeStates, 'minlinvelocity', minV);  %filters linpos.statematrix.traj by includestates
    if correctf == 2
        tempstate = 2* ones(size(linpos{epochs(i,1)}{epochs(i,2)}.statematrix.time));
    end
    
    %calccorrecttraj
    seqMatrix = calccorrecttraj(epochs(i,:), linpos, task, 'correctorder', correctorder);
    
    %calculate PassTimes
    wellchange = seqMatrix(:,1);
    startind = [wellchange(1:end-1)]; %vector of indexes where exit well and traj changes
    endind = [wellchange(2:end-1)-1; wellchange(end)];
    [outPass j]= find(linpos{epochs(i,1)}{epochs(i,2)}.statematrix.linearDistanceToWells < welldist);
    inPass = setdiff([1:length(linpos{epochs(i,1)}{epochs(i,2)}.statematrix.time)]', outPass); %indexes where animal on a pass
    passtimes = [];
    %calculate the time spent for each pass
    exclPass = zeros(length(startind),1);
    for p = 1:length(startind)
        [trajInPass] = inPass(find(inPass>=startind(p) & inPass<=endind(p)));
        if (~isempty(trajInPass))
            newStartInd = trajInPass(1);
            newEndInd = trajInPass(end);
            passtimes = [passtimes; (linpos{epochs(i,1)}{epochs(i,2)}.statematrix.time(newEndInd) - linpos{epochs(i,1)}{epochs(i,2)}.statematrix.time(newStartInd))];
        else
            passtimes = [passtimes; 0];
            exclPass(p) = 1;
        end
    end
    
    [r c] = size(seqMatrix);
    for k=1:(r-1) %length of sequence matrix /num rows
        if trajfollowing == 0 & nexttraj == 0
            if correctf == 1 || correctf == 2
                if seqMatrix(k,3) == 0 %if incorrect traj
                    tempstate( seqMatrix(k,1):(seqMatrix(k+1,1)-1) ) = -2;
                end
            elseif correctf==0
                if seqMatrix(k,3) ~= 0 %if correct traj
                    tempstate( seqMatrix(k,1):(seqMatrix(k+1,1)-1) ) = -2;
                end
            end
        elseif trajfollowing == 1 & k>1
            if seqMatrix(k,3) == 1 & seqMatrix(k-1,3) ==1 %current traj correct and past corret: C-C
                %leave
            elseif seqMatrix(k,3) == 1 & seqMatrix(k-1,3) ==0 %current traj correct and past incorret: IC-C
                tempstate( seqMatrix(k,1):(seqMatrix(k+1,1)-1) ) = -2;
            elseif seqMatrix(k,3) == 0 & seqMatrix(k-1,3) ==1 % C-IC
                tempstate( seqMatrix(k,1):(seqMatrix(k+1,1)-1) ) = -3;
            elseif seqMatrix(k,3) == 0 & seqMatrix(k-1,3) ==0 % IC-IC
                tempstate( seqMatrix(k,1):(seqMatrix(k+1,1)-1) ) = -4;
            end
        elseif nexttraj == 1 
            if seqMatrix(k,3) == 1 & seqMatrix(k+1,3) ==1 %current traj correct and nextcorret: C-C
                %leave
            elseif seqMatrix(k,3) == 1 & seqMatrix(k+1,3) ==0 %current traj correct and next incorret: C-IC
                tempstate( seqMatrix(k,1):(seqMatrix(k+1,1)-1) ) = -2;
            elseif seqMatrix(k,3) == 0 & seqMatrix(k+1,3) ==1 % IC-C
                tempstate( seqMatrix(k,1):(seqMatrix(k+1,1)-1) ) = -3;
            elseif seqMatrix(k,3) == 0 & seqMatrix(k+1,3) ==0 % IC-IC
                tempstate( seqMatrix(k,1):(seqMatrix(k+1,1)-1) ) = -4;
            end
        end
        
        if passtimes(k) > maxpasstime
            tempstate( seqMatrix(k,1):(seqMatrix(k+1,1)-1) ) = -1;
        end
    end
    %last row
    if seqMatrix(r,3) == 0 %if incorrect traj
        tempstate( seqMatrix(r,1):end ) = -2;
    end
    
    statevector = tempstate;
    % statevector(undefind)=-1;
    
    if mintrajl >0 && any(diff(statevector))  %find short trajectories and set state = -1
        l = find(diff(statevector));
        trajind = [ [1; l(1:end-1)] l]; %col 1 is start ind of each traj, col 2 is end ind
        ltraj = trajind(:,2)-trajind(:,1); %length of each trajectory
        short =  find(ltraj<mintrajl); %indexes back into trajind of trajectories that are too short
        for m = 1:length(short)
            statevector(trajind(short(m),1):trajind(short(m),2))=-1;
        end
    end
    
    %output
    out{epochs(i,1)}{epochs(i,2)}.state = statevector;
    out{epochs(i,1)}{epochs(i,2)}.time = linpos{epochs(i,1)}{epochs(i,2)}.statematrix.time;
    out{epochs(i,1)}{epochs(i,2)}.lindist = lindist;
end

end