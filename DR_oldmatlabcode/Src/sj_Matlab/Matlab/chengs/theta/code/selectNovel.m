function selectNovel(selectid)
%function selectNovel(selectid)
% 
% Select only place fields in either the novel arm (novelArm) or familiar
% satellite arm (famArm) during the novel exposure epochs (4,6, ....),
% or all (novel).
% OR the satellite arms in the familiar configuration (fam).
%
%   selectid:    'fam', 'novel', 'novelArm' or 'famArm'

% define satellite arm
%pad= 5;
%armlength= 65;
%minspikes= 200;
minspikes= 1;

% what to do if place field outside [minpos maxpos]  {'cut-off', 'reject'}; 
outside_arm= {'reject'}; 

switch selectid
case 'fam'
case 'novel'
case 'novelArm'
case 'famArm'
otherwise error('specify valid selectid as option');
end

lowcount= 0;
global fmaux select task info
loadFile(fmaux.select);
ncells=resetCellList;
n= 1;
nPF= [];
for n=1:ncells
    d= select.cellnum(n,1); e= select.cellnum(n,2); t= select.cellnum(n,3); c= select.cellnum(n,4);
    if n>length(select.a) | isempty(select.a{n}) 
        continue;
    end

    % save selection criteria and independent variable, so that they can be
    % duplicated (in case there's more than one place field per traj)
    olda= select.a{n};
    oldx= select.x{n};
    select.a{n}= {};
    select.x{n}= {};

    % select epochs
    if strcmp(selectid, 'fam')
        % only familiar epochs
        if e ~= 2 ; continue; end 
    else
        % only novel epochs
        if e < 4 | mod(e,2) ; continue; end 
    
        % determine which arm or trajs are novel
        % 1 is always home arm, 3 is familiar left satellite, 7 is familiar right
        loadVar(fmaux.datadir, 'task', 0, fmaux.prefix, 1);
        if task{d}{e}.task(2)== 3 & task{d}{e}.task(3)== 7
            error('should not get here: no novel arms in this epoch');
        end
    end
    global spikedata
    loadVar(fmaux.data2dir, 'spikedata', d);
    if length(spikedata{d}{e}{t}{c}.time) < minspikes; lowcount= lowcount+1; continue; end
    loadVar(fmaux.data2dir, 'info');


    switch selectid
    case 'novelArm'
        if task{d}{e}.task(2)~= 3
            trajsel= [0 1];
        else
            trajsel= [2 3];
        end
    case 'famArm'
        if task{d}{e}.task(3)~= 7
            trajsel= [0 1];
        else
            trajsel= [2 3];
        end
    case 'novel'
        trajsel= [0:3];
    case 'fam'
        trajsel= [0:3];
    end

%    minpos= info{d}{e}.centerlinpos+pad;
%    maxpos= minpos+armlength;
%    if(maxpos >= info{d}{e}.maxlinpos); error('inconsitency'); end

%@@
    minpos= info{d}{e}.centerlinpos;
    maxpos= info{d}{e}.maxlinpos;

%    disp([minpos, maxpos])
    nPFcell= 0;
    trajs=[]; ntrajs= 0; itrajs=[];

    for ia=1:length(olda) % loop over trajectories
        ntrajs= ntrajs+1;
        trajs(ntrajs)= olda{ia}.traj;
        itrajs(ntrajs)= ia;
    end
    
    nnewtraj= 0;
    newa=[]; newx=[];
    for it=1:ntrajs
        % select only requested traj
        if isempty(find(trajsel== trajs(it)))
            continue;
        end

        a= olda{itrajs(it)};
        x= oldx{itrajs(it)};

        minx= a.linpos(1);
        maxx= a.linpos(2);

        % deal with place fields outside requested area
        if strcmp(outside_arm,'reject') 
            % reject if placefields extend outside arm
            if minx < minpos | maxx > maxpos; continue; end
            % cut off placefield outside home arm [minpos maxpos]
            minx= max([minx minpos]);
            maxx= min([maxx maxpos]);
        end

        a.linpos=[minx; maxx];

        % save first traj
        a.traj= trajs(it);
        nnewtraj= nnewtraj+1;
        newa{nnewtraj}= a;
        newx{nnewtraj}= x;
    end
    % save new list
    select.a{n}= newa;
    select.x{n}= newx;
end
    
[select,totalana]= cleanSelect(select);
% save to disk
fname= [fmaux.select '-' selectid];
save(fname, 'select','totalana')
%showSelect([fmaux.selectid '-' selectid]);
lowcount
