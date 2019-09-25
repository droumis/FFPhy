function selectProRetro(pro_or_retro)
%function selectProRetro(pro_or_retro)

% pro_or_retro= {'Pro', 'Retro', 'HomeArm'}
% 
%  Select cells for pro-/retrospective coding analysis, i.e. only cells that have
%  a placefield on the home arm.  A placefield in on the home arm if the
%  linpos(edges) are between 10 and 70.
%  Also make sure that the corresponding trajectories are included in the
%  analysis.  When the two paired trajectories have different placefields then
%  take the union of both [@?].
%  
%
%  Use this function only after the placefields have been selected.

% define home arm
minpos= 10;
maxpos= 70;

% what to do if place field outside [minpos maxpos]  {'cut-off', 'reject'}; 
outside_homearm= {'reject'}; 

switch pro_or_retro
case 'Pro'
    traj_sel=[0 2];
case 'Retro'
    traj_sel=[1 3];
case 'HomeArm'
    traj_sel=[0:3];
otherwise error('specify ''Pro'', ''Retro'' or ''HomeArm'' as option');
end

traj_pairs= [2 3 0 1];

global fmaux select 
loadFile(fmaux.select);
ncells=resetCellList;
n= 1;
nPF= [];
for n=1:ncells

    if n>length(select.a) | isempty(select.a{n}) 
        continue;
    end

    % save selection criteria and independent variable, so that they can be
    % duplicated (in case there's more than one place field per traj)
    olda= select.a{n};
    oldx= select.x{n};
    select.a{n}= {};
    select.x{n}= {};

    if select.cellnum(n,2) ~= 2; continue; end % select only familiar epoch
    nPFcell= 0;
    trajs=[]; ntrajs= 0; itrajs=[];

    for ia=1:length(olda) % loop over trajectories
        ntrajs= ntrajs+1;
        trajs(ntrajs)= olda{ia}.traj;
        itrajs(ntrajs)= ia;
    end
    
    nnewtraj= 0;
    done= zeros(1,ntrajs);
    newa=[]; newx=[];
    for it=1:ntrajs
        if done(it) continue; end % do not include place field twice
        done(it)= 1;

        % select only requested traj
        if isempty(find(traj_sel== trajs(it)))
            continue;
        end

        a= olda{itrajs(it)};
        x= oldx{itrajs(it)};

        minx= a.linpos(1);
        maxx= a.linpos(2);

        % always analyze entire homearm
        if(minpos < minx & minx < maxpos) minx= minpos; end 
        if(minpos < maxx & maxx < maxpos) maxx= maxpos; end 

        pit= [];
        % check for trajectory pairs that were both included, where placefields
        % overlap
        for jt=it+1:ntrajs
            if trajs(jt)==traj_pairs(trajs(it)+1) | trajs(jt)== trajs(it)
                lo= olda{jt}.linpos(1);
                hi= olda{jt}.linpos(2);
                if (minx < lo & lo < maxx) | (minx < hi & hi < maxx)
                    pit(end+1)= itrajs(jt);
                end
            end
        end

        % calc union of placefields
        if ~isempty(pit)
            lo=[]; hi=[];
            for ip=1:length(pit)
                done(pit(ip))= 1; % mark place field as included
                lo= olda{pit(ip)}.linpos(1);
                hi= olda{pit(ip)}.linpos(2);
                minx= min([lo minx]);
                maxx= max([hi maxx]);
            end
        end

        % deal with place fields outside HomeArm
        if strcmp(outside_homearm,'reject') 
            % reject if placefields extend outside homearm
            if minx < minpos | maxx > maxpos; continue; end
            % cut off placefield outside home arm [minpos maxpos]
            minx= max([minx minpos]);
            maxx= min([maxx maxpos]);
        end

        % make sure the traj are sorted by number
        first_traj= traj_pairs(trajs(it)+1);
        if trajs(it) < first_traj; first_traj= trajs(it); end

        fprintf(1, '  n= %d, traj= %d : %.2f, %.2f\n', n, first_traj, a.linpos');

        a.linpos=[minx; maxx];

        % save first traj
        a.traj= first_traj;
        nnewtraj= nnewtraj+1;
        newa{nnewtraj}= a;
        newx{nnewtraj}= x;

        % save second traj
        a.traj= traj_pairs(first_traj+1);
        nnewtraj= nnewtraj+1;
        newa{nnewtraj}= a;
        newx{nnewtraj}= x;
    end
    % save new list
    select.a{n}= newa;
    select.x{n}= newx;

%    keyboard
    disp('---')
end
    
[select,totalana]= cleanSelect(select);
% save to disk
save([fmaux.select '-' pro_or_retro], 'select','totalana')

