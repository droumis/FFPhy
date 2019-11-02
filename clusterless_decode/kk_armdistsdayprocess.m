
calculate = 1;

if calculate
    %%% select data %%
    animals_tocalc = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
end

if calculate
    for aa = 1:length(animals_tocalc)
        animalname = animals_tocalc{aa};
        animalinfo = animaldef(animalname);
        animaldir = animalinfo{2};
        animalprefix = animalinfo{3};
        
        task = loaddata(animalname,'task');
        
        for day = 1:length(task)
            
            if isempty(task{day})
                continue
            end
            
            % determine run epochs
            runeps = [];
            for ee = 1:length(task{day})
                if ~isempty(task{day}{ee})
                    if isfield(task{day}{ee},'type')
                        if ~strcmp(task{day}{ee}.type,'sleep')
                            if isfield(task{day}{ee},'environment')
                                env = task{day}{ee}.environment;
                                if any(strcmp(env,{'WTrackA' 'WTrackB' 'TrackA' 'TrackB'}))
                                    runeps = [runeps ee];
                                end
                            end
                        end
                    end
                end
            end
            
            if isempty(runeps)
                continue
            end
            
            linpos = loaddatastruct(animaldir,animalprefix,'linpos',day);
            
            % initialize outputs
            armdists = {};        % two elements since up to two different W-tracks
            centerarmmax = [nan nan];
            rightarmmax = [nan nan];
            
            for hh = 1:length(runeps)
                
                ep = runeps(hh);
                postimevec = linpos{day}{ep}.statematrix.time;
                
                % Basic epoch data
                numpossamps = length(postimevec);
                starttime_ep = postimevec(1);
                endtime_ep = postimevec(end);
                
                % Identify W-track segment indices
                lindist = linpos{day}{ep}.statematrix.lindist;
                seg1 = linpos{day}{ep}.statematrix.segmentIndex==1;
                seg2 = linpos{day}{ep}.statematrix.segmentIndex==2;
                seg3 = linpos{day}{ep}.statematrix.segmentIndex==3;
                seg4 = linpos{day}{ep}.statematrix.segmentIndex==4;
                seg5 = linpos{day}{ep}.statematrix.segmentIndex==5;
                
                % output %%%%%%%%%%%%%%%%%%%
                armdists{day}{ep}.postimevec = postimevec;
                armdists{day}{ep}.armdists = nan(numpossamps,1) ;   % (horizontal vector)   +200 / +400 / +600 linear distance from center well
                        centerarmmax = max( lindist(seg1) );
                armdists{day}{ep}.centerarmmax =  centerarmmax;   % maximum linear distance of center arm
                        rightarmmax = max( lindist(seg4 | seg5) - centerarmmax );
                armdists{day}{ep}.rightarmmax  =  rightarmmax;
                
                armdists{day}{ep}.armdists(seg1)        = lindist(seg1)                 ;  % >200: center arm
                armdists{day}{ep}.armdists(seg4 | seg5) = lindist(seg4 | seg5)          ;  % >400: right arm
                armdists{day}{ep}.armdists(seg2 | seg3) = lindist(seg2 | seg3)    + rightarmmax;  % >600: left arm

            end
            
            
            % save
            cd(animaldir)
            save(sprintf('%sarmdists%02d.mat', animalprefix, day), 'armdists');
            
        end

    end
end







