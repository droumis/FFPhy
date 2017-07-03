datadir_superlin = '/opt/data13/kkay/Superlin_data/';
savedir = '/opt/data13/kkay/__WG';

runcollect = 1;
if runcollect
    placestate = 2;  % superlin timefilter # for place state
    ratethresh_field = [0.5  1  2  5];
    
    animals_torun = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
    animals_order = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander'};
end

if runcollect
   
    st = placestate;
    numthresh = length(ratethresh_field);
    
    adtc = [];
    actmatrix = {}; 
    
    % load one animal at a time, so divide up target_adtc into blocks
    for an = 1:length(animals_torun)
        
        animalname = animals_torun{an};
        animalnum = find(strcmp(animalname,animals_order));
        animalinfo = animaldef(animalname);
    
        disp(animalname)

        % load superlin file
        filename = dir(sprintf('%sSuperlin_%s*',datadir_superlin,animalname(1:3)));
        load([datadir_superlin filename(end).name],'superlin')
        
        % load cellinfo file
        spikes = loaddatastruct(animalinfo{2},animalinfo{3},'spikes');
        
        for reg = 1:5   % holdover from superlin, each hippocampal region (CA1, CA2, CA3, DG, Septal)
            
            if isempty(superlin.data{reg}.dtc_output)
                continue
            end
            
            for c = 1:length(superlin.data{reg}.dtc_output)
               
                adtc = [adtc ; animalnum  superlin.data{reg}.dtc_output(c).index];
                detc = superlin.data{reg}.dtc_output(c).indices;
                eps = superlin.data{reg}.dtc_output(c).indices(:,2)';
                    numeps = length(eps);
                
                outmat = nan(numeps,1 + numthresh) ;  
                    outmat(:,1) = eps';  % first column is the epoch #
                    
                % find traj data for this epoch
                for ee = 1:numeps
                    ep = eps(ee);
                    
                    % Verify that this unit was clustered in spikes data 
                    if isempty(spikes{detc(1,1)}{ep}{detc(1,3)}{detc(1,4)})
                        keyboard
                    end
                    
                    % Search for entry and look at each trajectory
                    for zz = 1:length(superlin.data{reg}.detc_output{st})
                        entryindex = superlin.data{reg}.detc_output{st}(zz).index;
                        if all (  detc(ee,:) ==  entryindex  )
                            data = superlin.data{reg}.detc_output{st}(zz).trajdata;
                             for th = 1:length(ratethresh_field) 
                                 thresh = ratethresh_field(th);
                                 active = 0;  
                                 for traj = 1:length(data)
                                     if ~isempty(data{traj})
                                         if any(data{traj}(:,5) > thresh)
                                             active = 1;
                                         end
                                     end
                                 end
                                 outmat(ee,th + 1) = active;
                             end
                             break
                        end
                    end
                end
                
                actmatrix = [actmatrix ; {outmat} ];
                
            end
        end
        
    end
    
    ACT.animals_order = animals_order;
    ACT.ratethresh_field = ratethresh_field;
    ACT.adtc = adtc;
    ACT.actmatrix = actmatrix;
    ACT.descript = 'list of units + corresponding W-track epochs + different place rate thresholds';
    
    cd(savedir)
    save('ACTIVEMATRIX','ACT','-v7.3')
    
end