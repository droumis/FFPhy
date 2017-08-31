% Cluster decoding SI on the W-track

datadir_1 = '/opt/data13/kkay/___Superfourcode_data/Decode/';
datadir_2 = '/opt/data13/kkay/___Superfourcode_data/Decode2/';

Calculate_1 = 1;   % Sorted
Calculate_2 = 0;   % Clusterless

    % Parameters
    animals_tocalc = {'Bond'}; %,'Government','Frank','Eight','Ten','Conley','Miles','Egypt','Dave','Corriander'}; 
    dayeps = [4 6];  %[6 2; 6 4; 6 6];  % leave empty if want to do all epochs
    remoteW = 0;   % set to 1 if want to decode in other W as well
    
        % calculate_1: sorted
    cellchoice =  1;  % look inside function to see what is specified
    winsize_1 = 0.125; % ** OBSOLETE ** (in seconds, decoding bin size)
    
        % calculate_2: clusterless
    winsize_2 = 0.005;  
    TETSET = [4]; %[1 2 3 4];  % 1: CA1 only, 2: CA2/3 only, 3: CA3 only, 4: All tets
        
if Calculate_1
    for aa = 1:length(animals_tocalc)
        animalname = animals_tocalc{aa};
            animalinfo = animaldef(animalname);
        kk_fourdecode(animalinfo{2},animalinfo{3},dayeps,animalname,...
            'savedir',datadir_1,...
            'cellchoice',cellchoice,...
            'remoteW',remoteW,...
            'winsize',winsize_1);
    end
end

if Calculate_2
    for aa = 1:length(animals_tocalc)
        animalname = animals_tocalc{aa};
            animalinfo = animaldef(animalname);
        for tetset = TETSET
            kk_fourdecode_2(animalinfo{2},animalinfo{3},dayeps,animalname,...
                'savedir',datadir_2,...
                'remoteW',remoteW,...
                'winsize',winsize_2,...
                'TETSET',tetset);
        end
    end
end
