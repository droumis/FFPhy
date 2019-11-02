% Cluster decoding SI on the W-track

Calculate = 1;   % Clusterless: outpos, finalized (uniform prior)

    % Parameters
    if 1
        animals_tocalc = {'Bond','Frank','Government','Dave','Corriander','Higgs','Miles','Ten','Conley','Eight','Egypt'};  %,'Chapati'
        dayeps = [] ;
        manual_time = []; 
        EXTYPE_torun = [ 2 1 3 4 ]; % 1: LR Pro, 2: LR Ret, 3: L inbound, 4: R inbound
        SKIP_SAVED_FILE = 1;  
    elseif 0
        animals_tocalc = {'Bond'}; % Poster Ex. 3
        dayeps = [6 4]; %[3 2];  %[4 6];  %[6 2; 6 4; 6 6];  % leave empty if want to do all epochs
        manual_time = [350 354.5];  
        EXTYPE_torun = 1;
        SKIP_SAVED_FILE = 0;        
    elseif 1
        animals_tocalc = {'Bond'}; % Poster Ex. 3
        dayeps = [8 4]; %[3 2];  %[4 6];  %[6 2; 6 4; 6 6];  % leave empty if want to do all epochs
        manual_time = [241 244.5]; 244.5 ; 
        EXTYPE_torun = 1;
        SKIP_SAVED_FILE = 0;
    end        

if Calculate
   
    % Encode params
    spikethresh     = 60;  % uV
    exclude_FS      = 0;
    msig            = 20 ;          % Gaus kernel, in uV
    xdel            = 1;            % cm, pos mark space
        xsig        = 8 * xdel;     % Gaus kernel, in cm    
    
    % Decode params
    if 0
        % Random walk (50x speedup)
        nummsbin        = 20;
        nummsoverlap    = 4;  % in ms, decoding bin size / overlap
        Priormodel      = 1;             % 0: uniform
        decodedir = '/opt/data50/kkay/__Decode/LR_decode_20_4_Randwalk_new/';
        if 1
            animals_tocalc = fliplr(animals_tocalc);
        end
    elseif 1
        % Uniform
        nummsbin        = 20;
        nummsoverlap    = 4;  % in ms, decoding bin size / overlap
        Priormodel      = 0;             % 0: uniform        
        decodedir = '/opt/data50/kkay/__Decode/LR_decode_20_4_Uniform/';
        if 0
            animals_tocalc = fliplr(animals_tocalc);
        end        
    end
    Excurbuffer     = 0.5;          % extra time added around excursions
    
    EXCISE = 1;
    
    for EXTYPE = EXTYPE_torun 
        for aa = 1:length(animals_tocalc)
            animalname = animals_tocalc{aa};
            animalinfo = animaldef(animalname);
            kk_clusterlessdecode17(animalinfo{2},animalinfo{3},dayeps,animalname,...
                'savedir',decodedir,'spikethresh',spikethresh,...
                'xdel',xdel,'xsig',xsig,'msig',msig,...
                'nummsbin',nummsbin,'nummsoverlap',nummsoverlap,...
                'SKIP_SAVED_FILE',SKIP_SAVED_FILE,'Priormodel',Priormodel,...
                'manual_time',manual_time,'EXTYPE',EXTYPE,...
                'exclude_FS',exclude_FS,'Excurbuffer',Excurbuffer,...
                'EXCISE',EXCISE);
        end
    end
    
end    
         
    