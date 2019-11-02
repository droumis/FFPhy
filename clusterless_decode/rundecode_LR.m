% Clusterless: outpos, finalized (uniform prior)

Calculate = 1;  

    allan = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander',...
                    'Eight','Ten','Conley','Miles'};
    core4   = {'Bond','Frank','Dave','Government'};
    anim1st = {'Government','Egypt','Chapati','Dave','Higgs','Frank'};         
    anim2nd = {'Bond','Corriander','Eight','Ten','Conley','Miles'};

if 1
    
    if 0
        animals_tocalc = {'Frank'};
    elseif 0
        animals_tocalc = core4;
    elseif 1
        animals_tocalc = allan; 
    elseif 1
        animals_tocalc = fliplr(allan);
    elseif 0
        animals_tocalc = fliplr(anim1st);
    elseif 0
        animals_tocalc = anim2nd;
    end
        
    dayeps          = [] ;
    manual_time     = []; 
    
    EXTYPE_torun    = [7 6]; %[7 6];  % 6: Pro, 7: Ret
    SKIP_SAVED_FILE = 0;  
    
    % THETA_BINS (0 through 4)
    % 0: if Uniform (sliding) decode
    % 1: if Theta phase (pi/2) based 1st-2nd decode
    % 2: if Theta phase (pi/2) based 1st encode (ALL ARMS) , Decode 1st and 2nd
    % 3: if Theta phase (pi/2) based 1st encode (OUTER ARMS ONLY), Decode 1st and 2nd
    % 4: if Theta phase (pi/2) based 1st encode (OUTER ARMS ONLY + CP REDUCED BY 20), Decode 1st and 2nd
    
    if 0      % Uniform (sliding) decode
        
        THETA_BINS      = 0;    % 0: if Uniform (sliding) decode
        savedir = '/opt/data50/kkay/__Decode/ProRet_decode_20_4_Uniform';      % Maxphase-based Thetabins  
        nummsbin        = 20        ;
        nummsoverlap    = 4         ;  % in ms, decoding bin size / overlap       
    
    else      % Selected thetaphase decode
    
        if 1
            THETA_BINS      = 1;    % 1: if Theta phase (pi/2) based 1st-2nd decode  ** re-run on 5.26.18, CA1 MUA theta + symmetric, center-half theta phase, kk
            savedir = '/opt/data50/kkay/__Decode/ProRet_decode_Theta_Uniform_Oppmax_Adj';      % Maxphase-based Thetabins
            disp('go')
        elseif 0
            THETA_BINS      = 2;     % kk runs on 3.29.18
            savedir = '/opt/data50/kkay/__Decode/ProRet_decode_Theta_1st_Encode';      % Maxphase-based Thetabins
        elseif 0
            THETA_BINS     = 3;
            savedir = '/opt/data50/kkay/__Decode/ProRet_decode_Theta_1st_Encode_Outeronly';      % Maxphase-based Thetabins
        elseif 0
            THETA_BINS      = 4;    % kk runs on 3.30.18
            savedir = '/opt/data50/kkay/__Decode/ProRet_decode_Theta_1st_Encode_Outeronly_CP20';      % Maxphase-based Thetabins
        end
        
        % Set this params to NaN if THETA_BINS > 0, as no longer relevant %
        nummsbin        = nan        ;
        nummsoverlap    = nan         ;  % in ms, decoding bin size / overlap
        
    end
    
elseif 1
    
    SKIP_SAVED_FILE = 1;  
    if 1
        animals_tocalc = {'Bond'}; % Poster Ex. 3
        dayeps = [8 4]; 
        manual_time = [241 244.5]; 244.5 ; 
        EXTYPE_torun = 6;
        SKIP_SAVED_FILE = 0;
        THETA_BINS = 1;
        savedir = '/opt/data50/kkay/__Decode/ProRet_decode_Theta_Uniform_Oppmax_Adj';      % Maxphase-based Thetabins
    
    end
    
end

% Calculate_20    % LR
if Calculate

    % Encode params
    spikethresh     = 60        ;  % uV
    exclude_FS      = 0         ;
    msig            = 20        ;  % Gaussian kernel, in uV
    xdel            = 1         ;  % cm, pos mark space
        xsig        = 8 * xdel  ;  % Gaussian kernel, in cm    

    % Decode params
    Priormodel      = 0         ;  % 0: uniform
    Excurbuffer     = 0.5       ;  % extra time added around excursions

    for EXTYPE = EXTYPE_torun 
        for aa = 1:length(animals_tocalc)
            animalname = animals_tocalc{aa};
            animalinfo = animaldef(animalname);
            kk_clusterlessdecode20(animalinfo{2},animalinfo{3},dayeps,animalname,...
                'savedir',savedir,'spikethresh',spikethresh,...
                'msig',msig,...
                'nummsbin',nummsbin,'nummsoverlap',nummsoverlap,...
                'SKIP_SAVED_FILE',SKIP_SAVED_FILE,'Priormodel',Priormodel,...
                'manual_time',manual_time,'EXTYPE',EXTYPE,...
                'exclude_FS',exclude_FS,'Excurbuffer',Excurbuffer,'THETA_BINS',THETA_BINS);
        end
    end

end   







