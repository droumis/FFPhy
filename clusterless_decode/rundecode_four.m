% Cluster decoding SI on the W-track

allan = {'Government','Egypt','Chapati','Dave','Higgs','Frank','Bond','Corriander',...
                'Eight','Ten','Conley','Miles'};
core4   = {'Bond','Frank','Dave','Government'};
anim1st = {'Government','Egypt','Chapati','Dave','Higgs','Frank'};         
anim2nd = {'Bond','Corriander','Eight','Ten','Conley','Miles'};


datadir_20 = '/opt/data13/kkay/___Superfourcode_data/Decode20/';   % kk 3.31.18
datadir_1 = '/opt/data13/kkay/___Superfourcode_data/Decode/';     % original location of Calc_1
datadir_2 = '/opt/data13/kkay/___Superfourcode_data/Decode2/';
datadir_3 = '/opt/data13/kkay/___Superfourcode_data/PCA/';
datadir_4 = '/opt/data50/kkay/__Decode';
datadir_6 = '/opt/data50/kkay/___Decode6';


Calculate_1 = 1;   % Sorted
Calculate_2 = 0;   % Clusterless
Calculate_3 = 0;   % PCA
Calculate_4 = 0;   % Clusterless: trajpos  (trajs 1-4 >> perhaps not the best approach (instead to C, L, R)
Calculate_5 = 0;   % Clusterless: pospos
Calculate_6 = 0;   % Clusterless: dirpos  (Wu Foster convention)
Calculate_7 = 0;   % Clusterless: outpos  
Calculate_8 = 0;   % Clusterless: outpos  (out/in separate) (Dec 2017)
Calculate_10 = 0;   % Clusterless: dironly  (Dec 2017)
Calculate_11 = 0;   % Clusterless: uniform prior, 6-arm  (Dec 2017)
Calculate_13 = 0;   % Clusterless: dironly (Jan 2018), finalized (uniform prior right now)
Calculate_14 = 0;   % Clusterless: sixpos, uniform prior (Jan 2018)
Calculate_12 = 0;   % Clusterless: outpos, under construction (uniform prior right now)

Calculate_17 = 0;   % Clusterless: outpos, finalized (uniform prior)


    EMPIRICAL_TRANSMAT = 0;
    PLOT_TRANSMATRIX = 1;
    DCOFFSET = 1e-10;
    SKIP_SAVED_FILE = 0; 
    
    % Parameters
    if 0
        CP_reduction = 20;
        animals_tocalc = {'Chapati'}; %{'Bond','Government','Frank','Dave'}; %{'Bond'}; %{'Chapati','Egypt','Corriander','Miles','Ten','Conley'}; 
        dayeps = [] ;
        manual_time = []; 
        EXTYPE_torun = [1 3 4];
        SKIP_SAVED_FILE = 1;  
    elseif 0
        animals_tocalc = {'Frank'}; %{'Bond'}; %{'Chapati','Egypt','Corriander','Miles','Ten','Conley'}; 
        dayeps = [4 2] ;%; 4 4; 4 6; 8 4; 8 2; 8 6; 6 4 ; 6 2; 6 6; 7 2; 7 4; 7 6; 9 2; 9 4; 9 6; 5 6; 5 4 ; 5 2; ...
                  %10 2; 10 4; 10 6]; %[3 2];  %[4 6];  %[6 2; 6 4; 6 6];  % leave empty if want to do all epochs
        dec_start = 813; % 134.9939  % this is the ppt start time     
        manualstoptime = 816.5 ; 
    elseif 0
        animals_tocalc = {'Egypt','Corriander','Eight','Ten','Conley','Miles'}; %{'Bond'}; %{'Chapati','Egypt','Corriander','Miles','Ten','Conley'}; 
        dayeps = [ ]; %[3 2];  %[4 6];  %[6 2; 6 4; 6 6];  % leave empty if want to do all epochs
        dec_start = 0;
    elseif 0
        animals_tocalc = {'Dave'}; %{'Bond'}; %{'Chapati','Egypt','Corriander','Miles','Ten','Conley'}; 
        dayeps = [  3 4; 3 2 ; 3 6]; %[3 2];  %[4 6];  %[6 2; 6 4; 6 6];  % leave empty if want to do all epochs
        dec_start = 0;        
    elseif 0
        animals_tocalc = {'Frank'}; % Modest Anti map example 
        dayeps = [6 2]; [4 6]; %[3 2];  %[4 6];  %[6 2; 6 4; 6 6];  % leave empty if want to do all epochs
        dec_start = 740;
    elseif 0
        animals_tocalc = {'Frank'}; % weak example?
        dayeps = [4 2]; %[3 2];  %[4 6];  %[6 2; 6 4; 6 6];  % leave empty if want to do all epochs
        dec_start = 293;     
    elseif 0
        animals_tocalc = {'Bond'}; % one Strong Anti sequence
        dayeps = [4 6]; %[3 2];  %[4 6];  %[6 2; 6 4; 6 6];  % leave empty if want to do all epochs
        dec_start = 145;
    elseif 1
        animals_tocalc = {'Dave'}; % Poster Ex. 1
        dayeps = [3 4]; %[3 2];  %[4 6];  %[6 2; 6 4; 6 6];  % leave empty if want to do all epochs
        manual_time = [ 702 705];
        EXTYPE_torun = 1;
    elseif 1
        animals_tocalc = {'Bond'}; % Poster Ex. 2
        dayeps = [6 4]; %[3 2];  %[4 6];  %[6 2; 6 4; 6 6];  % leave empty if want to do all epochs
        manual_time = [351.5 353]; %
        EXTYPE_torun = 1;
    elseif 1
        animals_tocalc = {'Bond'}; % Poster Ex. 3
        dayeps = [8 4]; %[3 2];  %[4 6];  %[6 2; 6 4; 6 6];  % leave empty if want to do all epochs
        manual_time = [241 244.5]; 244.5 ; 
        EXTYPE_torun = 1;
    end        

        % calculate_1: sorted
    cellchoice =  1;  % look inside function to see what is specified
    
        % calculate_2: clusterless
    winsize_2 = 0.005;  
    TETSET = [5]; %[1 2 3 4];  % 1: CA1 only, 2: CA2/3 only, 3: CA3 only, 4: All tets, 5: Tets w/ clustered principal units
    remoteW = 0;   % set to 1 if want to decode in other W as well


% Calculate_17    % outpos
if Calculate_17
   
    % Encode params
    spikethresh     = 60;  % uV
    exclude_FS      = 1;
    msig            = 20 ;        % Gaus kernel, in uV
    xdel            = 1;  % cm, pos mark space
        xsig        = 8 * xdel;    % Gaus kernel, in cm    
    
    % Decode params
    nummsbin        = 20;
    nummsoverlap    = 4;  % in ms, decoding bin size / overlap
    Priormodel      = 0;             % 0: uniform
    Excurbuffer     = 0.5;          % extra time added around excursions
    
    EXCISE = 1;
    
    for EXTYPE = EXTYPE_torun 
        for aa = 1:length(animals_tocalc)
            animalname = animals_tocalc{aa};
            animalinfo = animaldef(animalname);
            kk_clusterlessdecode17(animalinfo{2},animalinfo{3},dayeps,animalname,...
                'savedir',datadir_4,'spikethresh',spikethresh,...
                'xdel',xdel,'xsig',xsig,'msig',msig,...
                'nummsbin',nummsbin,'nummsoverlap',nummsoverlap,...
                'SKIP_SAVED_FILE',SKIP_SAVED_FILE,'Priormodel',Priormodel,...
                'manual_time',manual_time,'EXTYPE',EXTYPE,...
                'exclude_FS',exclude_FS,'Excurbuffer',Excurbuffer,...
                'EXCISE',EXCISE);
        end
    end
    
end    
         
    
% Calculate_14    % outpos
if Calculate_14
    
    TETSET = 5;
    nummsbin = 20;
    nummsoverlap = 4;
    plot_infunction = 0;
    % encoding parameters
    modelnum = 1;   % 1: moving, outbound, thresh, 2: + non_FS spikes
    spikethresh = 60;  % leave [] if want to use the amplitudes in the spike data
    dt = .001; %postimevec(2) - postimevec(1);
    mdel = 2;  % uV, spacing in amplitude mark space
    xdel = 1;  % cm, spacing in positional mark space
        asig =  10 * mdel ; %6 * mdel;  % gaussian kernel sigma in amplitude mark space
        xsig = 8 * xdel;   % gaussian kernel sigma in positional mark space
        
    % prior parameters
    Priormodel = 0;  % 1: random walk, 2: empirical, smoothed, 3: padded
    DCoffset = 1e-5; %10^-3.5; %1e-2 ; %1e-30;    
    sigma_randomwalk = 4; % sigma of random walk
    
    for aa = 1:length(animals_tocalc)
        animalname = animals_tocalc{aa};
            animalinfo = animaldef(animalname);
        kk_clusterlessdecode14(animalinfo{2},animalinfo{3},dayeps,animalname,...
            'savedir',datadir_4,'TETSET',TETSET,'modelnum',modelnum,'spikethresh',spikethresh,...
            'plot_powertraces',0,'xdel',xdel,'mdel',mdel,...
            'sxker',xsig,'smker',asig,'sigma_randomwalk',sigma_randomwalk,...
            'dt',dt,'plot_infunction',plot_infunction,'nummsbin',nummsbin,'nummsoverlap',nummsoverlap,'dec_start',dec_start,...
            'SKIPSAVED',SKIP_SAVED_FILE,'CELLMAX',CELLMAX,'DCoffset',DCoffset,'Priormodel',Priormodel,...
            'manualstoptime',manualstoptime);
    end
end       
    
         
    
% Calculate_13    % dironly
if Calculate_13
    TETSET = 5;
    nummsbin = 1;
    plot_infunction = 0;
    % encoding parameters
    modelnum = 1;   % 1: moving, outbound, thresh, 2: + non_FS spikes
    spikethresh = 60;  % leave [] if want to use the amplitudes in the spike data
    dt = .001; %postimevec(2) - postimevec(1);
    mdel = 2;  % uV, spacing in amplitude mark space
    xdel = 1;  % cm, spacing in positional mark space
        asig =  6 * mdel ; %6 * mdel;  % gaussian kernel sigma in amplitude mark space
        xsig = 3 * xdel;   % gaussian kernel sigma in positional mark space
    % prior parameters
    Priormodel = 0;  % 1: random walk, 2: empirical, smoothed, 3: padded
    DCoffset = 10^-3.5; %1e-2 ; %1e-30;    
    sigma_randomwalk = 1; % sigma of random walk
    
    for aa = 1:length(animals_tocalc)
        animalname = animals_tocalc{aa};
            animalinfo = animaldef(animalname);
        kk_clusterlessdecode13(animalinfo{2},animalinfo{3},dayeps,animalname,...
            'savedir',datadir_4,'TETSET',TETSET,'modelnum',modelnum,'spikethresh',spikethresh,...
            'plot_powertraces',0,'xdel',xdel,'mdel',mdel,...
            'sxker',xsig,'smker',asig,'sigma_randomwalk',sigma_randomwalk,...
            'dt',dt,'plot_infunction',plot_infunction,'nummsbin',nummsbin,'dec_start',dec_start,...
            'SKIPSAVED',SKIP_SAVED_FILE,'CELLMAX',CELLMAX,'DCoffset',DCoffset,'Priormodel',Priormodel,...
            'manualstoptime',manualstoptime);
    end
end    
    
  
    
% Calculate_7    % outpos
if Calculate_7
    TETSET = 5;
    nummsbin = 1;
    dec_start = 0;
    plot_infunction = 0;
    modelnum = 2;   % 1: all spikes, 2: trajencode, 3: exclude SWR only
    spikethresh = 60;  % leave [] if want to use the amplitudes in the spike data
    dt = .001; %postimevec(2) - postimevec(1);
    sigma_randomwalk = 1.5 ; % sigma for the gaussian used to smooth the transition matrix
    mdel = 2;  % uV, spacing in amplitude mark space
    xdel = 1;  % cm, spacing in positional mark space
        asig =  6 * mdel;  % gaussian kernel sigma in amplitude mark space
        xsig = 3 * xdel;   % gaussian kernel sigma in positional mark space    
    for aa = 1:length(animals_tocalc)
        animalname = animals_tocalc{aa};
            animalinfo = animaldef(animalname);
        kk_clusterlessdecode7(animalinfo{2},animalinfo{3},dayeps,animalname,...
            'savedir',datadir_4,'TETSET',TETSET,'modelnum',modelnum,'spikethresh',spikethresh,...
            'plot_powertraces',0,'xdel',xdel,'mdel',mdel,...
            'sxker',xsig,'smker',asig,'sigma_transmat',sigma_randomwalk,...
            'dt',dt,'plot_infunction',plot_infunction,'nummsbin',nummsbin,'dec_start',dec_start,...
            'SKIPSAVED',SKIP_SAVED_FILE,'CELLMAX',CELLMAX);
    end
end    
        
    
    % Calculate_8    %  dirpos
if Calculate_10
    
    TETSET = 5;
    nummsbin = 10;
    
    plot_infunction = 0;
    modelnum = 2;   % 1: all spikes, 2: trajencode, 3: exclude SWR only
    spikethresh = 60;  % leave [] if want to use the amplitudes in the spike data
    dt = .001; 0.0334; %postimevec(2) - postimevec(1);
    sigma_randomwalk = 1.5 ; % sigma for the gaussian used to smooth the transition matrix
    mdel = 2;  % uV, spacing in amplitude mark space
        asig =  6 * mdel;  % gaussian kernel sigma in amplitude mark space

    Tparms = 0; %[3 4];
    for Tparm = Tparms
        for aa = 1:length(animals_tocalc)
            animalname = animals_tocalc{aa};
            animalinfo = animaldef(animalname);
            kk_clusterlessdecode10(animalinfo{2},animalinfo{3},dayeps,animalname,...
                'savedir',datadir_4,'TETSET',TETSET,'modelnum',modelnum,'spikethresh',spikethresh,...
                'plot_powertraces',0,'mdel',mdel,...
                'smker',asig,'sigma_transmat',sigma_randomwalk,...
                'dt',dt,'plot_infunction',plot_infunction,'nummsbin',nummsbin,'dec_start',dec_start,...
                'SKIPSAVED',SKIP_SAVED_FILE,'CELLMAX',CELLMAX,...
                'EMPIRICAL_TRANSMAT',EMPIRICAL_TRANSMAT,'DCOFFSET',DCOFFSET,...
                'PLOT_TRANSMATRIX',PLOT_TRANSMATRIX,'Tparm',Tparm,'manualstoptime',manualstoptime);
        end
    end

end

% Calculate_8    %  dirpos
if Calculate_11
    
    TETSET = 5;
    nummsbin = 20;
    
    plot_infunction = 0;
    modelnum = 2;   % 1: all spikes, 2: trajencode, 3: exclude SWR only
    spikethresh = 60;  % leave [] if want to use the amplitudes in the spike data
    dt = .001; 0.0334; %postimevec(2) - postimevec(1);
    sigma_randomwalk = 1.5 ; % sigma for the gaussian used to smooth the transition matrix
    mdel = 2;  % uV, spacing in amplitude mark space
    xdel = 1;  % cm, spacing in positional mark space
        asig =  6 * mdel;  % gaussian kernel sigma in amplitude mark space
        xsig = 3 * xdel;   % gaussian kernel sigma in positional mark space    
    for aa = 1:length(animals_tocalc)
        animalname = animals_tocalc{aa};
            animalinfo = animaldef(animalname);
        kk_clusterlessdecode11(animalinfo{2},animalinfo{3},dayeps,animalname,...
            'savedir',datadir_6,'TETSET',TETSET,'modelnum',modelnum,'spikethresh',spikethresh,...
            'plot_powertraces',0,'xdel',xdel,'mdel',mdel,...
            'sxker',xsig,'smker',asig,'sigma_transmat',sigma_randomwalk,...
            'dt',dt,'plot_infunction',plot_infunction,'nummsbin',nummsbin,'dec_start',dec_start,...
            'SKIPSAVED',SKIP_SAVED_FILE,'CELLMAX',CELLMAX,...
            'EMPIRICAL_TRANSMAT',EMPIRICAL_TRANSMAT,'DCOFFSET',DCOFFSET,...
            'PLOT_TRANSMATRIX',PLOT_TRANSMATRIX);
    end
end

    
% Calculate_8    %  dirpos
if Calculate_8
    
    TETSET = 5;
    nummsbin = 20;
    
    plot_infunction = 0;
    modelnum = 2;   % 1: all spikes, 2: trajencode, 3: exclude SWR only
    spikethresh = 60;  % leave [] if want to use the amplitudes in the spike data
    dt = .001; 0.0334; %postimevec(2) - postimevec(1);
    sigma_randomwalk = 1.5 ; % sigma for the gaussian used to smooth the transition matrix
    mdel = 2;  % uV, spacing in amplitude mark space
    xdel = 1;  % cm, spacing in positional mark space
        asig =  6 * mdel;  % gaussian kernel sigma in amplitude mark space
        xsig = 3 * xdel;   % gaussian kernel sigma in positional mark space    
    for aa = 1:length(animals_tocalc)
        animalname = animals_tocalc{aa};
            animalinfo = animaldef(animalname);
        kk_clusterlessdecode9(animalinfo{2},animalinfo{3},dayeps,animalname,...
            'savedir',datadir_6,'TETSET',TETSET,'modelnum',modelnum,'spikethresh',spikethresh,...
            'plot_powertraces',0,'xdel',xdel,'mdel',mdel,...
            'sxker',xsig,'smker',asig,'sigma_transmat',sigma_randomwalk,...
            'dt',dt,'plot_infunction',plot_infunction,'nummsbin',nummsbin,'dec_start',dec_start,...
            'SKIPSAVED',SKIP_SAVED_FILE,'CELLMAX',CELLMAX,...
            'EMPIRICAL_TRANSMAT',EMPIRICAL_TRANSMAT,'DCOFFSET',DCOFFSET,...
            'PLOT_TRANSMATRIX',PLOT_TRANSMATRIX);
    end
end



% % Calculate_4
% if Calculate_4
%     TETSET = 5;
%     nummsbin = 1;
%     dec_start = 650;
%     plot_infunction = 0;
%     modelnum = 2;   % 1: all spikes, 2: trajencode, 3: exclude SWR only
%     spikethresh = 60;  % leave [] if want to use the amplitudes in the spike data
%     dt = .001; 0.0334; %postimevec(2) - postimevec(1);
%     sigma_transmat = 1.5 ; % sigma for the gaussian used to smooth the transition matrix
%     mdel = 2;  % uV, spacing in amplitude mark space
%     xdel = 1;  % cm, spacing in positional mark space
%         smker =  6 * mdel;  % gaussian kernel sigma in amplitude mark space
%         sxker = 3 * xdel;   % gaussian kernel sigma in positional mark space
%     for aa = 1:length(animals_tocalc)
%         animalname = animals_tocalc{aa};
%             animalinfo = animaldef(animalname);
%         kk_clusterlessdecode4(animalinfo{2},animalinfo{3},dayeps,animalname,...
%             'savedir',datadir_4,'TETSET',TETSET,'modelnum',modelnum,'spikethresh',spikethresh,...
%             'plot_powertraces',0,'xdel',xdel,'mdel',mdel,...
%             'sxker',sxker,'smker',smker,'sigma_transmat',sigma_transmat,...
%             'dt',dt,'plot_infunction',plot_infunction,'nummsbin',nummsbin,'dec_start',dec_start,...
%             'SKIPSAVED',SKIPSAVED,'CELLMAX',CELLMAX);
%     end
% end    


% % Calculate_5
% if Calculate_5
%     animals_tocalc = {'Bond'};
%     TETSET = 5;
%     dayeps = [4 6; 4 4; 4 2; 5 6; 9 4; 9 2; 5 4; 5 2]; 
%     nummsbin = 1;
%     dec_start = 0;
%     plot_infunction = 0;
%     modelnum = 2;   % 1: all spikes, 2: trajencode, 3: exclude SWR only
%     spikethresh = 60;  % leave [] if want to use the amplitudes in the spike data
%     dt = .001; 0.0334; %postimevec(2) - postimevec(1);
%     sigma_transmat = 5 ; % sigma for the gaussian used to smooth the transition matrix
%     mdel = 2;  % uV, spacing in amplitude mark space
%     xdel = 1;  % cm, spacing in positional mark space
%         smker =  6 * mdel;  % gaussian kernel sigma in amplitude mark space
%         sxker = 3 * xdel;   % gaussian kernel sigma in positional mark space    
%     for aa = 1:length(animals_tocalc)
%         animalname = animals_tocalc{aa};
%             animalinfo = animaldef(animalname);
%         kk_clusterlessdecode5(animalinfo{2},animalinfo{3},dayeps,animalname,...
%             'savedir',datadir_4,'TETSET',TETSET,'modelnum',modelnum,'spikethresh',spikethresh,...
%             'plot_powertraces',0,'xdel',xdel,'mdel',mdel,...
%             'sxker',sxker,'smker',smker,'sigma_transmat',sigma_transmat,...
%             'dt',dt,'plot_infunction',plot_infunction,'nummsbin',nummsbin,'dec_start',dec_start);
%     end
% end






        
        
if Calculate_1
    for aa = 1:length(animals_tocalc)
        animalname = animals_tocalc{aa};
            animalinfo = animaldef(animalname);
        kk_fourdecode(animalinfo{2},animalinfo{3},dayeps,animalname,...
            'savedir',datadir_20,...
            'cellchoice',cellchoice,'CP_reduction',CP_reduction);
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

if Calculate_3
    for aa = 1:length(animals_tocalc)
        animalname = animals_tocalc{aa};
            animalinfo = animaldef(animalname);
        for tetset = TETSET
            kk_fourdecode_pca(animalinfo{2},animalinfo{3},dayeps,animalname,...
                'savedir',datadir_2)
        end
    end
end
