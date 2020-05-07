%%%%% dayprocesses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% conventional %%%%
thetadayprocess(directoryname,fileprefix,days)    %% (!!) use the one in /usr/local/filtering -- NOT the one in my own personal filtering folder!
deltadayprocess(directoryname,fileprefix,days,'f','/usr/local/filtering/deltafilter.mat')
sj_rippledayprocess(directoryname,fileprefix,days);   

%%%% specialized %%%%%
kk_suprathetadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/suprathetafilter.mat']))   
thetagnddayprocess(directoryname,fileprefix,days)
lowgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/matlab/Filters/bonlowgammafilter.mat']));
lowgammadayprocess(directoryname,fileprefix,days,'f',(['/data12/kkay/Ann-Mitt-both-.1-1/lowgamma2040filter.mat']));
kk_highgammadayprocess(directoryname,fileprefix,days,'f','/home/kkay/matlab/Filters/bashighgammafilter.mat');
kk_fastfastgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/matlab/Filters/basfastfastgammafilter.mat']));
slowrippledayprocess(directoryname,fileprefix,5);
    % (!!) quickly modified slowrippledayprocess to change 'area' to 'area_tet', and also commented out area conditional
       
%for probe data
kk_filterdecimate(directoryname, fileprefix, days);       %% for 30 kHz data!


%%%%% event extracts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% ripple %%%%%
mindur = 0.015;                                
nstd = 3;                                             

for d = days
    extractripples(directoryname, fileprefix, d, tetrodes, mindur, nstd);       % modify extractripples directly to do daymean
    %kk_extractripples_nopos(directoryname, fileprefix, d, tetrodes, mindur, nstd);       % modify extractripples directly to do daymean
end

%%%% lowgamma + highgamma + fast-fast + slowripple %%%%%
nstd = 2; 
mindur_gamma = 0.025;
mindur_ripplelike = 0.015

for d = days
    extractlowgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd,'samethreshperday',1)
    extracthighgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma, nstd,'samethreshperday',1)
    %extractfastfastgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_ripplelike, nstd,'samethreshperday',1)
    %extractslowripples_derive(directoryname, fileprefix, d, tetrodes, mindur_ripplelike,nstd,'samethreshperday',1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%% ANIMALS %%%%%%%%%%%%%%%%%%%%%%%%%%%
% below are the last-used calls to dayprocess + event extracts


%% FRANK - Mattias
directoryname = '/data12/kkay/Fra';
fileprefix = 'fra';
tetrodes = 1:30;

% dayprocess
for days=1:12
            thetagnddayprocess(directoryname,fileprefix,days)
            kk_highgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/frahighgammafilter.mat']));
            lowgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/fralowgammafilter.mat']));

end

% events

% ripples were inherited from Maggie's Fra animdirect
 
nstd = 2; 
mindur_gamma = 0.025;                                
for d = 1:12
    extractlowgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd,'samethreshperday',0)
    extracthighgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd,'samethreshperday',0)
end





%% EGYPT - Kenny & Mari
directoryname = '/data12/mari/Egy/';
fileprefix = 'egy';
days = 1:12;
tetrodes = 1:21; % specify all tetrodes

% dayprocess
lowgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Matlab/Filters/bonlowgammafilter.mat']));
kk_highgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/bonhighgammafilter.mat']));

% ripple  
nstd = 2;
mindur_ripple = 0.015;    
mindur_gamma = 0.025; 

for d = days
    extractripples(directoryname, fileprefix, d, tetrodes, mindur_ripple, nstd);       % modify extractripples directly to do daymean
    extractlowgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd,'samethreshperday',0)
    extracthighgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd,'samethreshperday',0)
end


%% HIGGS - Kenny & Mari
directoryname = '/opt/data40/mari/Hig/';
fileprefix = 'hig';
days = [8 11 15];
tetrodes = 1:21;

nstd = 3;
mindur_ripple = 0.015;    
mindur_gamma = 0.025; 

for d = days
    extractripples(directoryname, fileprefix, d, tetrodes, mindur_ripple, nstd);       % modify extractripples directly to do daymean
%     extractlowgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd,'samethreshperday',0)
%     extracthighgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd,'samethreshperday',0)
end


%% CHAPATI - Kenny
directoryname = '/data12/kkay/Cha/';
fileprefix = 'cha';
days = 1:9;
tetrodes = 1:14; % specify all tetrodes

% dayprocesses
lowgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Matlab/Filters/bonlowgammafilter.mat']));
kk_highgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/bonhighgammafilter.mat']));

% events
nstd = 2;
mindur_ripple = 0.015;    
mindur_gamma = 0.025; 
                                           
for d = days
    extractripples(directoryname, fileprefix, d, tetrodes, mindur, nstd);       % modify extractripples directly to do daymean
    extractlowgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd,'samethreshperday',0)
    extracthighgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd,'samethreshperday',0)
end




%% CORRIANDER - Maggie
directoryname = '/data12/kkay/Cor/';
fileprefix = 'Cor';
days = 1:9;
tetrodes = 1:25; % specify all tetrodes

%%%% ripple %%%%%
 % ripples were inherited from Maggie's Fra animdirect, due to odd 
 
% low gamma
lowgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Matlab/Filters/Corlowgammafilter.mat']));
kk_highgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/Corhighgammafilter.mat']));

nstd = 2; 
mindur_gamma = 0.025;                             
for d = days
    extractlowgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd,'samethreshperday',0)
    extracthighgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd,'samethreshperday',0)
end
 

%% BOND - Mattias
directoryname = '/data12/kkay/Bon/';
fileprefix = 'bon';
days = 3:10;
tetrodes = 1:30; % specify all tetrodes

% dayprocess
% low gamma dayprocess eeg output was inherited from Maggie
kk_highgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/bonhighgammafilter.mat']));
    

% events
nstd = 2;
mindur_ripple = 0.015;   
mindur_gamma = 0.025;                                
for d = days
    extractripples(directoryname, fileprefix, d, tetrodes, mindur_ripple, nstd);       % modify extractripples directly to do daymean
    extractlowgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd,'samethreshperday',0)
    extracthighgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd,'samethreshperday',0)
end






%% BARACK - Annabelle
directoryname = '/data12/kkay/Bar/';
fileprefix = 'bar';
for days=1:22
            kk_suprathetadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/suprathetafilter.mat']))   
            deltadayprocess(directoryname,fileprefix,days,'f','/usr/local/filtering/deltafilter.mat')

end


%% CALVIN - Annabelle
directoryname = '/data12/kkay/Cal/';
fileprefix = 'cal';
for days=1:22
        try
            kk_suprathetadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/suprathetafilter.mat']))   
            deltadayprocess(directoryname,fileprefix,days,'f','/usr/local/filtering/deltafilter.mat')
        catch
        end
end


%% DWIGHT - Annabelle
directoryname = '/data12/kkay/Dwi/';
fileprefix = 'dwi';
for days=1:22
        try
            kk_suprathetadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/suprathetafilter.mat']))   
            deltadayprocess(directoryname,fileprefix,days,'f','/usr/local/filtering/deltafilter.mat')
        catch
        end
end

%% ARNOLD - Annabelle
directoryname = '/datatmp/kkay/Arn/';
fileprefix = 'arn';
for days=1:22
        try
            kk_suprathetadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/suprathetafilter.mat']))   
            deltadayprocess(directoryname,fileprefix,days,'f','/usr/local/filtering/deltafilter.mat')
        catch
        end
end



%% BASHIR - Kenny
directoryname = '/data12/kkay/Bas/';
fileprefix = 'bas';
days = 4:5;
tetrodes = 1:14; % specify all tetrodes

%% ANN & MITT - Anna

directoryname = '/data13/anna/Ann/';
fileprefix = 'ann';
days = 14;
tetrodes = 1:31; % specify all tetrodes

directoryname = '/data13/anna/Mitt/';
fileprefix = 'mit';
days = 13;
tetrodes = 1:31; % specify all tetrodes

%% FLORIDA - Kenny & Dan
directoryname = '/data12/kkay/Flo';
fileprefix = 'flo';
days = 1;
tetrodes = 1:28; % specify all tetrodes








