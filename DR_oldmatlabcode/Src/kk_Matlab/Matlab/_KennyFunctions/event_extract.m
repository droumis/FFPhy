%%%%% dayprocesses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% conventional %%%%
thetadayprocess(directoryname,fileprefix,days)    %% (!!) use the one in /usr/local/filtering -- NOT the one in my own personal filtering folder!
deltadayprocess(directoryname,fileprefix,days,'f','/usr/local/filtering/deltafilter.mat')
sj_rippledayprocess(directoryname,fileprefix,days);   

%%%% specialized %%%%%
kk_suprathetadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/suprathetafilter.mat']))   % 12-14 Hz
kk_spindledayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/spindlefilter.mat']))      % 12-18 Hz,  designeegfilt(1500,12,18)
thetagnddayprocess(directoryname,fileprefix,days)
lowgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/matlab/Filters/bonlowgammafilter.mat']));
lowgammadayprocess(directoryname,fileprefix,days,'f',(['/data12/kkay/Ann-Mitt-both-.1-1/lowgamma2040filter.mat']));
kk_lowgammadayprocess(directoryname,fileprefix,days,'reftognd',1,'f',(['/home/kkay/Src/Matlab/Filters/bonlowgammafilter.mat']));
kk_highgammadayprocess(directoryname,fileprefix,days,'f','/home/kkay/matlab/Filters/bashighgammafilter.mat');
kk_highgammadayprocess(directoryname,fileprefix,days,'reftognd',1,'f',(['/home/kkay/Src/Matlab/Filters/bonhighgammafilter.mat']));
kk_fastgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/basfastfastgammafilter.mat']));
slowrippledayprocess(directoryname,fileprefix,5);
    % (!!) quickly modified slowrippledayprocess to change 'area' to 'area_tet', and also commented out area conditional
       
%for probe data
kk_filterdecimate(directoryname, fileprefix, days);       %% for 30 kHz data!


%%%%% event extracts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% ripple %%%%%
mindur = 0.015;                                
nstd = 2;                                             

for d = days
    extractripples(directoryname, fileprefix, d, tetrodes, mindur, nstd);       % modify extractripples directly to do daymean
    %kk_extractripples_nopos(directoryname, fileprefix, d, tetrodes, mindur, nstd);       % modify extractripples directly to do daymean
end

%%%% lowgamma + highgamma + fast-fast + slowripple %%%%%
nstd = 2; 
mindur_gamma = 0.060;
mindur_ripplelike = 0.015

for d = days
    %extractlowgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd,'samethreshperday',1)
    %extracthighgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma, nstd,'samethreshperday',1)
    %extractfastfastgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_ripplelike, nstd,'samethreshperday',1)
    %extractslowripples_derive(directoryname, fileprefix, d, tetrodes, mindur_ripplelike,nstd,'samethreshperday',1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%% ANIMALS %%%%%%%%%%%%%%%%%%%%%%%%%%%
% below are the last-used calls to dayprocess + event extracts

%% EGYPT - Kenny & Mari
directoryname = '/data12/mari/Egy/';
fileprefix = 'egy';
days = 1:12;
tetrodes = 1:21; % specify all tetrodes

% dayprocess

kk_lowgammadayprocess(directoryname,fileprefix,days,'reftognd',1,'f',(['/home/kkay/Src/Matlab/Filters/bonlowgammafilter.mat']));
kk_highgammadayprocess(directoryname,fileprefix,days,'reftognd',1,'f',(['/home/kkay/Src/Matlab/Filters/bonhighgammafilter.mat']));
deltadayprocess(directoryname,fileprefix,days,'f','/usr/local/filtering/deltafilter.mat')
thetagnddayprocess(directoryname,fileprefix,days)
kk_spindledayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/spindlefilter.mat']))
kk_suprathetadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/suprathetafilter.mat']))
gammadayprocess(directoryname,fileprefix,days,'assignphase',0,'f','/usr/local/filtering/gammafilter.mat')
lowgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Matlab/Filters/bonlowgammafilter.mat']));
kk_highgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/bonhighgammafilter.mat']));
kk_fastgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/matlab/Filters/basfastfastgammafilter.mat']));

% ripple  
nstd_ripple = 3;
nstd_gamma = 2;
mindur_ripple = 0.015;    
mindur_gamma = 0.060; 

for d = days
    %extractripples(directoryname, fileprefix, d, tetrodes, mindur_ripple, nstd_ripple);       % modify extractripples directly to do daymean
    extractlowgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd_gamma,'samethreshperday',1)
    %extracthighgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd_gamma,'samethreshperday',1)
end





%% FRANK - Mattias
directoryname = '/mnt/vortexdata/kkay/Fra/';
fileprefix = 'fra';
tetrodes = 1:30;

% dayprocess
for days=1:12
            kk_lowgammadayprocess(directoryname,fileprefix,days,'reftognd',1,'f',(['/home/kkay/Src/Matlab/Filters/bonlowgammafilter.mat']));
            kk_highgammadayprocess(directoryname,fileprefix,days,'reftognd',1,'f',(['/home/kkay/Src/Matlab/Filters/bonhighgammafilter.mat']));
            %thetagnddayprocess(directoryname,fileprefix,days)
            %kk_highgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/frahighgammafilter.mat']));
            %lowgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/fralowgammafilter.mat']));

end

% events

% ripples were inherited from Maggie's Fra animdirect
 
nstd = 2; 
mindur = 0.060;                                
for d = days
    extractlowgamma_derive(directoryname, fileprefix, d, tetrodes, mindur,nstd,'samethreshperday',1)
    extracthighgamma_derive(directoryname, fileprefix, d, tetrodes, mindur,nstd,'samethreshperday',1)
end









%% CHAPATI - Kenny
directoryname = '/data12/kkay/Cha/';
fileprefix = 'cha';
days = 1:9;
tetrodes = 1:14; % specify all tetrodes

% dayprocesses
%kk_lowgammadayprocess(directoryname,fileprefix,days,'reftognd',1,'f',(['/home/kkay/Src/Matlab/Filters/bonlowgammafilter.mat']));
kk_highgammadayprocess(directoryname,fileprefix,days,'reftognd',1,'f',(['/home/kkay/Src/Matlab/Filters/bonhighgammafilter.mat']));
kk_spindledayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/spindlefilter.mat']))
gammadayprocess(directoryname,fileprefix,days,'assignphase',0,'f','/usr/local/filtering/gammafilter.mat')
lowgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Matlab/Filters/bonlowgammafilter.mat']));
kk_highgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/bonhighgammafilter.mat']));

% events
nstd = 2;
mindur_ripple = 0.015;    
mindur_gamma = 0.060; 
                                           
for d = days
    extractripples(directoryname, fileprefix, d, tetrodes, mindur, nstd);       % modify extractripples directly to do daymean
    extractlowgamma_derive(directoryname, fileprefix, d, tetrodes, mindur,nstd,'samethreshperday',1)
    extracthighgamma_derive(directoryname, fileprefix, d, tetrodes, mindur,nstd,'samethreshperday',1)
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
mindur = 0.060;                                
for d = days
    extractlowgamma_derive(directoryname, fileprefix, d, tetrodes, mindur,nstd,'samethreshperday',1)
    extracthighgamma_derive(directoryname, fileprefix, d, tetrodes, mindur,nstd,'samethreshperday',1)
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
mindur_gamma = 0.060;                                
for d = days
    extractripples(directoryname, fileprefix, d, tetrodes, mindur_ripple, nstd);       % modify extractripples directly to do daymean
    extractlowgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd,'samethreshperday',1)
    extracthighgamma_derive(directoryname, fileprefix, d, tetrodes, mindur_gamma,nstd,'samethreshperday',1)
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








