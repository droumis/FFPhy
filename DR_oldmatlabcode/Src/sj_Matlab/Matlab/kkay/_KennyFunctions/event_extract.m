

%% Barack
directoryname = '/data12/kkay/Bar/';
fileprefix = 'bar';
for days=1:22
            kk_suprathetadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/suprathetafilter.mat']))   
            deltadayprocess(directoryname,fileprefix,days,'f','/usr/local/filtering/deltafilter.mat')

end


%% Calvin
directoryname = '/data12/kkay/Cal/';
fileprefix = 'cal';
for days=1:22
        try
            kk_suprathetadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/suprathetafilter.mat']))   
            deltadayprocess(directoryname,fileprefix,days,'f','/usr/local/filtering/deltafilter.mat')
        catch
        end
end


%% Dwight
directoryname = '/data12/kkay/Dwi/';
fileprefix = 'dwi';
for days=1:22
        try
            kk_suprathetadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/suprathetafilter.mat']))   
            deltadayprocess(directoryname,fileprefix,days,'f','/usr/local/filtering/deltafilter.mat')
        catch
        end
end

%% Arnold - Annabelle
directoryname = '/datatmp/kkay/Arn/';
fileprefix = 'arn';
for days=1:22
        try
            kk_suprathetadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/suprathetafilter.mat']))   
            deltadayprocess(directoryname,fileprefix,days,'f','/usr/local/filtering/deltafilter.mat')
        catch
        end
end



%% Frank -- Mattias
directoryname = '/vortexdata/kkay/Fra';
fileprefix = 'fra';
for days=1:22
        try
            thetagnddayprocess(directoryname,fileprefix,days)
        catch
        end
end


%% Chapati
directoryname = '/data12/kkay/Cha/';
fileprefix = 'cha';
days = 1:9;
tetrodes = 1:14; % specify all tetrodes

%% Corriander
directoryname = '/data12/kkay/Cor/';
fileprefix = 'Cor';
days = 1:9;
tetrodes = 1:25; % specify all tetrodes

%% Bashir
directoryname = '/data12/kkay/Bas/';
fileprefix = 'bas';
days = 4:5;
tetrodes = 1:14; % specify all tetrodes

%% Bon
directoryname = '/data12/kkay/Bon/';
fileprefix = 'bon';
days = 6;
tetrodes = 1:30; % specify all tetrodes

%% Fra
directoryname = '/vortexdata/kkay/Fra';
fileprefix = 'fra';
days = 1:12;
tetrodes = 1:30; % specify all tetrodes


%% Ann & Mitt

directoryname = '/data13/anna/Ann/';
fileprefix = 'ann';
days = 14;
tetrodes = 1:31; % specify all tetrodes

directoryname = '/data13/anna/Mitt/';
fileprefix = 'mit';
days = 13;
tetrodes = 1:31; % specify all tetrodes


%% Egypt
directoryname = '/data12/mari/Egy/';
fileprefix = 'egy';
days = 1:12;
tetrodes = 1:21; % specify all tetrodes


%% Florida
directoryname = '/data12/kkay/Flo';
fileprefix = 'flo';
days = 1;
tetrodes = 1:28; % specify all tetrodes




%%% for probe data
kk_filterdecimate(directoryname, fileprefix, days);       %% for 30 kHz data!

%% dayprocesses
     
thetadayprocess(directoryname,fileprefix,days)    %% (!!) use the one in /usr/local/filtering -- NOT the one in my own personal filtering folder!

kk_suprathetadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/Src/Matlab/Filters/suprathetafilter.mat']))   
deltadayprocess(directoryname,fileprefix,days,'f','/usr/local/filtering/deltafilter.mat')
thetagnddayprocess(directoryname,fileprefix,days)

sj_rippledayprocess(directoryname,fileprefix,days);   
lowgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/matlab/Filters/bonlowgammafilter.mat']));
lowgammadayprocess(directoryname,fileprefix,days,'f',(['/data12/kkay/Ann-Mitt-both-.1-1/lowgamma2040filter.mat']));
kk_highgammadayprocess(directoryname,fileprefix,days,'f','/home/kkay/matlab/Filters/bashighgammafilter.mat');
kk_fastfastgammadayprocess(directoryname,fileprefix,days,'f',(['/home/kkay/matlab/Filters/basfastfastgammafilter.mat']));
slowrippledayprocess(directoryname,fileprefix,5);
    % (!!) quickly modified slowrippledayprocess to change 'area' to 'area_tet', and also commented out area conditional

%%

%% ripples
mindur = 0.015;                                
nstd = 3;                                             

for d = days
    extractripples(directoryname, fileprefix, d, tetrodes, mindur, nstd);       % modify extractripples directly to do daymean
end

%% ripples -- disregarding position
mindur = 0.015;                                
nstd = 7;                                             

for d = days
    kk_extractripples_nopos(directoryname, fileprefix, d, tetrodes, mindur, nstd);       % modify extractripples directly to do daymean
end

%% gamma                                         
for d = days
    gammadayprocess(directoryname,fileprefix,days,'assignphase',1)
end

%% lowgamma
nstd = 2; 
mindur = 0.030;                                
for d = days
    extractlowgamma_derive(directoryname, fileprefix, d, tetrodes, mindur,nstd,'samethreshperday',1)
end

%% highgamma
nstd = 2;
mindur = 0.030;  
for d = days
    extracthighgamma_derive(directoryname, fileprefix, d, tetrodes, mindur, nstd,'samethreshperday',1)
end

%% fastfastgamma
mindur=0.015;
nstd = 2; 
for d = days
    extractfastfastgamma_derive(directoryname, fileprefix, d, tetrodes, mindur, nstd,'samethreshperday',1)
end

%% slowripple
mindur = 0.015;                                
nstd = 2; 
for d = 5
    extractslowripples_derive(directoryname, fileprefix, d, tetrodes, mindur,nstd,'samethreshperday',1)
end


%% hightheta
nstd = 1; % 1 std
mindur=1;
maxpeak=1000;
for d = days
    extracthightheta(directoryname, fileprefix, d, tetrodes, mindur, nstd, maxpeak)
end




