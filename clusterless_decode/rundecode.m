
%%% select data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savedir = '/opt/data13/kkay/Superposteriors_data';
animalname = 'Bond';
    animalinfo = animaldef(animalname);
tetfilter = '( isequal($area, ''CA1'') || isequal($area, ''CA3'') )';
dayeps = [4 4];

%%% select what to decode %%%%%%%%%%%%%%%%%%%%%%%%
decodemode = 1; % 1: entire epoch, 2: SWRs
extratime = 500;  % if not doing decodemode 1, the ms around event start to plot


%%% select decoding parameters %%%%%%%%%%%%%%%%%%%%%
old_version = 1;
modelnum = 3;   % 1: all spikes, 2: trajencode, 3: exclude SWR only
spikethresh = 0;
dt = .001; 0.0334; %postimevec(2) - postimevec(1);
sigma_transmat = 5 ; % sigma for the gaussian used to smooth the transition matrix
mdel = 2;  % uV, spacing in amplitude mark space
xdel = 1;  % cm, spacing in positional mark space
    smker =  6 * mdel;  % gaussian kernel sigma in amplitude mark space
    sxker = 3 * xdel;   % gaussian kernel sigma in positional mark space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    



kk_clusterlessdecode(animalinfo{2},animalinfo{3},dayeps,animalname,...
                    'savedir',savedir,'decodemode',decodemode,'tetfilter',tetfilter,'model',model,'spikethresh',spikethresh,'old_version',old_version,...
                    'plot_powertraces',1,'calctraj',1,'xdel',xdel,'mdel',mdel,...
                    'sxker',sxker,'smker',smker,'extratime',extratime,'sigma_transmat',sigma_transmat,...
                    'dt',dt);