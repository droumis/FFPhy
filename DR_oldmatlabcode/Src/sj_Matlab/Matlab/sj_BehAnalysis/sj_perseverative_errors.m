
function [all_inbound_logic, all_inbound_wellstend, ntrajs_perday] = sj_perseverative_errors (animdirect, days,epochs,prefix,savedata1)

%% sj_perseverative_errors ('/data25/sjadhav/RippleInterruption/RE1_direct',[1,2],[2 4], 'RE1',1);
%% sj_perseverative_errors ('/data25/sjadhav/RippleInterruption/RE1_direct',[1:5],[2 4], 'RE1',1);
%% Calls sj_day_findinbound for EACH day, and returns vector of
%% inbound_logic (Correct incorrect), wells exited and entered, and number
%% of inbound trajectories per day, AND PERSEVERATIVE ERROR STATS
%% Saves In File For Current Animal
%% Shantanu Jadhav, 03/30/10

if nargin<5
    savedata1=1;
end

directoryname = animdirect;
if (animdirect(end) == '/')
    animdirect = animdirect(1:end-1);
end
cd(animdirect);
clr = {'b','g','c','m','y','k','r'};
%%directoryname = '/data25/sjadhav/RE1_direct';

%%Initialize
all_inbound_logic = [];
all_inbound_wellstend = [];
ntrajs_perday = [];

%% Loop over days
for i=1:length(days)
    
    filename = [prefix 'linpos' num2str(0) num2str(i)];
    load(filename);
    
    n_temp=0; day_nerr=0; day_nturnerr=0; day_nperserr=0; day_corr=0;
    
    for e = 1:length(epochs)
        epoch = epochs(e);
        [inbound_logic, inbound_wellstend] = sj_day_findinbound (linpos,i,epoch);
        
        %% Pers Err %%
        err_idx = find(inbound_wellstend(:,2)~=1);
        turnerr_idx = find(inbound_wellstend(:,1)==inbound_wellstend(:,2));   %% Start well = End well
        
        nerr_epoch(i,e) = length(err_idx);
        nturnerr_epoch(i,e) = length(turnerr_idx); frac_turnerr_epoch(i,e) = nturnerr_epoch(i,e)/length(inbound_logic);
        nperserr_epoch(i,e) = length(err_idx)-length(turnerr_idx); frac_perserr_epoch(i,e) = nperserr_epoch(i,e)/length(inbound_logic);
        
        
        day_nerr = day_nerr+length(err_idx);
        day_nturnerr = day_nturnerr+length(turnerr_idx);
        day_nperserr = day_nperserr+(length(err_idx)-length(turnerr_idx));
        day_corr = day_corr + length(find(inbound_logic==1));
        
        %% Add to cumulative count
        all_inbound_logic = [all_inbound_logic; inbound_logic];
        all_inbound_wellstend = [all_inbound_wellstend; inbound_wellstend];
        n_temp = n_temp + length(inbound_logic);
    end
    ntrajs_perday(i) = n_temp;
    ncorr_perday(i) = day_corr;
    frac_corr(i) = day_corr/n_temp;
    
    nerr(i) = day_nerr;
    nturnerr(i) = day_nturnerr; frac_turnerr(i) = day_nturnerr/n_temp;
    nperserr(i) = day_nperserr; frac_perserr(i) = day_nperserr/n_temp;
    
  
    
end

savedir = '/data25/sjadhav/RippleInterruption/ProcessedData/BehaviorSummIndivAnimals';
if savedata1==1
    %savefile = sprintf('%s/ProcessedData/%s_perserr.mat', animdirect, prefix);
    savefile = sprintf('%s/%s_perserr.mat', savedir, prefix);
    save(savefile);
end

%savefile = [prefix '_inbound'];
%save(savefile,'all_inbound_logic', 'all_inbound_wellstend', 'ntrajs_perday');


