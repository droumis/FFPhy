
% Sleep box, 13"~33cm, Novel W track, 30"~76cm. Lin track = 113 cm
% Familiar W-track is ~75cm

%%% Probe
% 
% cd '/data25/sjadhav/Probe/PRb/'
% prefix='PRb';
% dir='/data25/sjadhav/Probe/PRb/';
% dir2='/data25/sjadhav/Probe/PRb_direct/';
% 
% unit=0.6;
% sj_dayprocess_cont2_cmpix('01_082312', dir2, prefix, 1, 'cmperpix',[unit, unit]);

%%%%%%%%%%%%% Day Process, Ripple Day Process and makeparms - HPc %%%%%%%%%%%%%%%%%%%%%
% 
% cd '/mnt/data25/sjadhav/HPExpt/HPc/'
% prefix='HPc';
% dir='/mnt/data25/sjadhav/HPExpt/HPc/';
% dir2='/mnt/data25/sjadhav/HPExpt/HPc_direct/';
% unit = 0.64;  %Linear track and familiar W-track: 0.64 cm/pixel all days
% unit_sl = 0.17; % Sleep Box is 0.17 cm/pixel all days
% unit2 = 0.50; % Novel W-track is 0.50 cm/pixel all days
% 
% unit_day1 = [unit_sl,unit,unit_sl,unit,unit_sl,unit,unit_sl];
% unit_day2to5 = [unit_sl,unit,unit_sl,unit,unit_sl];
% unit_day6to8 = [unit_sl,unit,unit_sl,unit2,unit_sl];
% 
% %%% 
% %% SET UP TO AUTOFIND DIRECTORY NAME
%    sj_dayprocess_cont2_cmpix('01_120812', dir2, prefix, 1, 'cmperpix',unit_day1, 'epochs_run',[2 4 6]); 
% %    sj_dayprocess_cont2_cmpix('02_120912', dir2, prefix, 2, 'cmperpix',unit_day2to5);
% %    sj_dayprocess_cont2_cmpix('03_121012', dir2, prefix, 3, 'cmperpix',unit_day2to5);
% %    sj_dayprocess_cont2_cmpix('04_121112', dir2, prefix, 4, 'cmperpix',unit_day2to5);
% %    sj_dayprocess_cont2_cmpix('05_121212', dir2, prefix, 5, 'cmperpix',unit_day2to5);
% %    sj_dayprocess_cont2_cmpix('06_121312', dir2, prefix, 6, 'cmperpix',unit_day6to8);
% %    sj_dayprocess_cont2_cmpix('07_121412', dir2, prefix, 7, 'cmperpix',unit_day6to8);
% %    sj_dayprocess_cont2_cmpix('08_121512', dir2, prefix, 8, 'cmperpix',unit_day6to8);
% 
%    
% %% Task and Linear Day Process done for days 1:8. Added Linear track to Day 1
% %*****************************************************************************************
% %Days 2 to 8, W-track epochs = 2,4. Day1: EPoch2 is Lin Track, and Epochs 4 and 6 are W-track
% days = 1:8;
% for d = 1:length(days)
%     currday = days(d)
%     if currday == 1,
%         epochs_w = [currday 4; currday 6]
%         createtaskstruct(dir2,prefix,epochs_w,'getcoord_wtrack');
%         createtaskstruct(dir2,prefix,[currday 2],'getcoord_lineartrack','overwrite',0);
%     else
%         epochs_w = [currday 2; currday 4] 
%         createtaskstruct(dir2,prefix,epochs_w,'getcoord_wtrack');
%     end   
% end
% for d = 1:length(days)
%      currday = days(d)
%      sj_lineardayprocess(dir2,prefix,currday);
% end
% 
% %% Task updated for all days. linposstate also updated
% %***************************************************
% % Updating Task
% % --------------
% days = 1:8;
% for d = 1:length(days)
%     currday = days(d)
%     if currday == 1,
%         epochs_sl = [1,7]; epochs_run = [2 4 6]; epochs_rest = [3,5]; % type
%         epochs_wtr1 = [4 6]; epochs_lin = [2]; epochs_presl=1; epochs_postsl=7; % env
%         env = 'lin'; sj_updatetaskenv(dir2,prefix,currday,epochs_lin,env);
%     end
%     if currday>=2 && currday <=5
%         epochs_sl = [1,5]; epochs_run = [2 4]; epochs_rest = [3]; % type
%         epochs_wtr1 = [2 4]; epochs_presl=1; epochs_postsl=5; % env
%     end  
%     if currday>=6 
%         epochs_sl = [1,5]; epochs_run = [2 4]; epochs_rest = [3]; % type
%         epochs_wtr1 = [2]; epochs_wtr2 = [4]; epochs_presl=1; epochs_postsl=5; % env
%         env = 'wtr2'; sj_updatetaskenv(dir2,prefix,currday,epochs_wtr2,env);
%     end  
%     % Common for all days
%     % Type fields
%     type = 'sleep'; sj_updatetaskstruct(dir2,prefix,currday,epochs_sl,type);
%     type = 'run'; sj_updatetaskstruct(dir2,prefix,currday,epochs_run,type);
%     type = 'rest'; sj_updatetaskstruct(dir2,prefix,currday,epochs_rest,type);
%     % Environment fields
%     env = 'wtr1'; sj_updatetaskenv(dir2,prefix,currday,epochs_wtr1,env);
%     env = 'presleep'; sj_updatetaskenv(dir2,prefix,currday,epochs_presl,env);
%     env = 'postsleep'; sj_updatetaskenv(dir2,prefix,currday,epochs_postsl,env);   
% end
% 
% Updating linposstate
% % --------------------
% for n=1:8
%     sj_updatelinposstate(dir2,prefix,n);
% end
%   
% 
% %% EEG Should be updated for length and Gnd Added before filtering. 
% %****************************************************************** 
% %% EEG updated and ground added for all days 
% %***************************************************
% days = [1,2,3,4,5,6,7,8];
% epochs = [7,5,5,5,5,5,5,5];  % Tet 22 is split off Tet 21. Not useful
% updatetets = [15:22];
% alltets = [1:22];
% 
% for d = 1:length(days)
%     currday = days(d)
%     currepochs = 1:epochs(currday)
%     sj_eegupdate(prefix,currday, currepochs,1,updatetets);
%     sj_HPexpt_addgndtoeeg(prefix, currday, currepochs, alltets);
% end
% 
% 
% % Ripple band done day 2 and theta band and Thetagnd, done
% % %%%  delta and supratheta done NO
% %***************************************************
% days = [1:8];
% for d = 1:length(days)
%     n = days(d)
%     sj_rippledayprocess(dir2,prefix,n); 
%     sj_thetadayprocess(dir2,prefix,n); 
%     sj_thetadayprocess_gnd(dir2,prefix,n); 
%     %sj_deltadayprocess(dir2,prefix,n); 
%     %sj_suprathetadayprocess(dir2,prefix,n); 
% end
% 
% days = [1:8];
% for d = 1:length(days)
%     n = days(d)
%     sj_deltadayprocess(dir2,prefix,n); 
%     sj_suprathetadayprocess(dir2,prefix,n); 
% end
% 
% %% Do Spindlegnd and Deltagnd for PFC tets, What About Gamma
% %**********************************************************************************
% 
% %% Ripple, Spindle, etc Extract
% %**********************************************************************************
% 
% riptetlist = [1,2,3,4,5,6]; % Only use riptetlist since too much noise on reference.
% 
% for n=1:8
%    %n, sj_extractripples(dir2,prefix,n,-1,0.015,2); % Option -1: all tetrodes are processed
%    
%    % Using byepoch function (temporary), since I have to generate position files for "middle"
%    % rest sessions still.
%    if n==1,
%        epochs = [1,2,4,6,7];
%    else
%        epochs = [1,2,4,5];
%    end
%    n, sjHPexpt_extractripples_byepoch(dir2,prefix,n, riptetlist,0.015,2, epochs)
% end
%  
%  
% 
% %% Spikes
% %*********
%  for n=3:4;   % day number
%      n
%      disp('Making Spike Parameter File')
%     sj_makedayparms_dio(dir, prefix, n); % removes stim artifacts, does PC+extra parameters
%  end
% 
% 
%   
% keyboard;
% 
% 
% 
% %% Always update Tetinfo and Cellinfo after clustering and running dayprocess for spikes
% %************************************************************************************
% createtetinfostruct(dir2,prefix);
% createcellinfostruct(dir2,prefix); 
% 
% 
% %% Need to change to get ripple tets from cells condition. >=2 cells on Tet 
% % Instead, define ripple tets
% % -------------------------------------------------------------------------
% %riptetlist. dCA1  1:6 (7 ref. 2 had cells day 1, then below layer. 5 had cells day1, then MU: use 1,3,4,5,6)
% %            iCA1  None. Missed - in cortex (15 became ref. 16,18,20 good. 17 some cells MU: use 16,17,18,20 ) 
%  riptetlist = [1,2,3,4,5,6]; % 6 tets, 
%  sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet'); 
%  sj_addtetrodelocation(dir2,prefix,1:7,'CA1'); sj_addtetrodelocation(dir2,prefix,8:16,'PFC'); sj_addtetrodelocation(dir2,prefix,17:22,'Ctx');
%  sj_addtetrodedescription(dir2,prefix,[7],'CA1Ref'); % Ref 7 for CA1 throughout. Also, ref for PFC days 1-4
%  sj_addtetrodedescription(dir2,prefix,[17],'CtxRef'); % 17 for Ctx/iCA1 throughout
%  sj_addtetrodedescription_days(dir2,prefix,[13],'PFCRef',[5:8]); % 13 is PFC ref for days 5:8
%  sj_addcellinfotag2(dir2,prefix);  % Add tag2 based on tag or area and meanfiringrate 


%%



%%%%%%%%%%%%% Day Process, Ripple Day Process and makeparms - HPa %%%%%%%%%%%%%%%%%%%%%

% cd '/data25/sjadhav/HPExpt/HPa/'
% prefix='HPa';
% dir='/data25/sjadhav/HPExpt/HPa/';
% dir2='/data25/sjadhav/HPExpt/HPa_direct/';
% unit = 0.65;  %Linear track and familiar W-track: 0.65 cm/pixel all days
% unit_sl = 0.18; % Sleep Box is 0.18 cm/pixel all days
% unit2 = 0.69; % Novel W-track is 0.69 cm/pixel all days
% 
% unit_day1 = [unit_sl,unit,unit_sl,unit,unit_sl,unit,unit_sl];
% unit_day2to5 = [unit_sl,unit,unit_sl,unit,unit_sl];
% unit_day6to8 = [unit_sl,unit,unit_sl,unit2,unit_sl];

%%% 
%% SET UP TO AUTOFIND DIRECTORY NAME
%  sj_dayprocess_cont2_cmpix('01_020112', dir2, prefix, 1, 'cmperpix',unit_day1, 'epochs_run',[2 4 6]); 
%  sj_dayprocess_cont2_cmpix('02_020212', dir2, prefix, 2, 'cmperpix',unit_day2to5);
%  sj_dayprocess_cont2_cmpix('03_020312', dir2, prefix, 3, 'cmperpix',unit_day2to5);
%  sj_dayprocess_cont2_cmpix('04_020412', dir2, prefix, 4, 'cmperpix',unit_day2to5);
%  sj_dayprocess_cont2_cmpix('05_020512', dir2, prefix, 5, 'cmperpix',unit_day2to5);
%  sj_dayprocess_cont2_cmpix('06-020612', dir2, prefix, 6, 'cmperpix',unit_day6to8);
%  sj_dayprocess_cont2_cmpix('07_020712', dir2, prefix, 7, 'cmperpix',unit_day6to8);
%  sj_dayprocess_cont2_cmpix('08_020812', dir2, prefix, 8, 'cmperpix',unit_day6to8);

  
%%% Task and Linear Day Process done for days 1:8. Added Linear track to Day 1
% *****************************************************************************************
% Days 2 to 8, W-track epochs = 2,4. Day1: EPoch2 is Lin Track, and Epochs 4 and 6 are W-track
% days = 1:8;
% for d = 1:length(days)
%     currday = days(d)
%     if currday == 1,
%         epochs_w = [currday 4; currday 6]
%         createtaskstruct(dir2,prefix,epochs_w,'getcoord_wtrack');
%         createtaskstruct(dir2,prefix,[currday 2],'getcoord_lineartrack','overwrite',0);
%     else
%         epochs_w = [currday 2; currday 4] 
%         createtaskstruct(dir2,prefix,epochs_w,'getcoord_wtrack');
%     end   
%     sj_lineardayprocess(dir2,prefix,currday);
% end



%%% Task updated for all days. linposstate also updated
% ***************************************************
% Updating Task
% --------------
% days = 1:8;
% for d = 1:length(days)
%     currday = days(d)
%     if currday == 1,
%         epochs_sl = [1,7]; epochs_run = [2 4 6]; epochs_rest = [3,5]; % type
%         epochs_wtr1 = [4 6]; epochs_lin = [2]; epochs_presl=1; epochs_postsl=7; % env
%         env = 'lin'; sj_updatetaskenv(dir2,prefix,currday,epochs_lin,env);
%     end
%     if currday>=2 && currday <=5
%         epochs_sl = [1,5]; epochs_run = [2 4]; epochs_rest = [3]; % type
%         epochs_wtr1 = [2 4]; epochs_presl=1; epochs_postsl=5; % env
%     end  
%     if currday>=6 
%         epochs_sl = [1,5]; epochs_run = [2 4]; epochs_rest = [3]; % type
%         epochs_wtr1 = [2]; epochs_wtr2 = [4]; epochs_presl=1; epochs_postsl=5; % env
%         env = 'wtr2'; sj_updatetaskenv(dir2,prefix,currday,epochs_wtr2,env);
%     end  
%     % Common for all days
%     % Type fields
%     type = 'sleep'; sj_updatetaskstruct(dir2,prefix,currday,epochs_sl,type);
%     type = 'run'; sj_updatetaskstruct(dir2,prefix,currday,epochs_run,type);
%     type = 'rest'; sj_updatetaskstruct(dir2,prefix,currday,epochs_rest,type);
%     % Environment fields
%     env = 'wtr1'; sj_updatetaskenv(dir2,prefix,currday,epochs_wtr1,env);
%     env = 'presleep'; sj_updatetaskenv(dir2,prefix,currday,epochs_presl,env);
%     env = 'postsleep'; sj_updatetaskenv(dir2,prefix,currday,epochs_postsl,env);   
% end

% Updating linposstate
% --------------------
% for n=1:8
%     sj_updatelinposstate(dir2,prefix,n);
% end

%%% Sleepcmperpix is no longer needed
% -------------------------------------
% %  for n=1:8
% %      sj_updatesleepcmperpix(dir2,prefix,n, unit, unit_sl);
% %  end



%%% EEG Should be updated for length and Gnd Added before filtering. 
% ****************************************************************** 
%%% EEG updated and ground added for all days 1 to 8  
% ***************************************************
% days = [1:8];
% epochs = [7,5,5,5,5,5,5,5];
% updatetets = [15:20];
% alltets = [1:20];
% 
% for d = 1:length(days)
%     currday = days(d)
%     currepochs = 1:epochs(currday)
%     sj_eegupdate(prefix,currday, currepochs,1,updatetets);
%     sj_HPexpt_addgndtoeeg(prefix, currday, currepochs, alltets);
% end



%%% Ripple band and theta band done for all days 1 to 8 
%%% %%% Thetagnd, delta and supratheta done for all days
% ***************************************************
% days = [1:8];
% for d = 1:length(days)
%     n = days(d)
%     %sj_rippledayprocess(dir2,prefix,n); 
%     %sj_thetadayprocess(dir2,prefix,n); 
%     sj_thetadayprocess_gnd(dir2,prefix,n); 
%     sj_deltadayprocess(dir2,prefix,n); 
%     sj_suprathetadayprocess(dir2,prefix,n); 
%  end

%%% Do Spindlegnd and Deltagnd for PFC tets, What About Gamma
% **********************************************************************************

%%% Ripple, Spindle, etc Extract
% **********************************************************************************
%  for n=1:8
%    n, sj_extractripples(dir2,prefix,n,-1,0.015,2); % Option -1: all tetrodes are processed
%  end
%  



% %%% Spikes
% % *********
%  for n=7;   % day number
%      n
%      disp('Making Spike Parameter File')
%     sj_makedayparms_dio(dir, prefix, n); % removes stim artifacts, does PC+extra parameters
%  end


% Cluster qunatificn
% sj_clusterdayprocess(dir, dir2, prefix, 1:8, 2:5);


%%% Always update Tetinfo and Cellinfo after clustering and running dayprocess for spikes
% ************************************************************************
%   createtetinfostruct(dir2,prefix);
%   createcellinfostruct(dir2,prefix);

%%% Need to maybechange to get ripple tets from cells condition. >=2 cells
%%% on Tet. Instead, define ripple tets
% -------------------------------------------------------------------------------
% %riptetlist. dCA1  1,2,4,5,6,7 (3 is ref. 2 never had cells, but was layer. 6 had cells day 1 and MU: use 1,4,5,6,7)
% %             iCA1  8,9,11,12,14 (13 was ctx. all in layer - 12 had cells day 5 onwards: use 8,9,11,12,14) 
%  riptetlist = [1,4,5,6,7,8,9,11,12,14];  % 10 tets, 5+5
%  sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');
%  sj_addtetrodedescription(dir2,prefix,3,'CA1Ref'); % Ref is Tet 3 for everything throughout
%  sj_addtetrodelocation(dir2,prefix,1:7,'CA1'); sj_addtetrodelocation(dir2,prefix,8:14,'iCA1'); sj_addtetrodelocation(dir2,prefix,15:20,'PFC');
%  sj_addcellinfotag2(dir2,prefix);  % Add tag2 based on tag or area and meanfiringrate 



%%% Spectrogram 
% *************
%%% Specgram done for days 1,2,3 - 4,5,6,7,8 
% *****************************************************************



% days = [4:8];
% %days = [1:8];
% epochs = [7,5,5,5,5,5,5,5];
% dotets = 1:20;
% % 
% for d = 1:length(days)
%     
%     currday = days(d)
%     currepochs = 1:epochs(currday)
%     tets = dotets;
%     
%     sj_HPexpt_baselinespecgram(prefix, currday, currepochs, tets,0,'movingwin',[100 10]/1000,'fpass',[0 400]); % normal
%     sj_HPexpt_baselinespecgram_forgnd(prefix, currday, currepochs, tets,'movingwin',[100 10]/1000,'fpass',[0 400]);
%     sj_HPexpt_baselinespecgram(prefix, currday, currepochs, tets,0,'movingwin',[400 40]/1000,'fpass',[0 100]); % mid
%     sj_HPexpt_baselinespecgram_forgnd(prefix, currday, currepochs, tets,'movingwin',[400 40]/1000,'fpass',[0 100]);
%     sj_HPexpt_baselinespecgram(prefix, currday, currepochs, tets,0,'movingwin',[1000 100]/1000,'fpass',[0 40]); % low
%     sj_HPexpt_baselinespecgram_forgnd(prefix, currday, currepochs, tets,'movingwin',[1000 100]/1000,'fpass',[0 40]);
%     sj_HPexpt_baselinespecgram(prefix, currday, currepochs, tets,0,'movingwin',[8000 800]/1000,'fpass',[0 8]); % floor
%     sj_HPexpt_baselinespecgram_forgnd(prefix, currday, currepochs, tets,'movingwin',[8000 800]/1000,'fpass',[0 8]);
%      
% end




%%%%%%%%%%%% Day Process, Ripple Day Process and makeparms -HPb %%%%%%%%%%%%%%%%%%%%%

cd '/data25/sjadhav/HPExpt/HPb/'
prefix='HPb';
dir='/data25/sjadhav/HPExpt/HPb/';
dir2='/data25/sjadhav/HPExpt/HPb_direct/';
unit = 0.56;  %cmperpixel   Linear track and familiar W-track: 0.56 cm/pixel
unit_sl = 0.16; % Sleep Box is 0.16 cm/pixel 
unit2 = 0.70; % Novel W-track is 0.70 cm/pixel

unit_day1 = [unit_sl,unit,unit_sl,unit,unit_sl,unit,unit_sl];
unit_day2to5 = [unit_sl,unit,unit_sl,unit,unit_sl];
unit_day6to8 = [unit_sl,unit,unit_sl,unit2,unit_sl];

%% 

% SET UP TO AUTOFIND DIRECTORY NAME
  sj_dayprocess_cont2_cmpix('01_100412', dir2, prefix, 1, 'cmperpix',unit_day1, 'epochs_run',[2 4 6]); 
  sj_dayprocess_cont2_cmpix('02_100512', dir2, prefix, 2, 'cmperpix',unit_day2to5);
  sj_dayprocess_cont2_cmpix('03_100612', dir2, prefix, 3, 'cmperpix',unit_day2to5);
  sj_dayprocess_cont2_cmpix('04_100712', dir2, prefix, 4, 'cmperpix',unit_day2to5);
   sj_dayprocess_cont2_cmpix('05_100812', dir2, prefix, 5, 'cmperpix',unit_day2to5);
   sj_dayprocess_cont2_cmpix('06_100912', dir2, prefix, 6, 'cmperpix',unit_day6to8);
   sj_dayprocess_cont2_cmpix('07_101012', dir2, prefix, 7, 'cmperpix',unit_day6to8);
   sj_dayprocess_cont2_cmpix('08_101112', dir2, prefix, 8, 'cmperpix',unit_day6to8);



%% Task and Linear Day Process done for days 1:8. Added Linear track to Day 1
*****************************************************************************************
Days 2 to 8, W-track epochs = 2,4. Day1: EPoch2 is Lin Track, and Epochs 4 and 6 are W-track
days = 1:8;
for d = 1:length(days)
    currday = days(d)
    if currday == 1,
        epochs_w = [currday 4; currday 6]
        createtaskstruct(dir2,prefix,epochs_w,'getcoord_wtrack');
        createtaskstruct(dir2,prefix,[currday 2],'getcoord_lineartrack','overwrite',0);
    else
        epochs_w = [currday 2; currday 4] 
        createtaskstruct(dir2,prefix,epochs_w,'getcoord_wtrack');
    end   
    sj_lineardayprocess(dir2,prefix,currday);
end



%% Task updated for all days. linposstate also updated
***************************************************
% Updating Task
--------------
days = 1:8;
for d = 1:length(days)
    currday = days(d)
    if currday == 1,
        epochs_sl = [1,7]; epochs_run = [2 4 6]; epochs_rest = [3,5]; % type
        epochs_wtr1 = [4 6]; epochs_lin = [2]; epochs_presl=1; epochs_postsl=7; % env
        env = 'lin'; sj_updatetaskenv(dir2,prefix,currday,epochs_lin,env);
    end
    if currday>=2 && currday <=5
        epochs_sl = [1,5]; epochs_run = [2 4]; epochs_rest = [3]; % type
        epochs_wtr1 = [2 4]; epochs_presl=1; epochs_postsl=5; % env
    end  
    if currday>=6 
        epochs_sl = [1,5]; epochs_run = [2 4]; epochs_rest = [3]; % type
        epochs_wtr1 = [2]; epochs_wtr2 = [4]; epochs_presl=1; epochs_postsl=5; % env
        env = 'wtr2'; sj_updatetaskenv(dir2,prefix,currday,epochs_wtr2,env);
    end  
    % Common for all days
    % Type fields
    type = 'sleep'; sj_updatetaskstruct(dir2,prefix,currday,epochs_sl,type);
    type = 'run'; sj_updatetaskstruct(dir2,prefix,currday,epochs_run,type);
    type = 'rest'; sj_updatetaskstruct(dir2,prefix,currday,epochs_rest,type);
    % Environment fields
    env = 'wtr1'; sj_updatetaskenv(dir2,prefix,currday,epochs_wtr1,env);
    env = 'presleep'; sj_updatetaskenv(dir2,prefix,currday,epochs_presl,env);
    env = 'postsleep'; sj_updatetaskenv(dir2,prefix,currday,epochs_postsl,env);   
end

% Updating linposstate
--------------------
for n=1:8
    sj_updatelinposstate(dir2,prefix,n);
end

%% Sleepcmperpix is no longer needed
--------------------------------------
% for n=1:1
%     sj_updatesleepcmperpix(dir2,prefix,n, unit, unit_sl);
% end




%%% EEG Should be updated for length and Gnd Added before filtering. 
% ****************************************************************** 
%%% EEG updated and ground added for all days  1 to 8. 
%%% For Day7, Ep2and3 for Tets 17:20 - had to do manual for 1 idx 
% ***************************************************
% days = [1:8];
% epochs = [7,5,5,5,5,5,5,5];
% updatetets = [15:20];
% alltets = [1:20];
% for d = 1:length(days)
%     currday = days(d)
%     currepochs = 1:epochs(currday)
%     sj_eegupdate(prefix,currday, currepochs,1,updatetets);
%     sj_HPexpt_addgndtoeeg(prefix, currday, currepochs, alltets);
% end


%%% Ripple band and theta band done for all days 1:8 
%%% Thetagnd, delta and supratheta done for all days
% ***************************************************
% days = [1:8];
% for d = 1:length(days)
%     n = days(d)
%     %sj_rippledayprocess(dir2,prefix,n); 
%     %sj_thetadayprocess(dir2,prefix,n); 
%     sj_thetadayprocess_gnd(dir2,prefix,n); 
%     sj_deltadayprocess(dir2,prefix,n); 
%     sj_suprathetadayprocess(dir2,prefix,n); 
%     %%% gammadayprocess(dir2,prefix,n);
%     %%% sj_spindledayprocess(dir2,prefix,n,'dognd',1); 
%  end
 

%%% Do Spindlegnd and Deltagnd for PFC tets, What About Gamma
% **********************************************************************************

%%% Extract Ripples, Spindles, Etc
% **********************************************************************************
%   for n=1:8
%       n, sj_extractripples(dir2,prefix,n,-1,0.015,2); % Option -1: all tetrodes are processed
%   end
%spintetlist = [8,9,10,12,14];
%sj_extractspindles(dir2,prefix,1,spintetlist,0.1,2);


% Spikes
% --------

%  for n=8  % day number
%     n
%     disp('Making Spike Parameter File')
%    sj_makedayparms_dio(dir, prefix, n); % removes stim artifacts if asked for, does PC+extra parameters
%  end


% Cluster qunatificn
% sj_clusterdayprocess(dir, dir2, prefix, 1:8, 2:5);


%%% Always update Tetinfo and Cellinfo after clustering and running dayprocess for spikes
% ************************************************************************************
% createtetinfostruct(dir2,prefix);
% createcellinfostruct(dir2,prefix);

  

%%% Need to change to get ripple tets from cells condition. >=2 cells on Tet 
% Instead, define ripple tets
% -------------------------------------------------------------------------
% %riptetlist. dCA1  1,2,4,5,6 (7 ref. 2 had cells day 1, then below layer. 5 had cells day1, then MU: use 1,3,4,5,6)
% %            iCA1  8,9,11,12,14 (15 became ref. 16,18,20 good. 17 some cells MU: use 16,17,18,20 ) 
% riptetlist = [1,3,4,5,6,16,17,18,20]; % 9 tets, 5+4 
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet'); 
% sj_addtetrodelocation(dir2,prefix,1:7,'CA1'); sj_addtetrodelocation(dir2,prefix,8:14,'PFC'); sj_addtetrodelocation(dir2,prefix,15:20,'iCA1');
% sj_addtetrodedescription(dir2,prefix,[7],'CA1Ref'); % Ref 7 for CA1 throughout
% sj_addtetrodedescription_days(dir2,prefix,[15],'iCA1Ref',[3:8]); % 15 for iCA1 starts day3
% sj_addtetrodedescription_days(dir2,prefix,[11],'PFCRef',[1:6]); % 11 is PFC ref for days 1:6
% sj_addtetrodedescription_days(dir2,prefix,[13],'PFCRef',[7:8]); % 13 is PFC ref for days 7:8
% % But note that 11 stays PFC ref for Tet 12 on days 7 and 8 also.
% sj_addcellinfotag2(dir2,prefix);  % Add tag2 based on tag or area and meanfiringrate 



%%% Spectrogram
% *************
% 13 Jun onwards
% ---------------
%%% Specgram done for days 1:7. For 8, normal and mid range is done. 
% ******************************************************************************************************


% days = [1,2,3,5,6];
% days = 8;
% epochs = [7,5,5,5,5,5,5,5];
% dotets = 1:20;
% 
% for d = 1:length(days)
%     
%     currday = days(d)
%     currepochs = 1:epochs(currday)
%     tets = dotets;
%     
% %     sj_HPexpt_baselinespecgram(prefix, currday, currepochs, tets,0,'movingwin',[100 10]/1000,'fpass',[0 400]); % normal
% %     sj_HPexpt_baselinespecgram_forgnd(prefix, currday, currepochs, tets,'movingwin',[100 10]/1000,'fpass',[0 400]);
% %     sj_HPexpt_baselinespecgram(prefix, currday, currepochs, tets,0,'movingwin',[400 40]/1000,'fpass',[0 100]); % mid
% %     sj_HPexpt_baselinespecgram_forgnd(prefix, currday, currepochs, tets,'movingwin',[400 40]/1000,'fpass',[0 100]);
%     sj_HPexpt_baselinespecgram(prefix, currday, currepochs, tets,0,'movingwin',[1000 100]/1000,'fpass',[0 40]); % low
%     sj_HPexpt_baselinespecgram_forgnd(prefix, currday, currepochs, tets,'movingwin',[1000 100]/1000,'fpass',[0 40]);
%     sj_HPexpt_baselinespecgram(prefix, currday, currepochs, tets,0,'movingwin',[8000 800]/1000,'fpass',[0 8]); % floor
%     sj_HPexpt_baselinespecgram_forgnd(prefix, currday, currepochs, tets,'movingwin',[8000 800]/1000,'fpass',[0 8]);
%      
% end









% i=1;
% clear



















