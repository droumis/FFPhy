





% ------------------------- RIPPLE DISRUPTION ------------------------------



%%%%%%%%%%%%% Day Process, Ripple Day Process and makeparms %%%%%%%%%%%%%%%%%%%%%

% cd '/data25/sjadhav/RippleInterruption/REg'
% prefix='REg';
% dir='/data25/sjadhav/RippleInterruption/REg/';
% dir2='/data25/sjadhav/RippleInterruption/REg_direct/';
% unit=0.5;  %cmperpixel
% unit_sl=0.17; 
%%% 

% sj_dayprocess_cont2('01_031612', dir2, prefix, 1, 'cmperpix',unit); 
% sj_dayprocess_cont2('02_031712', dir2, prefix, 2, 'cmperpix',unit);
%  sj_dayprocess_cont2('03_031812', dir2, prefix, 3, 'cmperpix',unit);
%  sj_dayprocess_cont2('04_031912', dir2, prefix, 4, 'cmperpix',unit);
%  sj_dayprocess_cont2('05_032012', dir2, prefix, 5, 'cmperpix',unit);
%  sj_dayprocess_cont2('06_032112', dir2, prefix, 6, 'cmperpix',unit);
%  sj_dayprocess_cont2('07_032212', dir2, prefix, 7, 'cmperpix',unit);
%  sj_dayprocess_cont2('08_032312', dir2, prefix, 8, 'cmperpix',unit);

% Linearize position - eg. day 1 epoch 2 and 4
% -----------------------------------------------
% n = 8;   % day number
% createtaskstruct(dir2,prefix,[n 2; n 4],'getcoord_wtrack');
% sj_lineardayprocess(dir2,prefix,n);

% Ripple process and extract - add tetinfo
% -----------------------------------------
% riptetlist = [1,2,5,7,10,11];  riptetlist = [1,2,5,6,7]; 
% sj_rippledayprocess(dir2,prefix,n);
% for n=1:8
%   n, sj_extractripples(dir2,prefix,n,riptetlist,0.015,2);
% end
% createtetinfostruct(dir2,prefix);
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');

% Updating Data
% --------------
% epochs = [1,3,5]; type = 'sleep';
% for n=1:8
%     sj_updatetaskstruct(dir2,prefix,n,epochs,type);
% end
% for n=1:8
%     sj_updatelinposstate(dir2,prefix,n);
% end
%  for n=1:8
%      sj_updatesleepcmperpix(dir2,prefix,n, unit, unit_sl);
%  end
%

% Spikes
% --------

% n=1;   % day number
% sj_makedayparms_dio(dir, prefix, n); % removes stim artifacts, does PC+extra parameters
%
% 

% %After Spike Clustering is done - run dayprocess and then this
% createtetinfostruct(dir2,prefix);
% createcellinfostruct(dir2,prefix);
% riptetlist = [1,2,5,7,10,11]; 
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');

% for n=1:10
%     sj_rippledayprocess(dir2,prefix,n);
% end
%
% for n=1:10
%     thetadayprocess(dir2,prefix,n);
% end

% Cluster qunatificn
% sj_clusterdayprocess(dir, dir2, prefix, 1:8, 2:5);
%


%%%%%%%%%%%%% Day Process, Ripple Day Process and makeparms %%%%%%%%%%%%%%%%%%%%%

% cd '/data25/sjadhav/RippleInterruption/REh'
% prefix='REh';
% dir='/data25/sjadhav/RippleInterruption/REh/';
% dir2='/data25/sjadhav/RippleInterruption/REh_direct/';
% unit=0.5;  %cmperpixel
% unit_sl=0.17; 
%%% 

% sj_dayprocess_cont2('01_032312', dir2, prefix, 1, 'cmperpix',unit); 
% sj_dayprocess_cont2('02_032412', dir2, prefix, 2, 'cmperpix',unit);
% sj_dayprocess_cont2('03_032512', dir2, prefix, 3, 'cmperpix',unit);
% sj_dayprocess_cont2('04_032612', dir2, prefix, 4, 'cmperpix',unit);
% sj_dayprocess_cont2('05_032712', dir2, prefix, 5, 'cmperpix',unit);
% sj_dayprocess_cont2('06_032812', dir2, prefix, 6, 'cmperpix',unit);
% sj_dayprocess_cont2('07_032912', dir2, prefix, 7, 'cmperpix',unit);
% sj_dayprocess_cont2('08_033012', dir2, prefix, 8, 'cmperpix',unit);

% Linearize position - eg. day 1 epoch 2 and 4
% -----------------------------------------------
% n = 6;   % day number
% createtaskstruct(dir2,prefix,[n 2; n 4],'getcoord_wtrack');
% sj_lineardayprocess(dir2,prefix,n);

% Ripple process and extract - add tetinfo
% -----------------------------------------
% sj_rippledayprocess(dir2,prefix,n);
% riptetlist = [1,2,4,5,6]; 
%  for n=6:8
%    n,sj_extractripples(dir2,prefix,n,riptetlist,0.015,2);
%  end
%  createtetinfostruct(dir2,prefix);
%  sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');

% Updating Data
% --------------
%  epochs = [1,3,5]; type = 'sleep';
%  for n=1:8
%      sj_updatetaskstruct(dir2,prefix,n,epochs,type);
%  end
% % for n=1:8
% %     sj_updatelinposstate(dir2,prefix,n);
% % end
%  for n=1:8
%      sj_updatesleepcmperpix(dir2,prefix,n, unit, unit_sl);
%  end
%

% Spikes
% --------

% n=1;   % day number
% sj_makedayparms_dio(dir, prefix, n); % removes stim artifacts, does PC+extra parameters
%
% 

% %After Spike Clustering is done - run dayprocess and then this
% createtetinfostruct(dir2,prefix);
% createcellinfostruct(dir2,prefix);
% riptetlist = [1,2,4,5,6]; 
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');

% for n=1:10
%     sj_rippledayprocess(dir2,prefix,n);
% end
%
% for n=1:10
%     thetadayprocess(dir2,prefix,n);
% end

% Cluster qunatificn
% sj_clusterdayprocess(dir, dir2, prefix, 1:8, 2:5); 











%%%%%%%%%%%%% Day Process, Ripple Day Process and makeparms %%%%%%%%%%%%%%%%%%%%%


% cd '/data25/sjadhav/RippleInterruption/REc'
% prefix='REc';
% dir='/data25/sjadhav/RippleInterruption/REc/';
% dir2='/data25/sjadhav/RippleInterruption/REc_direct/';
% unit=0.46;  %cmperpixel
% unit_sl=0.15; 

% epochs = [1,3,5]; type = 'sleep';
% for n=1:8
%     sj_updatetaskstruct(dir2,prefix,n,epochs,type);
% end
% for n=1:8
%     sj_updatelinposstate(dir2,prefix,n);
% end
% for n=1:8
%     sj_updatesleepcmperpix(dir2,prefix,n, unit, unit_sl);
% end
% 
% riptetlist = [4,6,8,9,10]; 
% for n=1:8
%     %sj_rippledayprocess(dir2,prefix,n);
%     sj_extractripples(dir2,prefix,n,riptetlist,0.015,2);
% end
% createtetinfostruct(dir2,prefix);
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');


 
% sj_dayprocess_cont2('01_112310', dir2, prefix, 1, 'cmperpix',unit); 
% sj_dayprocess_cont2('02_112410', dir2, prefix, 2, 'cmperpix',unit);
% sj_dayprocess_cont2('03_112510', dir2, prefix, 3, 'cmperpix',unit);
% sj_dayprocess_cont2('04_112610', dir2, prefix, 4, 'cmperpix',unit);
% sj_dayprocess_cont2('05_112710', dir2, prefix, 5, 'cmperpix',unit);
% sj_dayprocess_cont2('06_112810', dir2, prefix, 6, 'cmperpix',unit);
% sj_dayprocess_cont2('07_112910', dir2, prefix, 7, 'cmperpix',unit);
% sj_dayprocess_cont2('08_113010', dir2, prefix, 8, 'cmperpix',unit);



% %After Spike Clustering is done
% createtetinfostruct(dir2,prefix);
% createcellinfostruct(dir2,prefix);
% riptetlist = [4,6,8,9,10]; 
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');





% for n=1:8
%     sj_makedayparms_dio(dir, prefix, n);
% end
%
% for n=1:8
%     sj_rippledayprocess(dir2,prefix,n);
% end




% % %------------------------------------------------------------------------
% % 

% cd '/data25/sjadhav/RippleInterruption/REd'
% prefix='REd';
% dir='/data25/sjadhav/RippleInterruption/REd/';
% dir2='/data25/sjadhav/RippleInterruption/REd_direct/';
% unit=0.45;  %cmperpixel
% unit_sl=0.17; 
% 
% epochs = [1,3,5]; type = 'sleep';
% for n=1:11
%     sj_updatetaskstruct(dir2,prefix,n,epochs,type);
% end
% for n=1:11
%     sj_updatelinposstate(dir2,prefix,n);
% end
% for n=1:11
%     sj_updatesleepcmperpix(dir2,prefix,n, unit, unit_sl);
% end
% 
% riptetlist = [3,4,5,6,10,11]; 
% for n=1:11
%     %sj_rippledayprocess(dir2,prefix,n);
%     sj_extractripples(dir2,prefix,n,riptetlist,0.015,2);
% end
% createtetinfostruct(dir2,prefix);
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');



% sj_dayprocess_cont2('01_121710', dir2, prefix, 1, 'cmperpix',unit); 
% sj_dayprocess_cont2('02_121810', dir2, prefix, 2, 'cmperpix',unit);
% sj_dayprocess_cont2('03_121910', dir2, prefix, 3, 'cmperpix',unit);
% sj_dayprocess_cont2('04_122010', dir2, prefix, 4, 'cmperpix',unit);
% sj_dayprocess_cont2('05_122110', dir2, prefix, 5, 'cmperpix',unit);
% sj_dayprocess_cont2('06_122210', dir2, prefix, 6, 'cmperpix',unit);
% sj_dayprocess_cont2('07_122310', dir2, prefix, 7, 'cmperpix',unit);
% sj_dayprocess_cont2('08_122410', dir2, prefix, 8, 'cmperpix',unit);
% sj_dayprocess_cont2('09_122510', dir2, prefix, 9, 'cmperpix',unit);
% sj_dayprocess_cont2('10_122610', dir2, prefix, 10, 'cmperpix',unit);
% sj_dayprocess_cont2('11_122710', dir2, prefix, 11, 'cmperpix',unit);



% %After Spike Clustering is done
% createtetinfostruct(dir2,prefix);
% createcellinfostruct(dir2,prefix);
% riptetlist = [3,4,5,6,10,11]; 
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');


% sj_clusterdayprocess(dir, dir2, prefix, 1:8, 2:5); 




% for n=1:11
%     sj_makedayparms_dio(dir, prefix, n);
% end
%
% for n=1:11
%     sj_rippledayprocess(dir2,prefix,n);
%     %sj_remart_rippledayprocess(dir2,prefix,n);
% end 

% %min_suprathresh_duration = 0.015; %sec
% %nstd=3;
% %tetrode=3; 
% %for n=1:6
%     %sj_extractripples_nostim(dir2, prefix, n, tetrode, min_suprathresh_duration, nstd);
% %end




% %------------------------------------------------------------------------
% 
% cd '/data25/sjadhav/RippleInterruption/REe'
% prefix='REe'
% dir='/data25/sjadhav/RippleInterruption/REe/';
% dir2='/data25/sjadhav/RippleInterruption/REe_direct/';
% unit=0.45;  %cmperpixel
% unit_sl=0.17; 
% 
% epochs = [1,3,5]; type = 'sleep';
% for n=1:10
%     sj_updatetaskstruct(dir2,prefix,n,epochs,type);
% end
% for n=1:10
%     sj_updatelinposstate(dir2,prefix,n);
% end
% for n=1:10
%     sj_updatesleepcmperpix(dir2,prefix,n, unit, unit_sl);
% end
% 
% riptetlist = [3,4,6,11,12,13]; 
% for n=1:10
%     %sj_rippledayprocess(dir2,prefix,n);
%     sj_extractripples(dir2,prefix,n,riptetlist,0.015,2);
% end
% createtetinfostruct(dir2,prefix);
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');



% sj_dayprocess_cont2('01_010511', dir2, prefix, 1, 'cmperpix',unit); 
% sj_dayprocess_cont2('02_010611', dir2, prefix, 2, 'cmperpix',unit);
% sj_dayprocess_cont2('03_010711', dir2, prefix, 3, 'cmperpix',unit);
% sj_dayprocess_cont2('04_010811', dir2, prefix, 4, 'cmperpix',unit);
% sj_dayprocess_cont2('05_010911', dir2, prefix, 5, 'cmperpix',unit);
% sj_dayprocess_cont2('06_011011', dir2, prefix, 6, 'cmperpix',unit);
% sj_dayprocess_cont2('07_011111', dir2, prefix, 7, 'cmperpix',unit);
% sj_dayprocess_cont2('08_011211', dir2, prefix, 8, 'cmperpix',unit);
% sj_dayprocess_cont2('09_011311', dir2, prefix, 9, 'cmperpix',unit);
% sj_dayprocess_cont2('10_011411', dir2, prefix, 10, 'cmperpix',unit);
% disp('Dayprocess done');



% %After Spike Clustering is done and dayprocess run for spikes
% createtetinfostruct(dir2,prefix);
% createcellinfostruct(dir2,prefix);
% riptetlist = [3,4,6,11,12,13]; 
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');


% sj_clusterdayprocess(dir, dir2, prefix, 1:8, 2:5); 


% for n=1:10
%     sj_makedayparms_dio(dir, prefix, n);
% end
%
% for n=1:10
%     sj_rippledayprocess(dir2,prefix,n);
% end
 

 


% %------------------------------------------------------------------------
 
% cd '/data25/sjadhav/RippleInterruption/REf'
% prefix='REf'
% dir='/data25/sjadhav/RippleInterruption/REf/';
% dir2='/data25/sjadhav/RippleInterruption/REf_direct/';
% unit=0.45;  %cmperpixel
% unit_sl=0.17; 
% % % 

% epochs = [1,3,5]; type = 'sleep';
% for n=1:10
%     sj_updatetaskstruct(dir2,prefix,n,epochs,type);
% end
% for n=1:10
%     sj_updatelinposstate(dir2,prefix,n);
% end
% for n=1:10
%     sj_updatesleepcmperpix(dir2,prefix,n, unit, unit_sl);
% end
%
% riptetlist = [1,5,9,10,11,12]; 
% for n=1:10
%     %sj_rippledayprocess(dir2,prefix,n);
%     sj_extractripples(dir2,prefix,n,riptetlist,0.015,2);
% end
% createtetinfostruct(dir2,prefix);
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');



% sj_dayprocess_cont2('01_013111', dir2, prefix, 1, 'cmperpix',unit); 
% sj_dayprocess_cont2('02_020111', dir2, prefix, 2, 'cmperpix',unit);
% sj_dayprocess_cont2('03_020211', dir2, prefix, 3, 'cmperpix',unit);
% sj_dayprocess_cont2('04_020311', dir2, prefix, 4, 'cmperpix',unit);
% sj_dayprocess_cont2('05_020411', dir2, prefix, 5, 'cmperpix',unit);
% sj_dayprocess_cont2('06_020511', dir2, prefix, 6, 'cmperpix',unit);
% sj_dayprocess_cont2('07_020611', dir2, prefix, 7, 'cmperpix',unit);
% sj_dayprocess_cont2('08_020711', dir2, prefix, 8, 'cmperpix',unit);
% sj_dayprocess_cont2('09_020811', dir2, prefix, 9, 'cmperpix',unit);
% sj_dayprocess_cont2('10_020911', dir2, prefix, 10, 'cmperpix',unit);



% %After Spike Clustering is done - run dayprocess and then this
% createtetinfostruct(dir2,prefix);
% createcellinfostruct(dir2,prefix);
% riptetlist = [1,5,9,10,11,12]; 
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');


% sj_clusterdayprocess(dir, dir2, prefix, 1:8, 2:5); 


% for n=1:10
%     sj_makedayparms_dio(dir, prefix, n);
% end
% 
% for n=1:10
%     sj_rippledayprocess(dir2,prefix,n);
% end
%
% for n=1:10
%     thetadayprocess(dir2,prefix,n);
% end





% %------------------------------------------------------------------------
% %------------------------------------------------------------------------
%% Control Stimulation Group
% %------------------------------------------------------------------------


% % 
% cd '/data25/sjadhav/RippleInterruption/RCd/'
% prefix='RCd'
% dir='/data25/sjadhav/RippleInterruption/RCd/';
% dir2='/data25/sjadhav/RippleInterruption/RCd_direct/';
% unit=0.45;  %cmperpixel
% unit_sl=0.17; 
% % % 

% epochs = [1,3,5]; type = 'sleep';
% for n=1:10
%     sj_updatetaskstruct(dir2,prefix,n,epochs,type);
% end
% for n=1:10
%     sj_updatelinposstate(dir2,prefix,n);
% end
% for n=1:10
%     sj_updatesleepcmperpix(dir2,prefix,n, unit, unit_sl);
% end
% 
% riptetlist = [1,2,3,4,5,6]; 
% for n=1:10
%     %sj_rippledayprocess(dir2,prefix,n);
%     sj_extractripples(dir2,prefix,n,riptetlist,0.015,2);
% end
% createtetinfostruct(dir2,prefix);
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');



% sj_dayprocess_cont2('01_042511', dir2, prefix, 1, 'cmperpix',unit); 
% sj_dayprocess_cont2('02_042611', dir2, prefix, 2, 'cmperpix',unit);
% sj_dayprocess_cont2('03_042711', dir2, prefix, 3, 'cmperpix',unit);
% sj_dayprocess_cont2('04_042811', dir2, prefix, 4, 'cmperpix',unit);
% sj_dayprocess_cont2('05_042911', dir2, prefix, 5, 'cmperpix',unit);
% sj_dayprocess_cont2('06_043011', dir2, prefix, 6, 'cmperpix',unit);
% sj_dayprocess_cont2('07_050111', dir2, prefix, 7, 'cmperpix',unit);
% sj_dayprocess_cont2('08_050211', dir2, prefix, 8, 'cmperpix',unit);
% sj_dayprocess_cont2('09_050311', dir2, prefix, 9, 'cmperpix',unit);
% sj_dayprocess_cont2('10_050411', dir2, prefix, 10, 'cmperpix',unit);



% %After Spike Clustering is done
% createtetinfostruct(dir2,prefix);
% createcellinfostruct(dir2,prefix);
% riptetlist = [1,2,3,4,5,6]; 
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');

% sj_clusterdayprocess(dir, dir2, prefix, [1:8], [2:5]); 



% for n=1:10
%     sj_makedayparms_dio(dir, prefix, n);
% end
%
% for n=1:10
%     sj_rippledayprocess(dir2,prefix,n);
% end




%-----------------------------------------------------------------------


% % 
% cd '/data25/sjadhav/RippleInterruption/RCc'
% prefix='RCc'
% dir='/data25/sjadhav/RippleInterruption/RCc';
% dir2='/data25/sjadhav/RippleInterruption/RCc_direct/';
% unit=0.45;  %cmperpixel
% unit_sl=0.17; 
% 
% epochs = [1,3,5]; type = 'sleep';
% for n=1:10
%     sj_updatetaskstruct(dir2,prefix,n,epochs,type);
% end
% for n=1:10
%     sj_updatelinposstate(dir2,prefix,n);
% end
% for n=1:10
%     sj_updatesleepcmperpix(dir2,prefix,n, unit, unit_sl);
% end
%
% riptetlist = [3,4,5,6,11,13]; 
% for n=1:10
%     %sj_rippledayprocess(dir2,prefix,n);
%     %sj_extractripples(dir2,prefix,n,riptetlist,0.015,2);
%     sj_getandsave_ripacrosstet(prefix,n);
% end
% createtetinfostruct(dir2,prefix);
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');



% sj_dayprocess_cont2('01_040311', dir2, prefix, 1, 'cmperpix',unit); 
% sj_dayprocess_cont2('02_040411', dir2, prefix, 2, 'cmperpix',unit);
% sj_dayprocess_cont2('03_040511', dir2, prefix, 3, 'cmperpix',unit);
% sj_dayprocess_cont2('04_040611', dir2, prefix, 4, 'cmperpix',unit);
% sj_dayprocess_cont2('05_040711', dir2, prefix, 5, 'cmperpix',unit);
% sj_dayprocess_cont2('06_040811', dir2, prefix, 6, 'cmperpix',unit);
% sj_dayprocess_cont2('07_040911', dir2, prefix, 7, 'cmperpix',unit);
% sj_dayprocess_cont2('08_041011', dir2, prefix, 8, 'cmperpix',unit);
% sj_dayprocess_cont2('09_041111', dir2, prefix, 9, 'cmperpix',unit);
% sj_dayprocess_cont2('10_041211', dir2, prefix, 10, 'cmperpix',unit);



% %After Spike Clustering is done
% createtetinfostruct(dir2,prefix);
% createcellinfostruct(dir2,prefix);
% riptetlist = [3,4,5,6,11,13]; 
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');

% sj_clusterdayprocess(dir, dir2, prefix, [1:8], [2:5]); 

% Day 6 is empty. In cellinfo, make each epoch empty
% cd(dir2);
% load RCccellinfo
% for ep=1:5
%     cellinfo{6}{ep}=[];
% end
% save RCccellinfo cellinfo;





% for n=1:10
%     sj_makedayparms_dio(dir, prefix, n);
% end
%
% for n=10:10
%     sj_rippledayprocess(dir2,prefix,n);
% end
 




%-----------------------------------------------------------------------

% 
% cd '/data25/sjadhav/RippleInterruption/RCb'
% prefix='RCb'
% dir='/data25/sjadhav/RippleInterruption/RCb/';
% dir2='/data25/sjadhav/RippleInterruption/RCb_direct/';
% unit=0.45;  %cmperpixel
% unit_sl=0.17; 
% 
% epochs = [1,3,5]; type = 'sleep';
% for n=1:10
%     sj_updatetaskstruct(dir2,prefix,n,epochs,type);
% end
% for n=1:10
%     sj_updatelinposstate(dir2,prefix,n);
% end
% for n=1:10
%     sj_updatesleepcmperpix(dir2,prefix,n, unit, unit_sl);
% end
% 
% riptetlist = [3,4,9,10,11,12]; 
% for n=1:10
%     %sj_rippledayprocess(dir2,prefix,n);
%     n
%     %sj_extractripples(dir2,prefix,n,riptetlist,0.015,2);
%     sj_getandsave_ripacrosstet(prefix,n);
% end
% createtetinfostruct(dir2,prefix);
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');



% sj_dayprocess_cont2('01_021811', dir2, prefix, 1, 'cmperpix',unit); 
% sj_dayprocess_cont2('02_021911', dir2, prefix, 2, 'cmperpix',unit);
% sj_dayprocess_cont2('03_022011', dir2, prefix, 3, 'cmperpix',unit);
% sj_dayprocess_cont2('04_022111', dir2, prefix, 4, 'cmperpix',unit);
% %% sj_dayprocess_cont2('05_022211', dir2, prefix, 5, 'cmperpix',unit); %Error this day
% sj_dayprocess_cont2('06_022311', dir2, prefix, 5, 'cmperpix',unit);
% sj_dayprocess_cont2('07_022411', dir2, prefix, 6, 'cmperpix',unit);
% sj_dayprocess_cont2('08_022511', dir2, prefix, 7, 'cmperpix',unit);
% sj_dayprocess_cont2('09_022611', dir2, prefix, 8, 'cmperpix',unit);
% sj_dayprocess_cont2('10_022711', dir2, prefix, 9, 'cmperpix',unit);
% sj_dayprocess_cont2('11_022811', dir2, prefix, 10, 'cmperpix',unit);

% sj_clusterdayprocess(dir, dir2, prefix, [1:8], [2:5]); 


% %After Spike Clustering is done
% createtetinfostruct(dir2,prefix);
% createcellinfostruct(dir2,prefix);
% riptetlist = [3,4,9,10,11,12]; 
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');





% for n=1:10
%     sj_makedayparms_dio(dir, prefix, n);
% end
% 
% for n=10:10
%     sj_rippledayprocess(dir2,prefix,n);
% end
 




%-----------------------------------------------------------------------

% 
% cd '/data25/sjadhav/RippleInterruption/RCa'
% prefix='RCa'
% dir='/data25/sjadhav/RippleInterruption/RCa/';
% dir2='/data25/sjadhav/RippleInterruption/RCa_direct/';
% unit=0.46;  %cmperpixel
% unit_sl=0.15;

% epochs = [1,3,5]; type = 'sleep';
% for n=1:8
%     sj_updatetaskstruct(dir2,prefix,n,epochs,type);
% end
% for n=1:8
%     sj_updatelinposstate(dir2,prefix,n);
% end
% for n=1:8
%     sj_updatesleepcmperpix(dir2,prefix,n, unit, unit_sl);
% end
%
% riptetlist = [2,3,4,6,9]; 
% for n=1:8
%     %sj_rippledayprocess(dir2,prefix,n);
%     sj_extractripples(dir2,prefix,n,riptetlist,0.015,2);
% end
% createtetinfostruct(dir2,prefix);
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');


% 
%  sj_dayprocess('01_052010', dir2, prefix, 1, 'cmperpix',unit); 
%  sj_dayprocess('02_052110', dir2, prefix, 2, 'cmperpix',unit);
%  sj_dayprocess('03_052210', dir2, prefix, 3, 'cmperpix',unit);
%  sj_dayprocess('04_052310', dir2, prefix, 4, 'cmperpix',unit);
%  sj_dayprocess('05_052410', dir2, prefix, 5, 'cmperpix',unit);
%  sj_dayprocess('06_052510', dir2, prefix, 6, 'cmperpix',unit);
%  sj_dayprocess('07_052610', dir2, prefix, 7, 'cmperpix',unit);
%  sj_dayprocess('08_052710', dir2, prefix, 8, 'cmperpix',unit);



% %After Spike Clustering is done
% createtetinfostruct(dir2,prefix);
% createcellinfostruct(dir2,prefix);
% riptetlist = [2,3,4,6,9]; 
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');

% %Day 3-8 are empty. In cellinfo, make each epoch empty
% cd(dir2);
% load RCacellinfo
% for d=3:8
%     for ep=1:5
%         cellinfo{d}{ep}=[];
%     end
% end
% save RCacellinfo cellinfo;


%sj_clusterdayprocess(dir, dir2, prefix, [1:8], [2:5]); 




% for n=1:8
%     sj_makedayparms_dio(dir, prefix, n);
% end
%
% for n=1:8
%     sj_rippledayprocess(dir2,prefix,n);
% end





%-----------------------------------------------------------------------

% cd '/data25/sjadhav/RippleInterruption/RE1'
% prefix='RE1'
% dir='/data25/sjadhav/RippleInterruption/RE1/';
% dir2='/data25/sjadhav/RippleInterruption/RE1_direct/';
% unit=0.67;  %cmperpixel
% unit_sl=0.27; 
% 

% epochs = [1,3,5]; type = 'sleep';
% for n=1:8
%     sj_updatetaskstruct(dir2,prefix,n,epochs,type);
% end
% for n=1:8
%     sj_updatelinposstate(dir2,prefix,n);
% end
% for n=1:8
%     sj_updatesleepcmperpix(dir2,prefix,n, unit, unit_sl);
% end
% 
% riptetlist = [1 5 7]; 
% for n=1:8
%     %sj_rippledayprocess(dir2,prefix,n);
%     sj_extractripples(dir2,prefix,n,riptetlist,0.015,2);
% end
% createtetinfostruct(dir2,prefix);
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');


% sj_dayprocess('01_022710', dir2, prefix, 1, 'cmperpix',unit); 
% sj_dayprocess('02_022810', dir2, prefix, 2, 'cmperpix',unit);
% sj_dayprocess('03_030110', dir2, prefix, 3, 'cmperpix',unit);
% sj_dayprocess('04_030210', dir2, prefix, 4, 'cmperpix',unit);
% sj_dayprocess('05_030310', dir2, prefix, 5, 'cmperpix',unit);
% sj_dayprocess('06_030410', dir2, prefix, 6, 'cmperpix',unit);
% sj_dayprocess('07_030510', dir2, prefix, 7, 'cmperpix',unit);
% sj_dayprocess('08_030610', dir2, prefix, 8, 'cmperpix',unit);
% sj_dayprocess('09_030710', dir2, prefix, 9, 'cmperpix',unit);
% sj_dayprocess('10_030810', dir2, prefix, 10, 'cmperpix',unit);

% %After Spike Clustering is done
% createtetinfostruct(dir2,prefix);
% createcellinfostruct(dir2,prefix);
% riptetlist = [1,5,7]; 
% sj_addtetrodedescription(dir2,prefix,riptetlist,'riptet');
% 

% sj_clusterdayprocess(dir, dir2, prefix, [1:8], [2:5]); 

% for n=1:10
%     sj_makedayparms_dio(dir, prefix, n);
% end
%
% for n=1:10
%     sj_rippledayprocess(dir2,prefix,n);
% end





%------------------------------------------------------------------------


% %------------------------------------------------------------------------
% %------------------------------------------------------------------------
%% Normal Stimulation Group
% %------------------------------------------------------------------------

% % 
% cd '/data25/sjadhav/RippleInterruption/RNa'
% 
% prefix='RNa';
% dir='/data25/sjadhav/RippleInterruption/RNa';
% dir2='/data25/sjadhav/RippleInterruption/RNa_direct/';
% unit=0.45;  %cmperpixel
% % unit_sl=0.17;


% sj_dayprocess_cont2('01_041311', dir2, prefix, 1, 'cmperpix',unit,'processeeg',0); 
% sj_dayprocess_cont2('02_041411', dir2, prefix, 2, 'cmperpix',unit,'processeeg',0);
% sj_dayprocess_cont2('03_041511', dir2, prefix, 3, 'cmperpix',unit,'processeeg',0);
% sj_dayprocess_cont2('04_041611', dir2, prefix, 4, 'cmperpix',unit,'processeeg',0);
% sj_dayprocess_cont2('05_041711', dir2, prefix, 5, 'cmperpix',unit,'processeeg',0);
% sj_dayprocess_cont2('06_041811', dir2, prefix, 6, 'cmperpix',unit,'processeeg',0);
% sj_dayprocess_cont2('07_041911', dir2, prefix, 7, 'cmperpix',unit,'processeeg',0);
% sj_dayprocess_cont2('08_042011', dir2, prefix, 8, 'cmperpix',unit,'processeeg',0);

%------------------------------------------------------------------------

% cd '/data25/sjadhav/RippleInterruption/RNb'
% 
% prefix='RNb';
% dir='/data25/sjadhav/RippleInterruption/RNb';
% dir2='/data25/sjadhav/RippleInterruption/RNb_direct/';
% unit=0.45;  %cmperpixel
% % unit_sl=0.17;

% makedayparms('08_051711'); 

% sj_dayprocess_cont2('01_051011', dir2, prefix, 1, 'cmperpix',unit,'processeeg',0); 
% sj_dayprocess_cont2('02_051111', dir2, prefix, 2, 'cmperpix',unit,'processeeg',0);
% sj_dayprocess_cont2('03_051211', dir2, prefix, 3, 'cmperpix',unit,'processeeg',0);
% sj_dayprocess_cont2('04_051311', dir2, prefix, 4, 'cmperpix',unit,'processeeg',0);
% sj_dayprocess_cont2('05_051411', dir2, prefix, 5, 'cmperpix',unit,'processeeg',0);
% sj_dayprocess_cont2('06_051511', dir2, prefix, 6, 'cmperpix',unit,'processeeg',0);
% sj_dayprocess_cont2('07_051611', dir2, prefix, 7, 'cmperpix',unit,'processeeg',0);
% sj_dayprocess_cont2('08_051711', dir2, prefix, 8, 'cmperpix',unit,'processeeg',0);


%------------------------------------------------------------------------

% cd '/data25/sjadhav/RippleInterruption/RNc'
% 
% prefix='RNc';
% dir='/data25/sjadhav/RippleInterruption/RNc';
% dir2='/data25/sjadhav/RippleInterruption/RNc_direct/';
% unit=0.45;  %cmperpixel
% % % unit_sl=0.17;
% 
%  sj_dayprocess_cont2('rn_m1', dir2, prefix, 1, 'cmperpix',unit,'processeeg',0); 
%  sj_dayprocess_cont2('rn_m2', dir2, prefix, 2, 'cmperpix',unit,'processeeg',0);
%  sj_dayprocess_cont2('rn_m3', dir2, prefix, 3, 'cmperpix',unit,'processeeg',0);
%  sj_dayprocess_cont2('rn_m4', dir2, prefix, 4, 'cmperpix',unit,'processeeg',0);
%  sj_dayprocess_cont2('rn_m5', dir2, prefix, 5, 'cmperpix',unit,'processeeg',0);
%  sj_dayprocess_cont2('rn_m6', dir2, prefix, 6, 'cmperpix',unit,'processeeg',0);
%  sj_dayprocess_cont2('rn_m7', dir2, prefix, 7, 'cmperpix',unit,'processeeg',0);
%  sj_dayprocess_cont2('rn_m8', dir2, prefix, 8, 'cmperpix',unit,'processeeg',0);

% cd(dir2)
% for n=1:8
%      sj_lineardayprocess(dir2,prefix,n);
% end

%------------------------------------------------------------------------

% cd '/data25/sjadhav/RippleInterruption/RNd'
% 
% prefix='RNd';
% dir='/data25/sjadhav/RippleInterruption/RNd';
% dir2='/data25/sjadhav/RippleInterruption/RNd_direct/';
% unit=0.45;  %cmperpixel
% % % unit_sl=0.17;
% 
%  sj_dayprocess_cont2('rn_n1', dir2, prefix, 1, 'cmperpix',unit,'processeeg',0); 
%  sj_dayprocess_cont2('rn_n2', dir2, prefix, 2, 'cmperpix',unit,'processeeg',0);
%  sj_dayprocess_cont2('rn_n3', dir2, prefix, 3, 'cmperpix',unit,'processeeg',0);
%  sj_dayprocess_cont2('rn_n4', dir2, prefix, 4, 'cmperpix',unit,'processeeg',0);
%  sj_dayprocess_cont2('rn_n5', dir2, prefix, 5, 'cmperpix',unit,'processeeg',0);
%  sj_dayprocess_cont2('rn_n6', dir2, prefix, 6, 'cmperpix',unit,'processeeg',0);
%  sj_dayprocess_cont2('rn_n7', dir2, prefix, 7, 'cmperpix',unit,'processeeg',0);
%  sj_dayprocess_cont2('rn_n8', dir2, prefix, 8, 'cmperpix',unit,'processeeg',0);

% cd(dir2)
% for n=1:8
%      sj_lineardayprocess(dir2,prefix,n);
% end
% 
 
 

%%%%%%%%%%%%% Linear Day Process %%%%%%%%%%%%%%%%%%%%%

% prefix='REc';
% dir2='/data25/sjadhav/RippleInterruption/REc_direct/';
% cd(dir2)
% for n=1:8
%     sj_lineardayprocess(dir2,prefix,n);
% end

% prefix='REd';
% dir2='/data25/sjadhav/RippleInterruption/REd_direct/';
% cd(dir2)
% for n=1:10
%     sj_lineardayprocess(dir2,prefix,n);
% end

% prefix='REe';
% dir2='/data25/sjadhav/RippleInterruption/REe_direct/';
% cd(dir2)
% for n=1:10
%     sj_lineardayprocess(dir2,prefix,n);
% end

% prefix='REf';
% dir2='/data25/sjadhav/RippleInterruption/REf_direct/';
% cd(dir2)
% for n=1:10
%     sj_lineardayprocess(dir2,prefix,n);
% end

% -------------------------------------------

%  CONTROL GROUP

% prefix='RCa';
% dir2='/data25/sjadhav/RippleInterruption/RCa_direct/';
% cd(dir2)
% for n=1:8
%     sj_lineardayprocess(dir2,prefix,n);
% end

% prefix='RCb';
% dir2='/data25/sjadhav/RippleInterruption/RCb_direct/';
% cd(dir2)
% for n=1:10
%     sj_lineardayprocess(dir2,prefix,n);
% end

% prefix='RCc';
% dir2='/data25/sjadhav/RippleInterruption/RCc_direct/';
% cd(dir2)
% for n=1:10
%     sj_lineardayprocess(dir2,prefix,n);
% end

% prefix='RCd';
% dir2='/data25/sjadhav/RippleInterruption/RCd_direct/';
% cd(dir2)
% for n=1:10
%     sj_lineardayprocess(dir2,prefix,n);
% end

% prefix='RE1';
% dir2='/data25/sjadhav/RippleInterruption/RE1_direct/';
% cd(dir2)
% for n=1:8
%     sj_lineardayprocess(dir2,prefix,n);
% end

% -------------------------------------------

%  OTH Control Data GROUP

%---
% Atleast Eight Days or More

% prefix='Cor';
% dir2='/data25/sjadhav/RippleInterruption/OthControlData/Cor/';
% cd(dir2)
% for n=4:9
%     sj_lineardayprocess_nodirn(dir2,prefix,n);
% end
% 
% prefix='Fiv';
% dir2='/data25/sjadhav/RippleInterruption/OthControlData/Fiv/';
% cd(dir2)
% for n=4:9
%     sj_lineardayprocess_nodirn(dir2,prefix,n);
% end
% 
% prefix='fra';
% dir2='/data25/sjadhav/RippleInterruption/OthControlData/Fra/';
% cd(dir2)
% for n=1:12
%     sj_lineardayprocess_nodirn(dir2,prefix,n);
% end
% 
% prefix='Sev';
% dir2='/data25/sjadhav/RippleInterruption/OthControlData/Sev/';
% cd(dir2)
% for n=1:9
%     sj_lineardayprocess_nodirn(dir2,prefix,n);
% end
% 
% %---
% % Less Than Eight Days
% 
% prefix='con';
% dir2='/data25/sjadhav/RippleInterruption/OthControlData/Con/';
% cd(dir2)
% for n=1:6
%     sj_lineardayprocess_nodirn(dir2,prefix,n);
% end
% 
%% Dudley has problems with day 6 - No completed trajectories
% % prefix='dud';
% % dir2='/data25/sjadhav/RippleInterruption/OthControlData/Dud/';
% % cd(dir2)
% % for n=1:5
% %     sj_lineardayprocess_nodirn(dir2,prefix,n);
% % end

% prefix='Eig';
% dir2='/data25/sjadhav/RippleInterruption/OthControlData/Eig/';
% cd(dir2)
% for n=1:7
%     sj_lineardayprocess_nodirn(dir2,prefix,n);
% end
% 
% prefix='Six';
% dir2='/data25/sjadhav/RippleInterruption/OthControlData/Six/';
% cd(dir2)
% for n=1:5
%     sj_lineardayprocess_nodirn(dir2,prefix,n);
% end

% prefix='ten';
% dir2='/data25/sjadhav/RippleInterruption/OthControlData/Ten';
% cd(dir2)
% for n=1:7
%     sj_lineardayprocess_nodirn(dir2,prefix,n);
% end

