
sj_dayprocess('01_102909_2','/data25/sjadhav/SJStimC_8arm','sjc',1,'cmperpix',0.75);


% sj_dayprocess('01_102909','/data25/sjadhav/SJStimC_direct','sjc',1,'cmperpix',0.75);
% sj_dayprocess('02_103009','/data25/sjadhav/SJStimC_direct','sjc',2,'cmperpix',0.75);
% sj_dayprocess('03_103109','/data25/sjadhav/SJStimC_direct','sjc',3,'cmperpix',0.75);
% sj_dayprocess('04_110109','/data25/sjadhav/SJStimC_direct','sjc',4,'cmperpix',0.75);
% sj_dayprocess('05_110209','/data25/sjadhav/SJStimC_direct','sjc',5,'cmperpix',0.75);
% sj_dayprocess('06_110309','/data25/sjadhav/SJStimC_direct','sjc',6,'cmperpix',0.75);
% sj_dayprocess('07_110409','/data25/sjadhav/SJStimC_direct','sjc',7,'cmperpix',0.75);

 
% lineardayprocess('/data25/sjadhav/SJStimC_direct/','sjc',1,'welldist',15);
% lineardayprocess('/data25/sjadhav/SJStimC_direct/','sjc',2,'welldist',15);
% lineardayprocess('/data25/sjadhav/SJStimC_direct/','sjc',3,'welldist',15);
% lineardayprocess('/data25/sjadhav/SJStimC_direct/','sjc',4,'welldist',15);
% lineardayprocess('/data25/sjadhav/SJStimC_direct/','sjc',5,'welldist',15);
% lineardayprocess('/data25/sjadhav/SJStimC_direct/','sjc',6,'welldist',15);
% lineardayprocess('/data25/sjadhav/SJStimC_direct/','sjc',7,'welldist',15);

% sj_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',3,'daytetlist',[3 2; 3 5; 3 6]);
% sj_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',4,'daytetlist',[4 2; 4 5; 4 6]);
% sj_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',5,'daytetlist',[5 2; 5 5; 5 6]);
% sj_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',6,'daytetlist',[6 2; 6 5; 6 6]);
% sj_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',7,'daytetlist',[7 2; 7 5; 7 6]);
% sj_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[1 2],'daytetlist',[1 2; 1 5; 1 6; 2 2; 2 5; 2 6]);


% sj_diodayprocess('07_110409','/data25/sjadhav/SJStimC_direct','sjc',7);
% sj_diodayprocess('06_110309','/data25/sjadhav/SJStimC_direct','sjc',6);
% sj_diodayprocess('05_110209','/data25/sjadhav/SJStimC_direct','sjc',5);
% sj_diodayprocess('04_110109','/data25/sjadhav/SJStimC_direct','sjc',4);
% sj_diodayprocess('03_103109','/data25/sjadhav/SJStimC_direct','sjc',3);
% sj_diodayprocess('02_103009','/data25/sjadhav/SJStimC_direct','sjc',2);
% sj_diodayprocess('01_102909','/data25/sjadhav/SJStimC_direct','sjc',1);
% 
% 
% 
% sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[3 4],'daytetlist',[3 2; 3 6; 4 2; 4 6]);
% sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[5 6 7],'daytetlist',[5 2; 5 6; 6 2; 6 6; 7 2; 7 6]);
% sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[1 2],'daytetlist',[1 2; 1 6; 2 2; 2 6]);


sj_extractripples_nostim('/data25/sjadhav/SJStimC_direct', 'sjc', 4, 6, 0.015, 3);






















%% 07 Jan

sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[3 4],'daytetlist',[3 2; 3 6; 4 2; 4 6]);

sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[5 6 7],'daytetlist',[5 2; 5 6; 6 2; 6 6; 7 2; 7 6]);

sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[1 2],'daytetlist',[1 2; 1 6; 2 2; 2 6]);


sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[ 7],'daytetlist',[7 2; 7 6]);













% sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[4],'daytetlist',[4, 6]);
% sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[4],'daytetlist',[4, 2]);
% 
% sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[3],'daytetlist',[3, 6]);
% sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[3],'daytetlist',[3, 2]);
% 
% sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[5],'daytetlist',[5, 6]);
% sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[6],'daytetlist',[6, 6]);
% sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[7],'daytetlist',[7, 6]);
% sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[5],'daytetlist',[5, 2]);
% sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[6],'daytetlist',[6, 2]);
% sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[7],'daytetlist',[7, 2]);
% 
% 
% 
% sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[1],'daytetlist',[1, 6]);
% sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[2],'daytetlist',[2, 6]);
% 
% sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[1],'daytetlist',[1, 2]);
% sj_remart_rippledayprocess('/data25/sjadhav/SJStimC_direct/','sjc',[2],'daytetlist',[2, 2]);
