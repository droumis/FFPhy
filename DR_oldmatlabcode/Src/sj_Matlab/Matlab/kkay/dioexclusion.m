

% run ms_diodayprocess
for d=2:8
    daydir=sprintf('dave0%d',d);
    ms_diodayprocess(daydir,'/data12/kkay/Dav/','dav',d)
end

% run Shantanu's parsing script (open it to fix the well #s)

open sj_findwellsfromdio1_Egypt

%% Chapati
[wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt('/datatmp/kkay/Cha','cha',1,[3 4 7 8],22:24);
[wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt('/datatmp/kkay/Cha','cha',2,[3 4 7 8],22:24);
[wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt('/datatmp/kkay/Cha','cha',3,[2 4],22:24);
[wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt('/datatmp/kkay/Cha','cha',4,[2 4],22:24);
[wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt('/datatmp/kkay/Cha','cha',5,[2 4],22:24);
[wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt('/datatmp/kkay/Cha','cha',6,[2 4 6],22:24);
[wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt('/datatmp/kkay/Cha','cha',7,[2 5],22:24);
[wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt('/datatmp/kkay/Cha','cha',8,[2 4 6],22:24);
[wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt('/datatmp/kkay/Cha','cha',9,[2 4 6],22:24);

%% Dave
[wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt_lin('/data12/kkay/Dav','dav',1,[2 4 6],22:24);
[wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt('/data12/kkay/Dav','dav',2,[2 4 6],22:24);
[wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt('/data12/kkay/Dav','dav',3,[2 4 6],22:24);
[wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt('/data12/kkay/Dav','dav',4,[2 4 6],22:24);
[wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt('/data12/kkay/Dav','dav',5,[3 5],22:24);
[wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt('/data12/kkay/Dav','dav',6,[2 4 6],22:24);
[wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt('/data12/kkay/Dav','dav',7,[2 3 5 7 9],22:24);
[wellsdio, rewardinfo] = sj_findwellsfromdio1_Egypt('/data12/kkay/Dav','dav',8,[2 3 5],22:24);




% run generate_rewarderror_times_Egypt to distinguish rewarded vs. error
open generate_rewarderror_times_Egypt

% verify with position plot
open plot_rewardtimes_position