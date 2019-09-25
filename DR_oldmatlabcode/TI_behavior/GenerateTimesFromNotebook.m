% 643830
% 13108176

ranges(1,:) = [0 13108177]; %start at 0 and go 1 unit past the last dio from the .dio file
ranges(2,:) = [0 13108177]; %start at 0 and go 1 unit past the last dio from the .dio file
names = {};
names{1} = '1>> All points';
names{2} = '2  Run1  00:01:00-00:22:15';
save('/data19/droumis/TransInf/t03today/times.mat', 'ranges', 'names');
 
diodayprocess('/data19/droumis/TransInf/t03today','/data19/droumis/TransInf/t03today',1,1)