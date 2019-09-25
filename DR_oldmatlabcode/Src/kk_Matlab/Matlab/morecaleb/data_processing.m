
if 0
%% These write to times.mat file directly
[names, ranges] = generateTimesFileFromCont();
[names, ranges] = modifyTimes(6,'run 1','0:23:31','0:37:20');
[names, ranges] = modifyTimes(10,'run 2 (novel 1)','1:24:30','1:39:41');

% Find cmperpix by running posrecon.run, extracting a good frame,
% then clicking on the outside corners of the track. Compare the
% diagonals to what we know they should be for 30 inch lengths.
centimeters_per_pixel = cell(size(ranges,1)-1,1);
centimeters_per_pixel{6} = 108/343.8; % 108 = sqrt(2)*2.54.*30; 343.8 is diagonal pixels mean
centimeters_per_pixel{10} = 108/414; % 108 = sqrt(2)*2.54.*30; 421/408 is diagonal pixels mean
end

daynum = 2;

stimdio = createstimstruct('pin',[16,17],'pindescriptions',{'ca3','ec'});
save stimdio stimdio
% (old) NOTE- pin "32" appears to correspond to "A" = CA3 = pin 16
% (old) NOTE- pin "31" appears to correspond to "B" = EC = pin 17

if 0
for ep = 1:length(ranges)-1
  task{ep}.epoch = ep;
  task{ep}.type = 'sleep';
  task{ep}.day = daynum;
end
task{6}.type = 'run';
task{10}.type = 'run';

task{6}.exposure = inf;
task{10}.exposure = 1;

save task task;

% Copy task and stim structures to processed directory
predayprocess(pwd,'/data27/data/monster/Uri','uri',daynum, ...
  'generateTetinfo', 0);
end

newdayprocess('/data27/data/monster/Uri','uri',daynum, ...
  'process_pos', 0, 'process_stim', 0, 'overwrite', 0, ...
  'process_stimeeg', 0, 'process_longstimeeg', 1, ...
  'longstimeeg_windowLength', 4, 'longstimeeg_winOffset', 3);
  % 'centimeters_per_pixel', centimeters_per_pixel);


