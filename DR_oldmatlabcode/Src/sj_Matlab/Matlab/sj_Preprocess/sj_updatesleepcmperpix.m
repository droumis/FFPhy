
function [pos] = sj_updatesleepcmperpix(animdirect,prefix,day, runcmperpix, newcmperpix)

% Shantanu Aug2011
% Similar to sj_updatelinposstate and sj_updatetaskstruct
% If sleep box cmperpix is wrongly assigned in the dayprocess file, update
% it with new cmperpix 

% 2 options: 
% 1) Get rawpos and do from scratch, OR
% 2) Just use ratio of wrong cmperpix and newcmperpix to update pos

% Use method 2

% pos file
posfile = sprintf('%s/%spos%02d.mat', animdirect, prefix, day);
load(posfile);

% Task file
taskfile = sprintf('%s/%stask%02d.mat', animdirect, prefix, day);
load(taskfile);

% For method 1
% % rawpos file
% rawposfile = sprintf('%s/%srawpos%02d.mat', animdirect, prefix, day);
% load(rawposfile);
% % Parameters similar to dayprocess
% diodepos = 0.5;
% reversex = 0;
% reversey = 0;

posfilt = gaussian(30*0.5, 60); % Similar to dayprocess

% Across all sleep epochs
for i = 1:length(task{day})
    if ((~isempty(task{day}{i})) && (strcmp(task{day}{i}.type,'sleep')) )
        
        disp(['Day ',num2str(day), ', Epoch ',num2str(i)])
        
        % Need to check if it is updated. Either compare the cmperpixel
        % field to runcmperpixel and newcmperpixel, or easier - make an
        % updateflag field in pos{day}{ep} or pos{day}{ep}.arg for the sleep epoch
        i;
        if ~isfield(pos{day}{i},'updateflag')
            % If updateflag doesnt exist, update with newcmperpixel
            pos{day}{i}.data(:,2:3) = pos{day}{i}.data(:,2:3)*(newcmperpix./runcmperpix);
            pos{day}{i}.cmperpixel = newcmperpix;
            % Update velocity: addvelocity will replace the original velocity
            % column with the new one by checking for the token "vel" in
            % pos{day}{epoch}.fields
            pos{day}{i} = addvelocity(pos{day}{i}, posfilt);
            % Add updateflag
            pos{day}{i}.updateflag=1;
        end
    end
end

% Save updated pos file
save(posfile,'pos');



