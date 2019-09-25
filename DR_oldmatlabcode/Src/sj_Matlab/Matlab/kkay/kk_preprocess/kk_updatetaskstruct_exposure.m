
function kk_updatetaskstruct_exposure(animdirect,prefix,dayepochexp)

% assigns exposure value to .exposure field to selected days and epochs
% ex. assigning 7th, 8th, and 9th, exposures for animal Ben for day 3
% epochs 2, 4, and 6
    % kk_updatetaskstruct_exposure('/data99/kkay/Ben/','ben',[3 2 7 ; 3 4 8 ; 3 6 9])
    
    % dayepochexp is a triplet of day, epoch, and exposure #

for p=1:size(dayepochexp,1)
    day=dayepochexp(p,1);
    epoch=dayepochexp(p,2);
    taskfile = sprintf('%s/%stask%02d.mat', animdirect, prefix, day);
    load(taskfile);
    task{day}{epoch}.exposure = dayepochexp(p,3);
    % Save updated linpos file
    save(taskfile,'task');
end

