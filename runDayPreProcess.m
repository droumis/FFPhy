

function runDayPreProcess()
% run an animal's preprocess log and call this at the bottom
tic
w = evalin('caller', 'whos');
for a = 1:length(w)
    dp.(w(a).name) = evalin('caller', w(a).name);
    pause(.01)
end
for i = 1:length(dp.days)
    fprintf('date %d day %d \n', dp.dates(i), dp.days(i));
    Trodes_dayprocess(dp.animal, dp.dates(i), dp.days(i), 'dp', dp);
end
fprintf('nextcoffeebreak = %0.0f seconds \n',toc);
end