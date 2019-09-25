function addtetrodethetamod(animdirect,fileprefix,tetrodes)
% This function adds the mean theta envelope to tetinfo.

load([animdirect, fileprefix,'tetinfo']);
o = cellfetch(tetinfo,'area');
target = o.index(find(ismember(o.index(:,3),tetrodes)),:);
for i = 1:size(target,1)
    if target(i,3) >= 10 
        load([animdirect,'EEG/',fileprefix,'theta0',num2str(target(i,1)),'-',num2str(target(i,2)),'-',num2str(target(i,3))]);
    elseif target(i,3) <10
        load([animdirect,'EEG/',fileprefix,'theta0',num2str(target(i,1)),'-',num2str(target(i,2)),'-0',num2str(target(i,3))]);
    end
    thetamod = mean(theta{target(i,1)}{target(i,2)}{target(i,3)}.data(:,3));
    tetinfo{target(i,1)}{target(i,2)}{target(i,3)}.thetamod = thetamod;
end

save([animdirect, fileprefix,'tetinfo'], 'tetinfo')

end

 