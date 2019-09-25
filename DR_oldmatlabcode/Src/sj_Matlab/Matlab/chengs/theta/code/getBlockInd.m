function getBlockInd
%function getBlocktimes
% figure out times, at which one pass ends

global fmaux behavdata
BlockInd= {};
%OnTraj= {};
oldd= -1; olde= -1; 
[d,e,t,c,ncells]= startCellList;
while ~isempty(d)
    if(d~=oldd) | (e~= olde)
        oldd= d; olde= e; 

        loadVar(fmaux.data2dir, 'behavdata', d);
        data= behavdata{d}{e};

%        figure; plot(data.traj); hold on

        T= size(data.time,1);
        for it=0:3
            ind= find(data.traj== it);
            di= diff(ind);
            j=[0; find(di> 1000); length(ind)]; % return to traj after > 2sec elsewhere
            nB= length(j)-1;
            BlockInd{d}{e}{it+1}= [ind(j(1:nB)+1),ind(j(2:end))];
%            for iB=1:nB;
%                OnTraj{d}{e}{it+1}{iB}= ind(j(iB)+1:j(iB+1));
%            end

%            plot(BlockInd{d}{e}{it+1}(:,1), it, 'sg')
%            plot(BlockInd{d}{e}{it+1}(:,2), it, 'or')
        end
    end
    [d,e,t,c]= getNextCell;
end

save([fmaux.data2dir '/BlockInd'], 'BlockInd');
%save([fmaux.data2dir '/OnTraj'], 'OnTraj');

