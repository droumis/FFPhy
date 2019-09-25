function out = JY_getbarrier(animaldir,animalprefix,epochs, varargin)
%
% Produces a cell structure with the fields:
% time, barrier (linearized)
%   EPOCHS - N by 2 matrix, columns are [day epoch]
% barrier = 1 means the barrier is present
% barrier = 0 means the barrier is absent
%
%


loaddays = unique(epochs(:,1));
data = loaddatastruct(animaldir, animalprefix, 'data', loaddays);
for i = 1:size(epochs,1)
    timeslist=data{epochs(i,1)}{epochs(i,2)}.Pos.correcteddata(:,1);
    barriertime=[];
    if ~isempty(data{epochs(i,1)}{epochs(i,2)}.Events.Barrier)
    barriertime=[data{epochs(i,1)}{epochs(i,2)}.Events.Barrier(1,1)...
        data{epochs(i,1)}{epochs(i,2)}.Events.Barrier(1,1)+data{epochs(i,1)}{epochs(i,2)}.Events.Barrier(1,4)];
    end
    barrier=zeros(size(timeslist));
    if ~isempty(barriertime)
        barrier=isExcluded(timeslist,barriertime./10000);
        
    end
        
    out{epochs(i,1)}{epochs(i,2)}.barrier=barrier;
    out{epochs(i,1)}{epochs(i,2)}.time=timeslist;
end

end