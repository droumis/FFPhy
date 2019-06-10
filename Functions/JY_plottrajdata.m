function out = JY_plottrajdata(index, excludetimes, spikes, linpos, includestates, minV, varargin)
%
% out = plottrajdata(index, excludetimes, spikes, linpos, includestates, minV, varargin)
%
% Options:
%   'appendindex', 1 or 0 -- set to 1 to append the cell index to the
%   output [tetrode cell value].  Default 0.
% OUTPUT:
%   []

appendindex = 1;
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
           otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

%calculate state and lindist
[state, lindist] = getbehavestate(linpos, index(1,1), index(1,2), includestates, 'minlinvelocity', minV);  %filters linpos.statematrix.traj by includestates

%apply exclude times to state
exclind = isExcluded(linpos{index(1,1)}{index(1,2)}.statematrix.time,excludetimes);  %find excluded time indeces
state(find(exclind)) = -1; %set excluded times to -1

%calculate trajdata
trajdata{index(1,1)}{index(1,2)}{index(1,3)}{index(1,4)} = calclinfields(spikes,state,lindist,linpos, index);

%plot trajdata
for trj = 1:length(trajdata{index(1,1)}{index(1,2)}{index(1,3)}{index(1,4)})
    if ~isempty(trajdata{index(1,1)}{index(1,2)}{index(1,3)}{index(1,4)}{trj})
        figure(trj)
       % subplot(length(trajdata{index(1,1)}{index(1,2)}{index(1,3)}{index(1,4)}),1,trj)
        plot(trajdata{index(1,1)}{index(1,2)}{index(1,3)}{index(1,4)}{trj}(:,1), trajdata{index(1,1)}{index(1,2)}{index(1,3)}{index(1,4)}{trj}(:,5));
        xlim([0 trajdata{index(1,1)}{index(1,2)}{index(1,3)}{index(1,4)}{trj}(end,1)])
        xlabel('Linear Position')
        ylabel('Occ Normd Firing Rate')
        title(['Index ', num2str(index), ' Traj', num2str(trj) ])
    end
end

pause

%output
tmpout = [];

%appendindex if output is not empty
if appendindex && ~isempty(tmpout)
    tmpout = [repmat(index, size(tmpout,1), 1) tmpout];
end

out = tmpout; %
end