function out = funcSwitchBox(index, excludeperiods, spikes, linpos,varargin)

out = [];
for i = 1:length(varargin)
    switch varargin{i}
        case {'calcoccnormmeanrate','calcpeakrate', 'calcoutfieldfiring','calcTrackActive'}
            if (i == 1)
                out = feval(varargin{i},index, excludeperiods, spikes, linpos,'appendindex',1);
            else
                out = [out feval(varargin{i},index, excludeperiods, spikes, linpos)];
            end
        case {'calctotalmeanrate','getmeanrate','getspikewidth'}
            if (i == 1)
                out = feval(varargin{i},index, excludeperiods, spikes, 'appendindex',1);
            else
                out = [out feval(varargin{i},index, excludeperiods, spikes)];
            end
        
        otherwise
            error(['Nothing to do for called function: ', varargin{i}]);
  
    end           
end