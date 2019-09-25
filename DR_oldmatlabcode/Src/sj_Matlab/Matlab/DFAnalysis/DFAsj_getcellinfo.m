function out = DFAsj_getcellinfo(index, excludetimes, cellinfo, varargin);
% out = DFAsj_getcellinfo(index, excludetimes, cellinfo, varargin);
% Returns cell properties from cellinfo field
%
% Options:
%   'appendindex', 1 or 0 -- set to 1 to append the cell index to the
%   output [tetrode cell value].  Default 0.
%


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


info = cellinfo{index(1)}{index(2)}{index(3)}{index(4)};

out.index = index;
out.meanrate=info.meanrate; % All values for epoch
out.spikewidth=info.spikewidth;
out.csi=info.csi;
out.propbursts=info.propbursts;
out.numspikes=info.numspikes;
out.tag=info.tag;


% if (appendindex)
%     out = [index rate]; %append the cell index to the value
% else
%     out = rate;
% end