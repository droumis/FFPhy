function out = dfa_template(idx, timeFilter, varargin)

% reward triggered swr and licks
% use with singledayanal
% DR_2020-10-13

fprintf('day %d \n',idx(1))
reqData = {'ca1rippleskons', 'pos','lick'};
check_required(reqData, varargin)

eventType = 'ca1rippleskons';
if ~isempty(varargin)
    assign(varargin{:});
end

out = init_out();
out.index = idx;
out.animal = animal;
day = idx(1);
eps = idx(2:end);




end

function out = init_out()
out.index = [];
out.animal = [];


end

