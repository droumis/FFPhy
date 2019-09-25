function recurrence = posrecurrence(smoothedpos);
%
% returns a sparse matrix which summarizes recurrence structure of position data
%

% for each position sample, find other position samples when rat's returns to
% within EPSILON of current position
EPSILON = 1;

numsamples = size(smoothedpos.data,1);

% how many recurrence events do we expect?
MAX_NUM_EVENTS = 0.01*numsamples^2;

% initialize the non-zero entries of the sparse matrix
row = numsamples*ones(MAX_NUM_EVENTS,1);
col = numsamples*ones(MAX_NUM_EVENTS,1);
val = zeros(MAX_NUM_EVENTS,1);

% keep a counter of the number of recurrence events we record
count = 0;
for i = 1:numsamples
    dist = hypot( ...
        smoothedpos.data(i,2) - smoothedpos.data(:,2), ...
        smoothedpos.data(i,3) - smoothedpos.data(:,3) );
    j = find(dist < EPSILON);
    j = j(j >= i);
    if ~isempty(j)
        newcount = count + numel(j);
        if newcount > MAX_NUM_EVENTS
            error('too many events! increase MAX_NUM_EVENTS');
        end
        row((count+1):newcount) = i*ones(size(j));
        col((count+1):newcount) = j;
        val((count+1):newcount) = ones(size(j));
        count = newcount;
    end
end

recurrence = sparse(row,col,val);

