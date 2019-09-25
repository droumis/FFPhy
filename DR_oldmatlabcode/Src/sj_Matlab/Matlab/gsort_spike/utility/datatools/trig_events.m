function events = trig_events(trig, trace, before, after);
%TRIG_EVENTS       Extracts events based on a trigger.
%   EVENTS = TRIG_EVENTS(TRIGGER, TRACE, BEFORE_SAMPLES, AFTER_SAMPLES);
%   Returns an N x (BEFORE+1+AFTER) matrix where each row contains the
%   values of TRACE lined up with respect to the (N) 0->~0 transitions of
%   TRIGGER; both TRIGGER and TRACE must be the same length.  BEFORE and
%   AFTER specify the number of samples to extract before and after the
%   time of the transition, respectively.
%
%   NOTE: When a trigger is within BEFORE samples of the start of the data
%   vector or within AFTER samples of the end, NaN values are used for the
%   invalid TRACE samples.

datalen = length(trace);
if (datalen ~= length(trig)), error('Trigger and data vector lengths do not match.');  end;

% Get 0->~0 transitions & make index matrix
marks = find(leading_edges(trig(:)));   
inds = repmat(marks, 1, before+after+1) + repmat([-before:after], length(marks), 1);

% Extract events
trace = double(trace);
trace(end+1) = NaN;                                     % store a dummy value here ...
inds([inds < 1] | [inds > datalen]) = datalen + 1;      %  ... and assign overflows to dummy pos
events = trace(inds);

% Screwy Matlab indexing: if inds is N x M, trace(inds) will be N x M;
%   Except if inds is 1 x M, in which case the dimensions of trace(inds)
%   are M x 1 if trace is a column vector.  So ...
if (size(events,2) == 1),  events = events';  end