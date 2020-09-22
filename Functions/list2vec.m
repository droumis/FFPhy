function vector = list2vec(list,timevec, varargin)


%  vector = list2vec(list,timevec)
% converts list of periods (i.e. [starttimes , endtimes])
%       into vector of 0-1 for corresponding times in timevec

%DR -  added a varargin 'marker' to input an array the length of list
%whose number at each row marker(x) gets input into the output array for all the
%vec indices between list(x, [start end]) times
% for example, turning a set of intervals into a vector that is labelled by
% the interval id num that each time was in

marker = ones(length(list(:,1)),1);
if (~isempty(varargin))
    assign(varargin{:});
end

timevec = timevec(:);

% initialize
vector = zeros(size(timevec));

for p = 1:size(list,1)
    startind = lookup(list(p,1),timevec);
    endind = lookup(list(p,2),timevec);
    vector(startind:endind) = marker(p);
end

end