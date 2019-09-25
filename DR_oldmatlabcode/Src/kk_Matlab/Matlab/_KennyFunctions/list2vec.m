function vector = list2vec(list,timevec)

%  vector = list2vec(list,timevec)
% converts list of periods (i.e. [starttimes , endtimes])
%       into vector of 0-1 for corresponding times in timevec

timevec = timevec(:);

% initialize
vector = zeros(size(timevec));

for p = 1:size(list,1)
    startind = lookup(list(p,1),timevec);
    endind = lookup(list(p,2),timevec);
    vector(startind:endind) = 1;
end

end