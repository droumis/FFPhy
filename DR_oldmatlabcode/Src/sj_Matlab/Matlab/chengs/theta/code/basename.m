function bn= basename(path)
%function bn= basename(path)
% strip directories from pathname

i= strfind(path, '/');
bn= path(i(end)+1:end);
