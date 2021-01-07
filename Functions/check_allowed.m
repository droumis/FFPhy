function check_allowed(allowed, passed)
varargs = passed(1:2:end);
if ~all(ismember(varargs, allowed))
    badies = varargs{~ismember(varargs, allowed)};
    error(sprintf('unexpected varargin: "%s" ', badies));
end