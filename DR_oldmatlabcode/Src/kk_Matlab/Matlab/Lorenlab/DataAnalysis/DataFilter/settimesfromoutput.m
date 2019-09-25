function f = settimesfromoutput(f, fieldname, varargin)
% function f = settimesfromoutput(f, fieldname)

for an = 1:length(f)
  if isempty(f(an).animal)
    error(['You must define an animal for the filter before filtering the tetrodes'])
  end
  if isempty(f(an).epochs)
    error(['You must define the desired epochs for the filter before filtering the tetrodes'])
  end

  for g = 1:length(f(an).epochs)
    if isempty(f(an).epochs{g}) || isempty(f(an).output{g})
      f(an).times{g} = [];
      continue;
    end
    for e = 1:size(f(an).epochs{g},1)
      excludeperiods = f(an).excludetime{g}{e};

      if isempty(fieldname) % assume output is a Nx1 or 1xN vector
        goodinds = find(~isExcluded(f(an).output{g}{e}, excludeperiods));
        f(an).output{g}{e} = f(an).output{g}{e}(goodinds);
        f(an).times{g}{e} = f(an).output{g}{e};
      else
        if ~iscell(fieldname)
          goodinds = find(~isExcluded(f(an).output{g}(e).(fieldname), excludeperiods));
          allfields = fieldnames(f(an).output{g}(e));
          for ff = 1:length(allfields)
            f(an).output{g}(e).(allfields{ff}) = f(an).output{g}(e).(allfields{ff})(goodinds);
          end
          f(an).times{g}{e} = f(an).output{g}(e).(fieldname);
        else
          goodinds = ~isExcluded(f(an).output{g}(e).(fieldname{1}), excludeperiods);
          for ff = 2:length(fieldname)
            goodinds = goodinds & ~isExcluded(f(an).output{g}(e).(fieldname{ff}), excludeperiods);
          end

          allfields = fieldnames(f(an).output{g}(e));
          for ff = 1:length(allfields)
            f(an).output{g}(e).(allfields{ff}) = f(an).output{g}(e).(allfields{ff})(goodinds);
          end

          f(an).times{g}{e} = f(an).output{g}(e).(fieldname{1});
          for ff = 2:length(fieldname)
            f(an).times{g}{e} = cat(2,f(an).times{g}{e},f(an).output{g}(e).(fieldname{ff}));
          end
        end


      end
    end
  end
end
