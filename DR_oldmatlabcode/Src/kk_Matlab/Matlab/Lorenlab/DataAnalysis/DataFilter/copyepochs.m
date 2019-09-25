function fout = copyepochs(fin)
% function fout = copyepochs(fin)
% A set of epochs are given in the output fields of fin.
% Generate fout from fin using these as the epochs and resetting
% the appropriate fields.


for an = 1:length(fin)
  fout(an).animal = fin(an).animal;
  fout(an).epochs = fin(an).output;
  for g = 1:length(fout(an).epochs)
    if isempty(fout(an).epochs{g})
      fout(an).data{g} = [];
      fout(an).excludetime{g} = [];
    else
      for e = 1:size(fout(an).epochs{g},1)
        fout(an).data{g}{e} = [];
        fout(an).excludetime{g}{e} = [];
      end
    end
  end
end

