function var = setcellstructfield(var, field, arg2, arg3, arg4)
% function var = setcellstructfields(var, field, epochs, values)
% function var = setcellstructfields(var, field, day, epochs, values)

if nargin==4
  day = arg2(:,1);
  epochs = arg2(:,2);
  values = arg3;
elseif nargin==5
  day = repmat(arg2,length(arg3),1);
  epochs = arg3;
  values = arg4;
else
  error('setfields: Arguments not understood (too many or too few)');
end

if length(values)==length(epochs) & iscell(values)
  for e = 1:length(epochs)
    var{day(e)}{epochs(e)}.(field) = values{e};
  end
else
  for e = 1:length(epochs)
    var{day(e)}{epochs(e)}.(field) = values;
  end
end

